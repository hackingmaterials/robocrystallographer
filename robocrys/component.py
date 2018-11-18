import numpy as np

from itertools import product

from pymatgen.core.structure import Structure
from pymatgen.core.surface import get_integer_index

from networkx import Graph, isolates, connected_components
from pymatgen.analysis.local_env import CrystalNN

cnn = CrystalNN()


class FragmentDescriber(object):

    def __init__(self, structure):
        self.structure = structure
        self.fragments = get_structure_fragments(structure)


def get_structure_fragments(structure):
    sites_to_translate = [i for i, site in enumerate(structure)
                          if 0 in site.frac_coords or 1 in site.frac_coords]
    structure.translate_sites(sites_to_translate, [1e-8, 1e-8, 1e-8],
                              frac_coords=False)

    # calculating nearest neighbours takes the most amount of time
    bonded_structure = cnn.get_bonded_structure(structure)
    supercell = bonded_structure * (3, 3, 3)

    # remove connectivity to outside supercell dimensions
    to_remove = [edge[:2] for edge in supercell.graph.edges(data=True)
                 if edge[2]['to_jimage'] != (0, 0, 0)]
    supercell.graph.remove_edges_from(to_remove)
    supercell.graph.remove_nodes_from(list(isolates(supercell.graph)))

    # find subgraphs
    supercell.graph = Graph(supercell.graph)
    all_subgraphs = (supercell.graph.subgraph(c)
                     for c in connected_components(supercell.graph))

    # extract fragments
    fragments = []
    for subgraph in all_subgraphs:
        fragment_coords = [supercell.structure[n].coords
                           for n in subgraph.nodes()]
        fragment_species = ['I'] * len(fragment_coords)

        fragment_structure = Structure(structure.lattice.matrix,
                                       fragment_species, fragment_coords,
                                       coords_are_cartesian=True)

        site_images = np.floor(fragment_structure.frac_coords).astype(int)
        site_images = tuple(map(tuple, site_images))

        # get the sites in the fragment that are in the unit cell only
        #  fragment_structure.translate_sites(
        #       range(fragment_structure.num_sites), [1, 1, 1])
        #  fragment_structure.merge_sites(tol=0.0001, mode='d')

        if (1, 1, 1) in site_images:
            dimensionality = _get_fragment_dimensionality(set(site_images))

            if dimensionality == 2 or dimensionality == 1:
                orientation = _get_fragment_orientation(
                    dimensionality, structure, supercell.structure,
                    fragment_coords)
            else:
                orientation = None

            fragments.append({'dimensionality': dimensionality,
                              'orientation': orientation,
                              'structure': fragment_structure})

    return fragments


def _get_fragment_dimensionality(site_images):
    reflections = _get_reflections(3)
    directions = np.array([np.array(x) - np.array(y) for x, y in reflections
                           if x in site_images and y in site_images])
    n_directions = len(directions)

    if n_directions < 2:
        return n_directions

    elif n_directions > 4:
        return 3

    # check if all directions are coplanar
    # by calculating scalar triple product
    if n_directions == 3:
        omega = np.linalg.det(directions.T)
    else:
        omega = np.linalg.det(directions[:3].T)
        omega += np.linalg.det(directions[1:].T)

    return 2 if omega == 0 else 3


def _get_fragment_orientation(dimensionality, structure, supercell,
                              fragment_coords, max_miller=3):
    supercell_fragment = Structure(supercell.lattice.matrix,
                                   ['X'] * len(fragment_coords),
                                   fragment_coords, coords_are_cartesian=True)
    n = max_miller * 2 + 1
    supercell_fragment.make_supercell((n, n, n))
    site_images = np.floor([structure.lattice.get_fractional_coords(c)
                            for c in supercell_fragment.cart_coords])

    site_images = site_images[::n**3]
    g = site_images.sum(axis=0) / site_images.shape[0]

    # run singular value decomposition
    u, s, vh = np.linalg.svd(site_images - g)

    # get direction (first column is best fit line, 3rd column is unitary norm)
    index = 2 if dimensionality == 2 else 0
    orientation = vh[index, :]

    return get_integer_index(orientation)


def _get_reflections(size):
    # size should always be odd
    # if size is 3
    # swap 2 and 0, leave 1 unchanged. E.g. (2, 1, 0) -> (0, 1, 2)
    # if size is 5
    # swap 4 and 0, swap 3 and 1, leave 2 unchanged. E.g. (4, 2, 1) -> (0, 2, 3)
    # effectively reflection symmetry through centre image
    reflections = set()
    n = size - 1
    for image in product(range(size), range(size), range(size)):
        if 0 or n in image:
            reflections.add(
                frozenset((image, tuple((-(x - n) for x in image)))))
    return reflections
