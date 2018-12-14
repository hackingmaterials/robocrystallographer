"""
This module provides functions for extracting information about site geometry.

TODO: distortion of geometry e.g. elongated along an axis
TODO: maybe have custom descriptions for octahedron, tetrahedron etc
TODO: handle the case where no geometry type is given, just the CN
"""
import copy
from collections import defaultdict
from typing import Dict, Any, List, Union

import numpy as np

from pymatgen import Composition
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.core.periodic_table import get_el_sp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.coord import get_angle
from pymatgen.util.string import formula_double_format
from robocrys.fingerprint import get_site_fingerprints
from robocrys.util import connected_geometries, get_el, defaultdict_to_dict


class SiteAnalyzer(object):
    """Class to extract information on site geometry and bonding.

    Args:
        bonded_structure: A bonded structure with nearest neighbor data
            included. For example generated using
            :class:`pymatgen.analysis.local_env.CrystalNN` or
            :class:`pymatgen.analysis.local_env.VoronoiNN`.
        use_symmetry_equivalent_sites: Whether to use symmetry to determine if
            sites are inequivalent. If ``False``, the site geometry and (next)
            nearest neighbor information will be used.
        symprec: The tolerance used when determining the symmetry of
            the structure. The symmetry can used both to determine if multiple
            sites are symmetrically equivalent and to obtain the symmetry labels
            for each site.
        use_iupac_formula (bool, optional): Whether to order formulas
            by the iupac "electronegativity" series, defined in
            Table VI of "Nomenclature of Inorganic Chemistry (IUPAC
            Recommendations 2005)". This ordering effectively follows
            the groups and rows of the periodic table, except the
            Lanthanides, Actanides and hydrogen. If set to ``False``, the
            elements will be ordered according to the electronegativity values.
    """

    def __init__(self,
                 bonded_structure: StructureGraph,
                 use_symmetry_equivalent_sites: bool = False,
                 symprec: float = 0.01,
                 use_iupac_formula: bool = True):
        self.bonded_structure = bonded_structure
        self.site_fingerprints = get_site_fingerprints(
            bonded_structure.structure)

        sga = SpacegroupAnalyzer(bonded_structure.structure,
                                 symprec=symprec)
        sym_data = sga.get_symmetry_dataset()
        if use_symmetry_equivalent_sites:
            self.equivalent_sites = sym_data['equivalent_atoms']

        else:
            self.equivalent_sites = self._calculate_equivalent_sites()

        self.symmetry_labels = self._calculate_symmetry_labels(
            sym_data['equivalent_atoms'])

        self.use_iupac_formula = use_iupac_formula

    def get_site_summary(self, site_index: int) -> Dict[str, Any]:
        element = str(self.bonded_structure.structure[site_index].specie)
        geometry = self.get_site_geometry(site_index)

        nn_sites = self.get_nearest_neighbors(
            site_index, inc_inequiv_id=True)
        nnn_sites = self.get_next_nearest_neighbors(
            site_index, inc_inequiv_id=True)

        nn_ids = [nn_site['inequiv_id'] for nn_site in nn_sites]
        nnn_ids = [nnn_site['inequiv_id'] for nnn_site in nnn_sites]

        nnn_geometries = [nnn_site['geometry'] for nnn_site in nnn_sites]

        equiv_sites = [i for i in range(len(self.equivalent_sites))
                       if self.equivalent_sites[i] == self.equivalent_sites[
                           site_index]]
        sym_labels = tuple(set([self.symmetry_labels[x] for x in equiv_sites]))

        def order_elements(el):
            if self.use_iupac_formula:
                return [get_el_sp(el).X, el]
            else:
                return [get_el_sp(el).iupac_ordering, el]

        # check if site is connected polyhedra
        poly_formula = None
        if (geometry['type'] in connected_geometries and
                any([nnn_geometry['type'] in connected_geometries
                     for nnn_geometry in nnn_geometries])):
            nn_els = [get_el(nn_site['element']) for nn_site in nn_sites]
            comp = Composition("".join(nn_els))
            el_amt_dict = comp.get_el_amt_dict()

            poly_formula = ""
            for e in sorted(el_amt_dict.keys(), key=order_elements):
                poly_formula += e
                poly_formula += formula_double_format(el_amt_dict[e])

        return {'element': element, 'geometry': geometry, 'nn': nn_ids,
                'nnn': nnn_ids, 'poly_formula': poly_formula,
                'sym_labels': sym_labels}

    def get_bond_summary(self, site_index: int
                         ) -> Dict[int, Dict[int, List[float]]]:
        bonds = defaultdict(lambda: defaultdict(list))
        for nn_site in self.get_nearest_neighbors(site_index):
            to_site = nn_site['inequiv_id']
            bonds[site_index][to_site].append(nn_site['dist'])

        return defaultdict_to_dict(bonds)

    def get_connectivity_summary(self, site_index):
        connectivities = defaultdict(
            lambda: defaultdict(lambda: defaultdict(list)))

        for nnn_site in self.get_next_nearest_neighbors(
                site_index, inc_inequiv_id=True):
            to_site = nnn_site['inequiv_id']
            connectivity = nnn_site['connectivity']

            connectivities[site_index][to_site][connectivity].extend(
                nnn_site['angles'])

        return defaultdict_to_dict(connectivities)

    def get_site_geometry(self, site_index: int) -> Dict[str, Any]:
        """Gets the bonding geometry of a site.

        For example, "octahedral" or "square-planar".

        Args:
            site_index: The site index (zero based).

        Returns:
            The site geometry information as a :obj:`dict` with keys "type"
            and "likeness", corresponding to the geometry type (e.g. octahedral)
            and order parameter, respectively.
        """
        # get fingerprint as a list of tuples, e.g. [("op name", val), ...]
        site_fingerprint = list(self.site_fingerprints[site_index].items())

        # get coordination number with largest weight, ignore op names with
        # just the coordination number weight (e.g. containing "wt")
        parameter = max(site_fingerprint,
                        key=lambda x: x[1] if "wt" not in x[0] else 0)

        # return the geometry type without the CN at the end, e.g.
        # "square co-planar CN_4" -> "square co-planar"
        geometry = " ".join(parameter[0].split()[:-1])

        return {'type': geometry, 'likeness': parameter[1]}

    def get_nearest_neighbor_data(self, site_index: int,
                                  split_into_groups: bool = True
                                  ) -> Dict[str, Any]:
        """Gets a summary of all the nearest neighbors to a site.

        Args:
            site_index: The site index (zero based).
            split_into_groups: Whether to split the nearest neighbors into
                groups of inequivalent atoms.

        Returns:
            A summary of the nearest neighbor information as a dict. If
            ``split_into_groups=True``, the data is formatted as::

                {
                    'Sn': {
                        'n_sites': 6,
                        'inequiv_groups': [
                            {
                                'n_sites': 4,
                                'inequiv_id': 0,
                                'dists': [1, 1, 2, 2]
                            },
                            {
                                'n_sites': 2,
                                'inequiv_id': 1,
                                'dists': [3, 3]
                            }
                        ]
                    }
                }

            If ``split_into_groups=False``, the data is formatted as::

                {
                    'Sn': {
                        'n_sites': 6,
                        'dists': [1, 1, 2, 2, 3, 3]
                    }
                }

        """
        nn_info = self.get_nearest_neighbors(
            site_index, inc_inequiv_id=split_into_groups)

        if split_into_groups:
            # first group nearest neighbors by element and inequiv_id
            # e.g. grouped_nn looks like {'el': {'inequiv_id': [sites]}}
            grouped_nn = defaultdict(lambda: defaultdict(list))
            for site in nn_info:
                grouped_nn[site['element']][site['inequiv_id']].append(site)

            data = {}
            for element, sym_data in grouped_nn.items():
                n_sites = sum([len(sites) for sites in sym_data.values()])
                sym_groups = [
                    {'n_sites': len(sites),
                     'inequiv_id': sym_id,
                     'dists': [x['dist'] for x in sites]
                     } for sym_id, sites in sym_data.items()]
                data[element] = {'n_sites': n_sites,
                                 'inequiv_groups': sym_groups}

        else:
            # first group nearest neighbors by element
            # e.g. grouped_nn looks like {'el': [sites]}
            grouped_nn = defaultdict(list)
            for site in nn_info:
                grouped_nn[site['element']].append(site)

            data = {element: {'n_sites': len(sites),
                              'dists': [x['dist'] for x in sites]}
                    for element, sites in grouped_nn.items()}

        return data

    def get_next_nearest_neighbor_data(self, site_index: int
                                       ) -> Dict[str, Any]:
        """Gets a summary of the next nearest neighbor connectivity.

        Args:
            site_index: The site index (zero based).

        Returns:
            A summary of the next nearest neighbor information as a dict.
            Formatted as::

                {
                    'Sn': {
                        'octahedral': {
                            'corner-sharing': {
                                'n_sites': 8,
                                'angles': [180, 180, 180, ...]
                            }
                        }
                    }
                }

        """
        nnn_info = self.get_next_nearest_neighbors(site_index)

        # group next nearest neighbors by element, connectivity and geometry
        # e.g. grouped_nnn looks like {el: {connectivity: {geometry: [sites]}}}
        grouped_nnn = defaultdict(
            lambda: defaultdict(lambda: defaultdict(list)))

        for site in nnn_info:
            grouped_nnn[site['element']][
                site['geometry']['type']][site['connectivity']].append(site)

        nnn_data = {}
        for element, geom_data in grouped_nnn.items():
            nnn_el_data = {}
            for geometry, con_data in geom_data.items():
                nnn_geom_data = {}
                for connectivity, sites in con_data.items():
                    nnn_geom_data[connectivity] = {
                        'n_sites': len(sites),
                        'angles': [angle for site in sites
                                   for angle in site['angles']]}
                nnn_el_data[geometry] = nnn_geom_data
            nnn_data[element] = nnn_el_data
        return nnn_data

    def get_inequivalent_site_ids(self, site_ids: List[int]) -> List[int]:
        """Gets the inequivalent sites from a list of site indices.

        Args:
            site_ids: The site indices.

        Returns:
            The inequivalent site indices. For example, if a structure has 4
            sites where the first two are equivalent and the last two are
            inequivalent. If  ``site_ids=[0, 1, 2, 3]`` the output will be::

                [0, 2, 3]

        """
        return list(set(self.equivalent_sites[i] for i in site_ids))

    def get_nearest_neighbors(self, site_index: int,
                              inc_inequiv_id: bool = True
                              ) -> List[Dict[str, Any]]:
        """Gets information about the bonded nearest neighbors.

        Args:

            site_index: The site index (zero based).
            inc_inequiv_id: Whether to include the inequivalent site ids.

        Returns:
            For each site bonded to ``site_index``, returns a :obj:`dict`
            with the format::

                {'element': el, 'dist': distance}

            If ``inc_inequiv_id=True``, the data will have an additional key
            ``'inequiv_id'`` corresponding to the inequivalent site index. E.g.
            if two sites are structurally/symmetrically equivalent (depending
            on the value of ``self.use_symmetry`` then they will have the same
            ``inequiv_id``.
        """

        nn_sites = self.bonded_structure.get_connected_sites(site_index)

        if inc_inequiv_id:
            return [{'element': str(site.site.specie),
                     'inequiv_id': self.equivalent_sites[site.index],
                     'dist': site.dist} for site in nn_sites]
        else:
            return [{'element': str(site.site.specie),
                     'dist': site.dist} for site in nn_sites]

    def get_next_nearest_neighbors(self, site_index: int,
                                   inc_inequiv_id: bool = False
                                   ) -> List[Dict[str, Any]]:
        """Gets information about the bonded next nearest neighbors.

        Args:
            site_index: The site index (zero based).
            inc_inequiv_id: Whether to include the inequivalent site ids.

        Returns:
            A list of the next nearest neighbor information. For each next
            nearest neighbor site, returns a :obj:`dict` with the format::

                {'element': el, 'connectivity': con, 'geometry': geom,
                 'angles': angles}

            The ``connectivity`` property is the connectivity type to the
            next nearest neighbor, e.g. "face-sharing", "corner-sharing", or
            "edge-sharing". The ``geometry`` property gives the geometry of the
            next nearest neighbor site. See the ``get_site_geometry`` method for
            the format of this data. The ``angles`` property gives the bond
            angles between the site and the next nearest neighbour. Returned as
            a :obj:`list` of :obj:`int`. Multiple bond angles are given when
            the two sites share more than nearest neighbor (e.g. if they are
            face-sharing or edge-sharing).
        """

        def get_coords(a_site_index, a_site_image):
            return np.asarray(
                self.bonded_structure.structure.lattice.get_cartesian_coords(
                    self.bonded_structure.structure.frac_coords[a_site_index] +
                    a_site_image))

        nn_sites = self.bonded_structure.get_connected_sites(site_index)
        next_nn_sites = [site for nn_site in nn_sites for site in
                         self.bonded_structure.get_connected_sites(
                             nn_site.index, jimage=nn_site.jimage)]

        nn_sites_set = set((site.index, site.jimage) for site in nn_sites)

        seen_nnn_sites = set()
        next_nn_summary = []
        for nnn_site in next_nn_sites:
            if (nnn_site.index == site_index and nnn_site.jimage == (0, 0, 0)
                    or (nnn_site.index, nnn_site.jimage) in seen_nnn_sites):
                # skip the nnn site if it is the original atom of interest
                continue

            seen_nnn_sites.add((nnn_site.index, nnn_site.jimage))

            sites = set((site.index, site.jimage) for site in
                        self.bonded_structure.get_connected_sites(
                            nnn_site.index, jimage=nnn_site.jimage))
            shared_sites = nn_sites_set.intersection(sites)
            n_shared_atoms = len(shared_sites)

            if n_shared_atoms == 1:
                connectivity = 'corner-sharing'
            elif n_shared_atoms == 2:
                connectivity = 'edge-sharing'
            else:
                connectivity = 'face-sharing'

            site_coords = get_coords(site_index, (0, 0, 0))
            nnn_site_coords = get_coords(nnn_site.index, nnn_site.jimage)
            nn_site_coords = [get_coords(nn_site_index, nn_site_image)
                              for nn_site_index, nn_site_image in shared_sites]

            # can't just use Structure.get_angles to calculate angles as it
            # doesn't take into account the site image
            angles = [get_angle(site_coords - x, nnn_site_coords - x)
                      for x in nn_site_coords]

            geometry = self.get_site_geometry(nnn_site.index)

            summary = {'element': str(nnn_site.site.specie),
                       'connectivity': connectivity,
                       'geometry': geometry,
                       'angles': angles}

            if inc_inequiv_id:
                summary['inequiv_id'] = self.equivalent_sites[nnn_site.index]

            next_nn_summary.append(summary)

        return next_nn_summary

    def _calculate_equivalent_sites(self,
                                    geometry_likeness_tol: float = 0.1,
                                    bond_dist_tol: float = 0.1,
                                    bond_angle_tol: float = 0.1
                                    ) -> List[int]:
        """Determines the indices of the structurally inequivalent sites.

        Args:
            geometry_likeness_tol: The tolerance to use when determining whether
                two geometry likeness parameters are the same.
            bond_dist_tol: The tolerance used to determine if two bond lengths
                are the same.
            bond_angle_tol: The tolerance used to determine if two bond angles
                are the same.

        Two sites are considered equivalent if they are the same element,
        geometry, and (next) nearest neighbor summaries.

        Returns:
            A :obj:`list` of indices mapping each site in the structure to a
            structurally equivalent site. For example, if the first two sites
            are equivalent and the last two are both inequivalent, the data will
            be formatted as::

                [0, 0, 2, 3]

        """
        # TODO: Consider using site fingerprint rather than geometry type.
        inequiv_sites = {}
        equivalent_sites = []

        for site_id, site in enumerate(self.bonded_structure.structure):
            element = get_el_sp(site.specie)
            geometry = self.get_site_geometry(site_id)
            nn_sites = self.get_nearest_neighbors(
                site_id, inc_inequiv_id=False)
            nnn_sites = self.get_next_nearest_neighbors(
                site_id, inc_inequiv_id=False)

            matched = False
            for inequiv_id, inequiv_site in inequiv_sites.items():
                elem_match = element == inequiv_site['element']
                geom_match = geometries_match(
                    geometry, inequiv_site['geometry'],
                    likeness_tol=geometry_likeness_tol)
                nn_match = nn_summaries_match(
                    nn_sites, inequiv_site['nn_sites'],
                    bond_dist_tol=bond_dist_tol)
                nnn_match = nnn_summaries_match(
                    nnn_sites, inequiv_site['nnn_sites'],
                    bond_angle_tol=bond_angle_tol)

                if elem_match and geom_match and nn_match and nnn_match:
                    equivalent_sites.append(inequiv_id)
                    matched = True
                    break

            if not matched:
                # no matches therefore store original site id
                equivalent_sites.append(site_id)
                site_data = {'element': element,
                             'geometry': geometry,
                             'nn_sites': nn_sites,
                             'nnn_sites': nnn_sites}
                inequiv_sites[site_id] = site_data

        return equivalent_sites

    def _calculate_symmetry_labels(self, sym_equivalent_atoms: List[int]
                                   ) -> Dict[int, int]:
        """Calculates the symmetry labels for all sites in the structure.

        The symmetry labels number the sites in the structure. If two sites
        are symmetrically equivalent they share the same symmetry label. The
        numbering begins at 1 for each element in the structure.

        Args:
            sym_equivalent_atoms: A :obj:`list` of indices mapping each site in
                the structure to a symmetrically equivalent site. The data
                should be formatted as given by the ``equivalent_atoms`` key in
                :meth`SpacegroupAnalyzer.get_symmetry_dataset()`.

        Returns:
            A mapping between the site index and symmetry label for that site.
        """
        symmetry_labels = defaultdict(dict)
        for specie in self.bonded_structure.structure.species:
            el_indices = self.bonded_structure.structure.indices_from_symbol(
                get_el(specie))
            equiv_indices = [sym_equivalent_atoms[x] for x in el_indices]

            count = 1
            equiv_id_to_sym_label = {}
            for el_id, equiv_id in zip(el_indices, equiv_indices):
                if equiv_id in equiv_id_to_sym_label:
                    symmetry_labels[el_id] = equiv_id_to_sym_label[equiv_id]
                else:
                    equiv_id_to_sym_label[equiv_id] = count
                    symmetry_labels[el_id] = count
                    count += 1

        return dict(symmetry_labels)


def geometries_match(geometry_a: Dict[str, Any], geometry_b: Dict[str, Any],
                     likeness_tol: float = 0.01) -> bool:
    """Determine whether two site geometries match.

    Geometry data should be formatted the same as produced by
    :meth:`robocrys.site.SiteAnalyzer.get_site_geometry`.

    Args:
        geometry_a: The first set of geometry data.
        geometry_b: The second set of geometry data.
        likeness_tol: The tolerance to use when determining whether two
            geometry likeness parameters are the same.

    Returns:
        Whether the two geometries are the same.
    """
    # TODO: Handle case when geometry types are both None
    return (geometry_a['type'] == geometry_b['type'] and
            abs(geometry_a['likeness'] - geometry_b['likeness']) < likeness_tol)


def nn_summaries_match(nn_sites_a: List[Dict[str, Union[int, str]]],
                       nn_sites_b: List[Dict[str, Union[int, str]]],
                       bond_dist_tol: float = 0.1,
                       match_bond_dists: bool = True) -> bool:
    """Determine whether two sets of nearest neighbors match.

    Nearest neighbor data should be formatted the same as produced by
    :meth:`robocrys.site.SiteAnalyzer.get_nearest_neighbors`.

    Args:
        nn_sites_a: The first set of nearest neighbors.
        nn_sites_b: The second set of nearest neighbors.
        bond_dist_tol: The tolerance used to determine if two bond lengths
            are the same.
        match_bond_dists: Whether to consider bond distances when matching.

    Returns:
        Whether the two sets of nearest neighbors match.
    """

    def nn_sites_order(nn_site):
        return [nn_site['element'], nn_site['dist']]

    if len(nn_sites_a) != len(nn_sites_b):
        return False

    nn_sites_a = sorted(nn_sites_a, key=nn_sites_order)
    nn_sites_b = sorted(nn_sites_b, key=nn_sites_order)

    dists_match = [abs(site_a['dist'] - site_b['dist']) < bond_dist_tol
                   if match_bond_dists else True
                   for site_a, site_b in zip(nn_sites_a, nn_sites_b)]
    elements_match = [site_a['element'] == site_b['element']
                      for site_a, site_b in zip(nn_sites_a, nn_sites_b)]

    return all(d and e for d, e in zip(dists_match, elements_match))


def nnn_summaries_match(nnn_sites_a: List[Dict[str, Any]],
                        nnn_sites_b: List[Dict[str, Any]],
                        bond_angle_tol: float = 0.1,
                        match_bond_angles: bool = True):
    """Determine whether two sets of next nearest neighbors match.

    Next nearest neighbor data should be formatted the same as produced by
    :meth:`robocrys.site.SiteAnalyzer.get_next_nearest_neighbors`.

    Args:
        nnn_sites_a: The first set of next nearest neighbors.
        nnn_sites_b: The second set of next nearest neighbors.
        bond_angle_tol: The tolerance used to determine if two bond angles
            are the same.
        match_bond_angles: Whether to consider bond angles when matching.

    Returns:
        Whether the two sets of next nearest neighbors match.
    """

    def nnn_sites_order(nnn_site):
        return [nnn_site['element'], nnn_site['geometry']['type'],
                nnn_site['connectivity'], sorted(nnn_site['angles'])]

    if len(nnn_sites_a) != len(nnn_sites_b):
        return False

    nnn_sites_a = sorted(nnn_sites_a, key=nnn_sites_order)
    nnn_sites_b = sorted(nnn_sites_b, key=nnn_sites_order)

    elements_match = [site_a['element'] == site_b['element']
                      for site_a, site_b in zip(nnn_sites_a, nnn_sites_b)]
    cons_match = [site_a['connectivity'] == site_b['connectivity']
                  for site_a, site_b in zip(nnn_sites_a, nnn_sites_b)]
    geoms_match = [site_a['geometry']['type'] == site_b['geometry']['type']
                   for site_a, site_b in zip(nnn_sites_a, nnn_sites_b)]
    angles_match = [all([abs(a_a - a_b) < bond_angle_tol for a_a, a_b in
                         zip(sorted(site_a['angles']),
                             sorted(site_b['angles']))])
                    if match_bond_angles else True
                    for site_a, site_b in zip(nnn_sites_a, nnn_sites_b)]

    return all(e and c and g and a for e, c, g, a in
               zip(elements_match, cons_match, geoms_match, angles_match))


def merge_similar_sites(sites: List[Dict[str, Any]]):
    """Merges sites with the same properties except bond angles and distances.

    Args:
        sites: A list of sites. Each site is formatted as a :ob:`dict` with the
            keys:

            - ``'element'`` (``str``): The element of the site.
            - ``'geometry'`` (``dict``): The geometry, as output by
                :meth:`SiteAnalyzer.get_site_geometry`.
            - ``'nn_data'`` (``dict``): The nearest neighbor data, as output by
                :meth:`SiteAnalyzer.get_nearest_neighbor_data`.
            - ``'nnn_data'`` (``dict``): The next nearest neighbor data, as
                given by :meth:`SiteAnalyzer.get_next_nearest_neighbor_data`.

    Returns:
        A list of merged sites with the same format as above. Merged sites
        have a different ``nn_data`` format than unmerged sites. For example,
        ``nn_data`` in unmerged sites is formatted as::

            {
                'Sn': {
                    'n_sites': 6,
                    'inequiv_groups': [
                        {
                            'n_sites': 4,
                            'inequiv_id': 0,
                            'dists': [1, 1, 2, 2]
                        },
                        {
                            'n_sites': 2,
                            'inequiv_id': 1,
                            'dists': [3, 3]
                        }
                    ]
                }
            }

        Merged sites do not contain an ``inequiv_groups`` key and are instead
        formatted as::

            {
                'n_sites': 6
                'dists': [1, 1, 1, 2, 2, 2, 2, 3, 3]
                )
            }

        Note that there are now more distances than there are number of sites.
        This is because n_sites gives the number of bonds to a specific site,
        whereas the distances are for the complete set of distances for all
        similar (merged) sites. Similarly, merged next nearest neighbor
        data can contain more angles than number of sites, however, the
        general format of the ``nnn_data`` dict is unaltered.
    """
    sites = copy.deepcopy(sites)
    new_sites = []

    for site in sites:

        matched = False
        for new_site in new_sites:
            elem_match = site['element'] == new_site['element']
            geom_match = geometries_match(
                site['geometry'], new_site['geometry'], likeness_tol=1)
            nn_match = nn_summaries_match(
                site['nn_data'], new_site['nn_data'],
                match_bond_dists=False)
            nnn_match = nnn_summaries_match(
                site['nnn_data'], new_site['nnn_data'], match_bond_angles=False)

            if elem_match and geom_match and nn_match and nnn_match:
                new_site['nn_data'] = _merge_nn_data(site['nn_data'],
                                                     new_site['nn_data'])
                new_site['nnn_data'] = _merge_nnn_data(site['nnn_data'],
                                                       new_site['nnn_data'])
                matched = True
                break

        if not matched:
            # no matches therefore store original site id
            new_sites.append(site)

    return new_sites


def _merge_nn_data(site_nn_data, new_site_nn_data):
    """Utility function to merge nearest neighbor data.

    See the ``merge_similar_sites`` docstring for information on the format of
    the merged data.

    Note an error will be thrown if this function is called on two sites that do
    not have matching nearest neighbor summaries (ignoring bond distances).
    """

    for el in site_nn_data:
        site_dists = [dist for group in
                      site_nn_data[el]['inequiv_groups']
                      for dist in group['dists']]

        if 'inequiv_groups' in new_site_nn_data[el]:
            # remove inequiv_groups key and group all distances
            # together
            groups = new_site_nn_data[el].pop('inequiv_groups')
            dists = [dist for dist_set in groups
                     for dist in dist_set['dists']]
            new_site_nn_data[el]['dists'] = dists + site_dists
        else:
            new_site_nn_data[el]['dists'] += site_dists

    return new_site_nn_data


def _merge_nnn_data(site_nnn_data, new_site_nnn_data):
    """Utility function to merge next nearest neighbor data.

    See the ``merge_similar_sites`` docstring for information on the format of
    the merged data.

    Note an error will be thrown if this function is called on two sites that do
    not have matching next nearest neighbor summaries (ignoring bond angles).
    """
    for el in site_nnn_data:
        for geometry in site_nnn_data[el]:
            for connectivity in site_nnn_data[el][geometry]:
                new_site_nnn_data[el][geometry][connectivity]['angles'].extend(
                    site_nnn_data[el][geometry][connectivity]['angles'])

    return new_site_nnn_data
