"""
This module provides functions for generating site geometry descriptions.

TODO: make functions for
 - distortion of geometry e.g. elongated along an axis
 - connectivity e.g. edge-sharing
 - maybe have custom descriptions for octahedrons, tetrahedron etc
 - handle the case where no geometry type is given, just the CN
"""

import inflect

from collections import defaultdict

from pymatgen.analysis.local_env import CrystalNN
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from robocrys.fingerprint import get_site_fingerprints

en = inflect.engine()


class SiteDescriber(object):
    """Class to generate descriptions of site geometry and bonding.

    Args:
        structure (Structure): A pymatgen `Structure` object.
        bonded_structure (NearNeighbors, optional): A structure object with
            nearest neighbour data included. One should use a subclasses of
            `NearNeighbors` such as `CrystalNN` or `VoronoiNN`. If no
            bonded structure provided, the bonding will be calculated using
            the `CrystalNN` class.
        symprec (float): The tolerance used when determining the symmetry of
            the structure. The symmetry is used to determine if multiple
            sites are symmetrically equivalent.
    """

    def __init__(self, structure, bonded_structure=None, symprec=0.01):
        self.structure = structure
        self.site_fingerprints = get_site_fingerprints(structure)

        if not bonded_structure:
            cnn = CrystalNN()
            bonded_structure = cnn.get_bonded_structure(structure)

        sga = SpacegroupAnalyzer(structure, symprec=symprec)
        equivalent_sites = sga.get_symmetry_dataset()['equivalent_atoms']

        self.nearest_neighbour_data = []
        for site_index in range(len(structure)):
            con_sites = bonded_structure.get_connected_sites(site_index)
            data = tuple({'element': x.site.specie.name,
                          'sym_id': equivalent_sites[x.index],
                          'dist': x.dist} for x in con_sites)
            self.nearest_neighbour_data.append(data)

        self._processed_neighbour_data = {}

    def get_site_description(self, site_index, describe_bond_lengths=True):
        """Gets a description of the geometry and bonding of a site.

        Args:
            site_index (int): The site index (zero based).
            describe_bond_lengths (bool, optional): Whether to provide a
                description of the bond lengths. Default to `True`.

        Returns:
            (str): A description of the geometry and bonding of a site.
        """
        geometry = self.get_site_geometry(site_index)
        nn_info = self._get_processed_nearest_neighbour_info(site_index)
        element = self.structure[site_index].specie.name

        desc = "{} is bonded in {} geometry to ".format(
                element, en.a(geometry))

        # First tackle the case that the bonding is only to a single element
        if len(nn_info) == 1:
            bond_element, bond_data = list(nn_info.items())[0]

            if bond_data['n_sites'] == 1:
                desc += "one {} atom. ".format(bond_element)

            elif len(bond_data['sym_groups']) == 1:
                desc += "{} symmetrically equivalent {} atoms. ".format(
                    en.number_to_words(bond_data['sym_groups'][0]['n_sites']),
                    bond_element)

            else:
                desc += "{} {} atoms. ".format(
                    en.number_to_words(bond_data['n_sites']), bond_element)

            if describe_bond_lengths:
                desc += self.get_bond_length_description(
                    site_index, bond_element)

            return desc

        # tackle the case where the bonding is to multiple elements
        bonding_atoms = ["{} {}".format(en.number_to_words(data['n_sites']), el)
                         for el, data in nn_info.items()]
        desc += "{} atoms. ".format(en.join(bonding_atoms))

        for i, (bond_element, bond_data) in enumerate(nn_info.items()):

            desc += "Of these, the " if i == 0 else "The "

            if len(bond_data['sym_groups']) == 1 and bond_data['n_sites'] > 1:
                desc += "{} atoms are symmetrically equivalent. ".format(
                    bond_element)

            elif bond_data['n_site'] > 1:
                desc += ("{} atoms are found in {} symmetry distinct "
                         "environments. ").format(
                    bond_element,
                    en.number_to_words(len(bond_data['sym_groups'])))

            if describe_bond_lengths:
                desc += self.get_bond_length_description(
                    site_index, bond_element)

        return desc

    def get_site_geometry(self, site_index, distorted_tol=0.6):
        """Gets the bonding geometry of a site.

        For example, "octahedral" or "square-planar". If the order parameter
        is less than  `distorted_tol`, "distorted" will be added to the
        geometry description.

        Args:
            site_index (int): The site index (zero based).
            distorted_tol (float): The value under which the site geometry will
            be classified as distorted.

        Returns:
            (str): The site bonding geometry as a string.
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

        distorted = '' if parameter[1] > distorted_tol else 'distorted '
        return distorted + geometry

    def get_bond_length_description(self, site_index, bond_element):
        """Gets a description of the bonding of a site to an element.

        Bonding description based on the nearest neighbour data in
        `self.nearest_neighbour_data`. If you ask for the bonding description
        for site_index and element that are not bonded an error will be thrown.

        Args:
            site_index (int): The site index (zero based).
            bond_element (str): The element too which bonding will be described.

        Returns:
            (str): A description of the bond lengths between `site_index` and
            any `bond_element` atoms.
        """
        bond_data = self._get_processed_nearest_neighbour_info(
            site_index)[bond_element]
        element = self.structure[site_index].specie.name

        dists = sum([x['dists'] for x in bond_data['sym_groups']], [])

        # if only one bond length
        if len(dists) == 1:
            return "The {}–{} bond length is {}. ".format(
                element, bond_element, _distance_to_string(dists[0]))

        discrete_bond_lengths = _rounded_bond_lengths(dists)

        # if multiple bond lengths but they are all the same
        if len(set(discrete_bond_lengths)) == 1:
            intro = "Both" if len(discrete_bond_lengths) == 2 else "All"
            return "{} {}–{} bond lengths are {}. ".format(
                intro, element, bond_element, _distance_to_string(dists[0]))

        # if two sets of bond lengths
        if len(set(discrete_bond_lengths)) == 2:
            small = min(discrete_bond_lengths)
            small_count = en.number_to_words(discrete_bond_lengths.count(small))
            big = max(discrete_bond_lengths)
            big_count = en.number_to_words(discrete_bond_lengths.count(big))

            return ("In this arrangement, there {} {} shorter ({}) and {} "
                    "longer ({}) {}–{} bond lengths. ").format(
                en.plural_verb('is', small_count), small_count,
                _distance_to_string(small), big_count,
                _distance_to_string(big), element, bond_element)

        # otherwise just detail the spread of bond lengths
        return ("There is a spread of {}–{} bond distances, ranging from "
                "{}. ").format(
            element, bond_element,
            _distance_range_to_string(min(discrete_bond_lengths),
                                      max(discrete_bond_lengths)))

    def get_nearest_neighbour_info(self, site_index):
        """Gets information about the bonded nearest neighbours.

        Args:
            site_index (int): The site index (zero based).

        Returns:
            (tuple of dict): For each site bonded to `site_index`, returns a
            dictionary of::

                ({'element': el, 'sym_id': i, 'dist': d}, ...)

            The `sym_id` property is the symmetry index for the site. E.g. if
            two sites are symmetrically  equivalent then they will have the
            same `sym_id`.
        """
        return self.nearest_neighbour_data[site_index]

    def _get_processed_nearest_neighbour_info(self, site_index):
        """Utility method for repacking the nearest neighbour information.

        Args:
            site_index (int): The site index (zero based).

        Returns:
            (dict): A summary of the nearest neighbour information as a dict.
            Formatted as::

                {
                    'Sn': {
                        'n_sites': 6
                        'sym_groups': (
                            {
                                'n_sites': 4
                                'sym_id': 0
                                'dists': [1, 1, 1, 2, 2, 2]
                            },
                            {
                                'n_sites': 2
                                'sym_id': 1
                                'dists': [3, 3]
                            }
                        )
                    }
                }
        """
        if site_index in self._processed_neighbour_data:
            return self._processed_neighbour_data[site_index]

        nn_info = self.get_nearest_neighbour_info(site_index)

        # first group nearest neighbours by element and sym_id
        # e.g. grouped_nn looks like {'el': {'sym_id': [sites]}}
        grouped_nn = defaultdict(lambda: defaultdict(list))
        for site in nn_info:
            grouped_nn[site['element']][site['sym_id']].append(site)

        data = {}
        for element, sym_data in grouped_nn.items():
            n_sites = sum([len(sites) for sites in sym_data.values()])
            sym_groups = tuple({
                'n_sites': len(sites),
                'sym_id': sym_id,
                'dists': [x['dist'] for x in sites]
            } for sym_id, sites in sym_data.items())
            data[element] = {'n_sites': n_sites, 'sym_groups': sym_groups}

        return data


def _rounded_bond_lengths(data, decimal_places=3):
    """Utility function to round bond lengths to a number of decimal places."""
    return tuple(float("{:.{}f}".format(x, decimal_places)) for x in data)


def _distance_to_string(distance, decimal_places=2):
    """Utility function to round a distance and add an Angstrom symbol."""
    return "{:.{}f} Å".format(distance, decimal_places)


def _distance_range_to_string(dist_a, dist_b, decimal_places=2):
    """Utility function to format a range of distances."""
    return "{:.{}f}–{:.{}f} Å".format(dist_a, decimal_places, dist_b,
                                      decimal_places)
