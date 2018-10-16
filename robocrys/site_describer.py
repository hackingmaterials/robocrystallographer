"""
This module provides functions for getting site geometry descriptions

TODO:

Make functions for
 - distortion of geometry e.g. elongated along an axis
 - neighbours e.g. bonded to six symmetrically equivalent oxygen atoms
 - connectivity e.g. edge-sharing
"""

import numpy as np

import inflect

from collections import defaultdict

from pymatgen.analysis.local_env import CrystalNN
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from robocrys.fingerprint import get_site_fingerprints

en = inflect.engine()


class SiteDescriber(object):

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
            data = [{'element': x.site.specie.name,
                     'sym_id': equivalent_sites[x.index],
                     'dist': x.dist} for x in con_sites]
            self.nearest_neighbour_data.append(data)

        self._processed_neighbour_data = {}

    def get_site_description(self, site_index, bond_lengths=True):
        """
        Sn is bonded in an octahedral geometry to six symmetrically equivalent
        O atoms. All Sn-O bonds are the same length
        or In this arragnemnt there are two longer and four shorter
        Sn-O bonds.

        Sn is bonded in an octahedral geometry to 3 Br, 2 Cl, and 1 I atoms.
        Of these, the Br atoms are all symmetrically equivalent.
        All Sn-Br bonds are the same length (4.93 A).
        The Cl atoms are symmetrically equivalent.
        Both Sn-Cl bonds are the same length (3.92 A).
        The Sn-I bond length is 1.82 A.
        """
        geometry = self.get_site_geometry(site_index)
        nn_info = self._get_processed_nearest_neighbour_info(site_index)
        element = self.structure[site_index].specie.name

        desc = "{} is bonded in {} geometry ".format(
            element, en.a(geometry))

        if len(nn_info) == 1:
            bond_element, bond_data = list(nn_info.items())[0]

            if bond_data['n_sites'] == 1:
                desc += "to one {} atom. ".format(bond_element)

            elif len(bond_data['sym_groups']) == 1:
                desc += "to {} symmetrically equivalent {} atoms. ".format(
                    en.number_to_words(bond_data['sym_groups'][0]['n_sites']),
                    bond_element)

            else:
                desc += "to {} {} atoms. ".format(
                    en.number_to_words(bond_data['n_sites']), bond_element)

            if bond_lengths:
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
            return "All {}–{} bond lengths are {}. ".format(
                element, bond_element, _distance_to_string(dists[0]))

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
        return ("There is a spread of {}–{} bond distances, ranging from"
                "{}–{}. ").format(
            element, bond_element,
            _distance_range_to_string(min(discrete_bond_lengths),
                                      max(discrete_bond_lengths)))

    def get_nearest_neighbour_info(self, site_index):
        """Gets information about the bonded nearest neighbours.

        Args:
            site_index (int): The site index (zero based).

        Returns:
            (list of dict): For each site bonded to `site_index`, returns a
            dictionary of::

                {'element': el, 'sym_id': i, 'dist': d}

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
                        'sym_groups': [
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
                        ]
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
            sym_groups = [{
                'n_sites': len(sites),
                'sym_id': sym_id,
                'dists': [x['dist'] for x in sites]
            } for sym_id, sites in sym_data.items()]
            data[element] = {'n_sites': n_sites, 'sym_groups': sym_groups}

        return data




def _rounded_bond_lengths(data, decimal_places=3):
    rounded_data = [float("{:.{}f}".format(x, decimal_places)) for x in data]
    return rounded_data


def _distance_to_string(distance, decimal_places=2):
    return "{:.{}f} Å".format(distance, decimal_places)


def _distance_range_to_string(dist_a, dist_b, decimal_places=2):
    return "{:.{}f}–{:.{}f} Å".format(dist_a, decimal_places, dist_b,
                                      decimal_places)
