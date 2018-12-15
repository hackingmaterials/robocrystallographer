"""
This module provides a class for generating descriptions of condensed structure
data.

TODO: Handle distortion in connected polyhedra description.
"""

from collections import defaultdict
from typing import Dict, Any, Tuple, List

from robocrys.describe.adapter import DescriptionAdapter
from robocrys.util import (connected_geometries, geometry_to_polyhedra,
                           dimensionality_to_shape, get_el)
import inflect

en = inflect.engine()


class Describer(object):

    def __init__(self, distorted_tol: float = 0.6,
                 describe_mineral: bool = True,
                 describe_component_dimensionaly: bool = False,
                 describe_components: bool = False,
                 describe_oxidation_states: bool = True,
                 only_describe_cation_polyhedra_connectivity: bool = True):
        """

        Args:
            distorted_tol: The value under which the site geometry will be
                classified as distorted.
        """
        self.distorted_tol = distorted_tol
        self.describe_mineral = describe_mineral
        self.describe_component_dimensionality = describe_component_dimensionaly
        self.describe_components = describe_components
        self.describe_oxidation_state = describe_oxidation_states
        self.cation_polyhedra_only = \
            only_describe_cation_polyhedra_connectivity
        self._da = None

    def describe(self, condensed_structure) -> str:
        self._da = DescriptionAdapter(condensed_structure)

        description = list()

        if self.describe_mineral:
            description.append(self._get_mineral_description())

        if self.describe_component_dimensionality:
            description.append(self._get_component_dimensionality_description())

        if self.describe_components:
            description.append(self._get_component_descriptions())

        return " ".join(description)

    def _get_mineral_description(self) -> str:
        """Gets the mineral name and space group description.

        If the structure is a perfect match for a known prototype (e.g.
        the distance parameter is -1, the mineral name is the prototype name.
        If a structure is not a perfect match but similar to a known mineral,
        "-like" will be added to the mineral name. If the structure is a good
        match to a mineral but contains a different number of element types than
        the mineral prototype, "-derived" will be added to the mineral name.

        Returns:
            The description of the mineral name.
        """
        # TODO: indicate when molecules have been simplified
        if self._da.mineral['type']:
            if not self._da.mineral['n_species_type_match']:
                suffix = "-derived"
            elif self._da.mineral['distance'] >= 0:
                suffix = "-like"
            else:
                suffix = ""

            mineral_name = "{}{}".format(self._da.mineral['type'], suffix)

            desc = "{} is {} structured and".format(
                self._da.formula, mineral_name)

        else:
            desc = "{}".format(self._da.formula)

        desc += " crystallizes in the {} {} space group.".format(
            self._da.crystal_system, self._da.spg_symbol)
        return desc

    def _get_component_dimensionality_description(self):
        desc = "The structure is {}-dimensional".format(
            en.number_to_words(self._da.dimensionality))

        if len(self._da.components) == 1:
            desc += "."

        else:
            desc += " and consists of "

            dimen_component_descriptions = []
            for comp_dimen, formula_data in self._da.components.items():
                for formula, comp in formula_data.items():

                    count = en.number_to_words(comp['count'])

                    if 'molecule_name' in comp:
                        molecules = en.plural("molecule", comp['count'])
                        comp_desc = "{} {} {}".format(
                            count, comp['molecule_name'], molecules)

                    else:
                        shape = en.plural(dimensionality_to_shape[comp_dimen],
                                          count)
                        comp_desc = "{} {} {}".format(count, formula, shape)

                    if comp_dimen in [1, 2]:
                        orientations = set([c['orientation'] for c in
                                            comp['inequiv_components']])
                        dirs = en.plural("direction", len(orientations))
                        orientations = en.join([o for o in orientations])
                        comp_desc += " oriented in the {} {}".format(
                            orientations, dirs)

                    dimen_component_descriptions.append(comp_desc)

            desc += en.join(dimen_component_descriptions) + "."
        return desc

    def _get_component_descriptions(self, component_data: Dict[int, Any],
                                    n_components: int) -> str:
        if n_components == 1:
            # very messy way of getting the only component out of a nested
            # dict series where we don't know any of the keys in advance.
            comp = list(list(component_data.values())[0].values())[0][
                'inequiv_components'][0]

            if _component_contains_connected_polyhedra(
                    comp, self.cation_polyhedra_only):
                desc = "The structure contains "
            else:
                desc = ""

            desc += self._get_component_description(comp)

        else:
            dimen_component_descriptions = []
            for comp_dimen, formula_data in component_data.items():
                for formula, comp in formula_data.items():

                    # don't describe known molecules
                    if "molecule_name" in comp:
                        continue

                    if len(comp['inequiv_components']) == 1:
                        inequiv_comp = comp['inequiv_components'][0]
                        count = inequiv_comp['count']
                        prefix = "the" if count == 1 else "each"
                        shape = dimensionality_to_shape[comp_dimen]

                        if _component_contains_connected_polyhedra(
                                inequiv_comp, self.cation_polyhedra_only):
                            comp_desc = "{} {} {} contains ".format(
                                prefix.capitalize(), formula, shape)

                        else:
                            comp_desc = "In {} {} {}, ".format(
                                prefix, formula, shape)

                        comp_desc += self._get_component_description(
                            inequiv_comp)

                    else:
                        comp_desc = ""
                        for inequiv_comp in comp['inequiv_components']:
                            count = en.number_to_words(inequiv_comp['count'])
                            shape = en.plural(
                                dimensionality_to_shape[comp_dimen],
                                2)

                            if _component_contains_connected_polyhedra(
                                    inequiv_comp, self.cation_polyhedra_only):
                                comp_desc += "{} of the {} {}".format(
                                    count.capitalize(), formula, shape)

                            else:
                                comp_desc += "In {} of the {} {}, ".format(
                                    count, formula, shape)

                            comp_desc += self._get_component_description(
                                inequiv_comp)

                    dimen_component_descriptions.append(comp_desc)

            desc = " ".join(dimen_component_descriptions)
        return desc

    def _get_component_description(self, component: Dict[str, Any]) -> str:
        conn_sites, other_sites = _filter_connected_polyhedra(
            component, cation_polyhedra_only=self.cation_polyhedra_only)

        desc = []
        if conn_sites:
            desc.append(
                self._get_connected_polyhedral_sites_description(conn_sites))

        for site in other_sites:
            desc.append(self._get_site_description(
                site['element'], site['geometry'], nn_data=site['nn_data'],
                describe_bond_lengths=True))

        desc = " ".join(desc)
        return desc

    def _get_connected_polyhedral_sites_description(self,
                                                    connected_sites: List[dict]
                                                    ) -> str:
        """Gets a summary of the connected polyhedra sites."""

        # first reorganise the data slightly.
        # we only want to extract the connectivities for bonding to other
        # tet or oct sites. Also we want to separate the oct and tet sites so
        # that we can describe them separately.
        # futhermore, the tilt angles are only nicely defined for
        # corner sharing octahedra so we only report those here.

        connected_data = {}
        ordered_sites = []
        octahedral_tilts = []

        for geometry in connected_geometries:
            geometry_data = defaultdict(list)

            for site in connected_sites:
                if geometry in site['geometry']['type']:
                    connectivities = []
                    for el_data in site['nnn_data'].values():
                        for to_geometry in connected_geometries:
                            if to_geometry in el_data:
                                # store the different types of connectivity
                                connectivities.extend(
                                    el_data[to_geometry].keys())

                                # if corner sharing octahedra store the angles
                                if (geometry == 'octahedral' and
                                        'corner-sharing' in el_data[
                                            to_geometry]):
                                    octahedral_tilts.extend(
                                        el_data[to_geometry]['corner-sharing'][
                                            'angles'])

                    if connectivities:
                        connectivities = list(set(connectivities))
                        connectivities = [c.replace("sharing", "")
                                          for c in connectivities]

                        if len(connectivities) > 1:
                            connectivities_str = "a mixture of "
                        else:
                            connectivities_str = ""
                        connectivities_str += en.join(
                            connectivities) + "sharing"

                        if (site['polyhedra_formula'] not in
                                geometry_data[connectivities_str]):
                            geometry_data[connectivities_str].append(
                                site['polyhedra_formula'])
                            ordered_sites.append(site)

            if geometry_data:
                connected_data[geometry] = geometry_data

        desc = []

        for geometry in connected_geometries:
            if geometry in connected_data:
                for con_str, formulas in connected_data[geometry].items():
                    formulas_str = en.join(formulas)
                    polyhedra = geometry_to_polyhedra[geometry]
                    desc.append("{} {} {}".format(con_str, formulas_str,
                                                  polyhedra))

        desc[-1] = desc[-1] + "."
        all_sites_desc = [en.join(desc)]

        if octahedral_tilts:
            octahedral_tilts = [180 - t for t in octahedral_tilts]

            tilts = list(set(_rounded_bond_lengths(
                octahedral_tilts, decimal_places=1)))
            tilt_desc = "The corner-sharing octahedra have tilt angles "

            if len(tilts) == 1:
                tilt_desc += "of {}°."
            else:
                tilt_desc += "between {}–{}°.".format(min(tilts), max(tilts))

            all_sites_desc.append(tilt_desc)

        for site in ordered_sites:
            for bond_element, bond_data in site['nn_data'].items():
                all_sites_desc.append(self._get_bond_length_description(
                    site['element'], bond_element, bond_data))

        return " ".join(all_sites_desc)

    def _get_site_description(self, element: str, geometry: dict, nn_data: dict,
                              describe_bond_lengths: bool = True) -> str:
        """Gets a description of the geometry of a site.

        If the site likeness (order parameter) is less than ``distorted_tol``,
        "distorted" will be added to the geometry description.

        Args:
            element: The element of the species at the site.
            geometry: The site geometry as a :obj:`dict` with keys "geometry"
                and "likeness", corresponding to the geometry type (e.g.
                octahedral) and order parameter, respectively.
            nn_data: The nearest neighbour data for the site. This should have
                the same format as returned by
                :obj:`robocrys.site.SiteAnalyzer.get_nearest_neighbor_data`.
            distorted_tol: The value under which the site geometry will be
                classified as distorted.
            describe_bond_lengths: Whether to provide a description of the
                bond lengths. Defaults to ``True``.

        Returns:
            A description of the geometry and bonding of a site.
        """
        if not self.describe_oxidation_state:
            element = get_el(element)

        if geometry['likeness'] < self.distorted_tol:
            geometry_desc = "distorted"
        else:
            geometry_desc = ""
        geometry_desc += geometry['type']

        desc = "{} is bonded in {} geometry to ".format(
            element, en.a(geometry_desc))

        # tackle the case that the bonding is only to a single element
        if len(nn_data) == 1:
            return desc + self._get_single_element_bonding_description(
                element, nn_data, describe_bond_lengths)
        else:
            return desc + self._get_multi_element_bonding_description(
                element, nn_data, describe_bond_lengths)

    @staticmethod
    def _get_bond_length_description(element: str, bond_element: str,
                                     bond_data: dict) -> str:
        """Gets a description of the bonding of a site to an element.

        Bonding description based on the nearest neighbour data in
        `self.nearest_neighbour_data`. If you ask for the bonding description
        for site_index and element that are not bonded an error will be thrown.

        Args:
            element: The central element name (e.g. element of species at site
                to which bonds are made).
            bond_element (str): The element too which bonding will be described.
            bond_data: The bonding information as a :obj:`dict` with the
                format::

                    {
                        'n_sites': 6
                        'inequiv_groups': (
                            {
                                'n_sites': 4
                                'inequiv_id': 0
                                'dists': [1, 1, 2, 2]
                            },
                            {
                                'n_sites': 2
                                'inequiv_id': 1
                                'dists': [3, 3]
                            }
                        )
                    }

                This is the same as output by
                :obj:`robocrys.site.SiteAnalyzer.get_nearest_neighbor_data`.
                Alternatively, if the information for several sites has been
                merged together, the data can be formatted as::

                    {
                        'n_sites': 6
                        'dists': [1, 1, 1, 2, 2, 2, 2, 3, 3]
                        )
                    }

                Note that there are now more distances than there are number of
                sites. This is because n_sites gives the number of bonds to a
                specific site, whereas the distances are the complete set of
                distances for all similar (merged) sites.

        Returns:
            A description of the bond lengths.
        """
        element = get_el(element)
        bond_element = get_el(bond_element)

        if 'dists' in bond_data:
            dists = bond_data['dists']
        else:
            # combine all the inequivalent distances
            dists = sum([x['dists'] for x in bond_data['inequiv_groups']], [])

        # if only one bond length
        if len(dists) == 1:
            return "The {}–{} bond length is {}.".format(
                element, bond_element, _distance_to_string(dists[0]))

        discrete_bond_lengths = _rounded_bond_lengths(dists)

        # if multiple bond lengths but they are all the same
        if len(set(discrete_bond_lengths)) == 1:
            intro = "Both" if len(discrete_bond_lengths) == 2 else "All"
            return "{} {}–{} bond lengths are {}.".format(
                intro, element, bond_element, _distance_to_string(dists[0]))

        # if two sets of bond lengths
        if len(set(discrete_bond_lengths)) == 2:
            small = min(discrete_bond_lengths)
            small_count = en.number_to_words(discrete_bond_lengths.count(small))
            big = max(discrete_bond_lengths)
            big_count = en.number_to_words(discrete_bond_lengths.count(big))

            # length will only be singular when there is 1 small and 1 big len
            length = en.plural('length', int((small + big) / 2))

            return ("In this arrangement, there {} {} shorter ({}) and {} "
                    "longer ({}) {}–{} bond {}.").format(
                en.plural_verb('is', small_count), small_count,
                _distance_to_string(small), big_count,
                _distance_to_string(big), element, bond_element, length)

        # otherwise just detail the spread of bond lengths
        return ("There are a spread of {}–{} bond distances ranging from "
                "{}.").format(
            element, bond_element,
            _distance_range_to_string(min(discrete_bond_lengths),
                                      max(discrete_bond_lengths)))

    def _get_single_element_bonding_description(self, element, nn_data,
                                                describe_bond_lengths) -> str:
        """Utility func to get description of a site bonded to a one element."""
        desc = ""
        bond_element, bond_data = list(nn_data.items())[0]

        if not self.describe_oxidation_state:
            bond_element = get_el(bond_element)

        if bond_data['n_sites'] == 1:
            desc += "one {} atom. ".format(bond_element)

        else:
            if ('inequiv_groups' in bond_data and
                    len(bond_data['inequiv_groups']) == 1):
                equivalent = " equivalent "
            else:
                equivalent = " "

            desc += "{}{}{} atoms.".format(
                en.number_to_words(bond_data['n_sites']), equivalent,
                bond_element)

        if describe_bond_lengths:
            desc += " "
            desc += self._get_bond_length_description(element, bond_element,
                                                      bond_data)
        return desc

    def _get_multi_element_bonding_description(self, element, nn_data,
                                               describe_bond_lengths) -> str:
        """Utility func to get description of a site bonded to multi elems."""
        desc = ""

        bonding_atoms = []
        for el, data in nn_data.items():

            if not self.describe_oxidation_state:
                el = get_el(el)

            count = en.number_to_words(data['n_sites'])

            if ('inequiv_groups' in data and len(data['inequiv_groups']) == 1
                    and data['n_sites'] > 1):
                equiv = " equivalent "
            else:
                equiv = " "

            bonding_atoms.append("{}{}{}".format(count, equiv, el))

        desc += "{} atoms.".format(en.join(bonding_atoms))

        intro = None
        # only describe the equivalent groups if inequiv_groups key is in the
        # data. if it is not in the data this means the information for several
        # similar sites has been merged. See the merge_similar_sites method in
        # structure.py for more information.
        for i, (bond_element, bond_data) in enumerate(nn_data.items()):

            if not self.describe_oxidation_state:
                bond_element = get_el(bond_element)

            if ('inequiv_groups' in bond_data and
                    len(bond_data['inequiv_groups']) == 1 and
                    bond_data['n_sites'] > 1):
                # already described these previously
                pass

            elif 'inequiv_groups' in bond_data and bond_data['n_sites'] > 1:
                intro = "Of these, the" if not intro else "The"

                desc += ("{} {} atoms are found in {} distinct "
                         "environments.").format(
                    intro, bond_element,
                    en.number_to_words(len(bond_data['inequiv_groups'])))

            if describe_bond_lengths:
                desc += " "
                desc += self._get_bond_length_description(element, bond_element,
                                                          bond_data)
        return desc


def _component_contains_connected_polyhedra(component: Dict[str, Any],
                                            cation_polyhedra_only: bool = False
                                            ) -> bool:
    """Check if a component contains connected polyhedra."""
    return any(['polyhedra_formula' in s for s in component['sites']
                if not cation_polyhedra_only or '+' in s['element']])


def _filter_connected_polyhedra(component: Dict[str, Any],
                                cation_polyhedra_only: bool = False
                                ) -> Tuple[List, List]:
    """Function to separate out connected and non-connected sites."""
    conn = [s for s in component['sites'] if 'polyhedra_formula' in s and
            (not cation_polyhedra_only or '+' in s['element'])]
    other = [s for s in component['sites'] if
             'polyhedra_formula' not in s or
             (cation_polyhedra_only and '+' not in s['element'])]
    return conn, other


def _rounded_bond_lengths(data, decimal_places=3) -> Tuple[float]:
    """Utility function to round bond lengths to a number of decimal places."""
    return tuple(float("{:.{}f}".format(x, decimal_places)) for x in data)


def _distance_to_string(distance, decimal_places=2) -> str:
    """Utility function to round a distance and add an Angstrom symbol."""
    return "{:.{}f} Å".format(distance, decimal_places)


def _distance_range_to_string(dist_a, dist_b, decimal_places=2) -> str:
    """Utility function to format a range of distances."""
    return "{:.{}f}–{:.{}f} Å".format(dist_a, decimal_places, dist_b,
                                      decimal_places)
