"""
This module provides a class for generating descriptions of condensed structure
data.

TODO: Add option to ignore equivalent groups.
"""

from typing import Dict, Text, Any

import inflect

from robocrys.component import Component

en = inflect.engine()

dimensionality_to_shape = {
    3: 'framework', 2: 'sheet', 1: 'ribbon', 0: 'cluster'}


class Describer(object):

    def __init__(self, distorted_tol: float = 0.6,
                 describe_mineral: bool = True,
                 describe_component_dimensionaly: bool = True,
                 describe_components: bool = True):
        """

        Args:
            distorted_tol: The value under which the site geometry will be
                classified as distorted.
        """
        self.distored_tol = distorted_tol
        self.describe_mineral = describe_mineral
        self.describe_component_dimensionality = describe_component_dimensionaly
        self.describe_components = describe_components

    def describe(self, condensed_structure) -> str:
        description = list()

        if self.describe_mineral:
            mineral_desc = get_mineral_description(
                condensed_structure['mineral'], condensed_structure['formula'],
                condensed_structure['spg_symbol'],
                condensed_structure['crystal_system'])
            description.append(mineral_desc)

        if self.describe_component_dimensionality:
            dimen_desc = get_component_dimensionality_description(
                condensed_structure['dimensionality'],
                condensed_structure['components'],
                condensed_structure['n_components'])
            description.append(dimen_desc)

        if self.describe_components:
            comp_descs = get_component_descriptions(
                condensed_structure['components'],
                condensed_structure['n_components'])
            description.append(comp_descs)

        return " ".join(description)


def get_mineral_description(mineral_data: dict, formula: str,
                            spg_symbol: str, crystal_system: str) -> str:
    """Gets the mineral name and space group description.

    If the structure is a perfect match for a known prototype (e.g.
    the distance parameter is -1, the mineral name is the prototype name.
    If a structure is not a perfect match but similar to a known mineral,
    "-like" will be added to the mineral name. If the structure is a good
    match to a mineral but contains a different number of element types than
    the mineral prototype, "-derived" will be added to the mineral name.

    Args:
        mineral_data: The mineral information as a :obj:`dict` with the keys
            "mineral", "distance", "n_species_types_match", corresponding to the
            mineral name, the fingerprint distance between the prototype and
            known mineral, and whether the number of species types in the
            structure matches the number in the known prototype, respectively.
            If no mineral match, mineral_data will be ``None``.
        formula: The formula of the structure.
        spg_symbol: The space group symbol.
        crystal_system: The crystal system.

    Returns:
        The description of the mineral name.
    """
    # TODO: indicate when molecules have been simplified
    if mineral_data['type']:
        if not mineral_data['n_species_type_match']:
            suffix = "-derived"
        elif mineral_data['distance'] >= 0:
            suffix = "-like"
        else:
            suffix = ""

        mineral_name = "{}{}".format(mineral_data['type'], suffix)

        desc = "{} is {} structured and".format(formula, mineral_name)

    else:
        desc = "{}".format(formula)

    desc += " crystallizes in the {} {} space group.".format(
        crystal_system, spg_symbol)
    return desc


def get_component_dimensionality_description(dimensionality: int,
                                             component_data: Dict[int, Any],
                                             n_components: int):
    desc = "The structure is {}-dimensional".format(
        en.number_to_words(dimensionality))

    if n_components == 1:
        desc += "."

    else:
        desc += " and consists of "

        dimen_component_descriptions = []
        for comp_dimen, formula_data in component_data.items():
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


def get_component_descriptions(component_data: Dict[int, Any],
                               n_components: int):
    if n_components == 1:
        # very messy way of getting the only component out of a nested
        # dict series where we don't know any of the keys in advance.
        comp = list(list(component_data.values())[0].values())[0][
            'inequiv_components'][0]
        desc = get_component_description(comp)

    else:
        dimen_component_descriptions = []
        for comp_dimen, formula_data in component_data.items():
            for formula, comp in formula_data.items():

                # don't describe known molecules
                if "molecule_name" in comp:
                    continue

                if len(comp['inequiv_components']) == 1:
                    count = comp['inequiv_components'][0]['count']
                    prefix = "the" if count == 1 else "each"
                    shape = dimensionality_to_shape[comp_dimen]

                    comp_desc = "In {} {} {}, ".format(prefix, formula, shape)
                    comp_desc += get_component_description(
                        comp['inequiv_components'][0])

                else:
                    comp_desc = ""
                    for inequiv_comp in comp['inequiv_components']:
                        count = en.number_to_words(inequiv_comp['count'])
                        shape = en.plural(dimensionality_to_shape[comp_dimen],
                                          2)
                        comp_desc += "In {} of the {} {}, ".format(
                            count, formula, shape)
                        comp_desc += get_component_description(inequiv_comp)

                dimen_component_descriptions.append(comp_desc)

        desc = " ".join(dimen_component_descriptions)
    return desc


def get_component_description(component: Component):
    desc = []

    for site in component['sites']:
        desc.append(get_site_description(
            site['element'], site['geometry'], nn_data=site['nn_data'],
            describe_bond_lengths=True))

    desc = " ".join(desc)
    return desc


def get_site_description(element: str, geometry: dict, nn_data: dict,
                         distorted_tol: float = 0.6,
                         describe_bond_lengths: bool = True) -> Text:
    """Gets a description of the geometry of a site.

    If the site likeness (order parameter) is less than ``distorted_tol``,
    "distorted" will be added to the geometry description.

    Args:
        element: The element of the species at the site.
        geometry: The site geometry as a :obj:`dict` with keys "geometry" and
            "likeness", corresponding to the geometry type (e.g. octahedral)
            and order parameter, respectively.
        nn_data: The nearest neighbour data for the site. This should have the
            same format as returned by
            :obj:`robocrys.site.SiteAnalyzer.get_nearest_neighbor_summary`.
        distorted_tol: The value under which the site geometry will be
            classified as distorted.
        describe_bond_lengths: Whether to provide a description of the
            bond lengths. Defaults to ``True``.

    Returns:
        A description of the geometry and bonding of a site.
    """
    geometry_desc = '' if geometry['likeness'] > distorted_tol else 'distorted '
    geometry_desc += geometry['type']

    desc = "{} is bonded in {} geometry to ".format(
        element, en.a(geometry_desc))

    # handle tetrahedral and octahedral geometries separately

    # tackle the case that the bonding is only to a single element
    if len(nn_data) == 1:
        return desc + _get_single_element_bonding_description(
            element, nn_data, describe_bond_lengths)
    else:
        return desc + _get_multi_element_bonding_description(
            element, nn_data, describe_bond_lengths)


def get_bond_length_description(element: str, bond_element: str,
                                bond_data: dict) -> str:
    """Gets a description of the bonding of a site to an element.

    Bonding description based on the nearest neighbour data in
    `self.nearest_neighbour_data`. If you ask for the bonding description
    for site_index and element that are not bonded an error will be thrown.

    Args:
        element: The central element name (e.g. element of species at site to
            which bonds are made).
        bond_element (str): The element too which bonding will be described.
        bond_data: The bonding information as a :obj:`dict` with the format::

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
            :obj:`robocrys.site.SiteAnalyzer.get_nearest_neighbor_summary`.
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

        # length will only be singular when there is 1 small and 1 big length
        length = en.plural('length', int((small + big)/2))

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


def _get_single_element_bonding_description(element, nn_data,
                                            describe_bond_lengths):
    """Utility func to get description of a site bonded to a single element."""
    desc = ""
    bond_element, bond_data = list(nn_data.items())[0]

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
        desc += get_bond_length_description(element, bond_element,
                                            bond_data)
    return desc


def _get_multi_element_bonding_description(element, nn_data,
                                           describe_bond_lengths):
    """Utility func to get description of a site bonded to multiple elements."""
    desc = ""

    bonding_atoms = []
    for el, data in nn_data.items():
        count = en.number_to_words(data['n_sites'])

        if ('inequiv_groups' in data and len(data['inequiv_groups']) == 1
                and data['n_sites'] > 1):
            equiv = " equivalent "
        else:
            equiv = " "

        bonding_atoms.append("{}{}{}".format(count, equiv, el))

    desc += "{} atoms.".format(en.join(bonding_atoms))

    intro = None
    # only describe the equivalent groups if inequiv_groups key is in the data.
    # if it is not in the data this means the information for several
    # similar sites has been merged. See the _merge_similar_sites method in
    # structure.py for more information.
    for i, (bond_element, bond_data) in enumerate(nn_data.items()):

        if ('inequiv_groups' in bond_data and
                len(bond_data['inequiv_groups']) == 1 and
                bond_data['n_sites'] > 1):
            # already described these previously
            continue

        elif 'inequiv_groups' in bond_data and bond_data['n_sites'] > 1:
            intro = "Of these, the" if not intro else "The"

            desc += ("{} {} atoms are found in {} distinct "
                     "environments.").format(
                intro, bond_element,
                en.number_to_words(len(bond_data['inequiv_groups'])))

        if describe_bond_lengths:
            desc += " "
            desc += get_bond_length_description(element, bond_element,
                                                bond_data)
    return desc


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
