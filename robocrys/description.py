import inflect

en = inflect.engine()

dimensionality_to_shape = {2: 'sheet', 1: 'ribbon', 0: 'cluster'}

class Describer(object):

    def __init__(self, distorted_tol: float=0.6):
        """

        Args:
            distorted_tol: The value under which the site geometry will be
                classified as distorted.
        """
        self.distored_tol = distorted_tol

    def describe(self, condensed_structure):
        description = list()

        try:
            description.append(get_mineral_description(
                condensed_structure['mineral'], condensed_structure['formula']))
        except ValueError as e:
            raise e

        dimen_desc = " The structure is {} dimensional".format(
            en.number_to_words(condensed_structure['dimensionality']))

        if len(condensed_structure['components']) == 1:
            dimen_desc += "."

        else:
            dimen_desc += " and consists of"

            dimen_component_descriptions = []
            for comp in condensed_structure['components']:
                count = en.number_to_words(comp['count'])
                shape = en.plural(
                    dimensionality_to_shape[comp['dimensionality']], count)
                formula = comp['formula']

                desc = " {} {} {}".format(
                    count, formula, shape)

                if comp['dimensionality'] in [1, 2]:
                    desc += " oriented in the {} direction".format(
                        comp['orientation'])

                dimen_component_descriptions.append(desc)

            dimen_desc += en.join(dimen_component_descriptions)
            dimen_desc += "."

        description.append(dimen_desc)

        #     desc += "In each {}, ".format(shape)
        # else:
        #     desc += ". "
        #
        # logging.info(desc)
        #
        # sga = SpacegroupAnalyzer(structure)
        # structure = sga.get_symmetrized_structure()
        #
        # site_describer = SiteAnalyzer(structure)
        # for i, list_sites in enumerate(structure.equivalent_indices):
        #     # very rough way of not overloading with information about bond lengths
        #     bond_lengths = i == len(structure.equivalent_indices) - 1
        #     logging.info(site_describer.get_site_description(
        #         list_sites[0], describe_bond_lengths=bond_lengths))

        return "".join(description)


def get_mineral_description(mineral_data: dict, formula: str) -> str:
    """Gets the mineral name description.

    If the structure is a perfect match for a known prototype (e.g.
    the distance parameter is 1, the mineral name is the prototype name.
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

    Returns:
        (str): The description of the mineral name.
    """

    if not mineral_data['type']:
        raise ValueError("No mineral name in mineral_data, cannot provide "
                         "description.")

    if not mineral_data['n_species_type_match']:
        suffix = "-derived"
    elif mineral_data['distance'] < 1:
        suffix = "-like"
    else:
        suffix = ""

    mineral_name = "{}{}".format(mineral_data['type'], suffix)

    desc = "{} is {} structured.".format(formula, mineral_name)
    return desc


def get_site_description(element: str, geometry: dict, nn_data: dict,
                         distorted_tol: float=0.6,
                         describe_bond_lengths: bool=True) -> str:
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
        (str): A description of the geometry and bonding of a site.
    """
    geometry = geometry['geometry']
    geometry += '' if geometry['likeness'] > distorted_tol else 'distorted'

    desc = "{} is bonded in {} geometry to ".format(
        element, en.a(geometry))

    # First tackle the case that the bonding is only to a single element
    if len(nn_data) == 1:
        bond_element, bond_data = list(nn_data.items())[0]

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
            desc += get_bond_length_description(element, bond_element,
                                                bond_data)
        return desc

    # tackle the case where the bonding is to multiple elements
    bonding_atoms = ["{} {}".format(en.number_to_words(data['n_sites']), el)
                     for el, data in nn_data.items()]
    desc += "{} atoms. ".format(en.join(bonding_atoms))

    intro = None
    for i, (bond_element, bond_data) in enumerate(nn_data.items()):

        if len(bond_data['sym_groups']) == 1 and bond_data['n_sites'] > 1:
            intro = "Of these, the" if not intro else "The"

            desc += "{} {} atoms are symmetrically equivalent. ".format(
                intro, bond_element)

        elif bond_data['n_sites'] > 1:
            intro = "Of these, the" if not intro else "The"

            desc += ("{} {} atoms are found in {} symmetry distinct "
                     "environments. ").format(
                intro, bond_element,
                en.number_to_words(len(bond_data['sym_groups'])))

        if describe_bond_lengths:
            desc += get_bond_length_description(element, bond_element,
                                                bond_data)

    return desc


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

            This is the same as output by
            :obj:`robocrys.site.SiteAnalyzer.get_nearest_neighbor_summary`.

    Returns:
        A description of the bond lengths.
    """

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
