"""
This module provides a class for generating descriptions of condensed structure
data.

TODO:
    * Indicate when molecules have been simplified in the mineral description.
    * Handle distortion in connected polyhedra description.
"""

from typing import Dict, Any, Tuple, List

from pymatgen.core.periodic_table import get_el_sp, Specie
from robocrys.describe.adapter import DescriptionAdapter
from robocrys.util import (connected_geometries, geometry_to_polyhedra,
                           dimensionality_to_shape, get_el)
import inflect

en = inflect.engine()


class Describer(object):

    def __init__(self,
                 describe_mineral: bool = True,
                 describe_component_makeup: bool = True,
                 describe_components: bool = False,
                 describe_symmetry_labels: bool = True,
                 describe_oxidation_states: bool = True,
                 describe_bond_lengths: bool = True,
                 bond_length_decimal_places: int = 2,
                 distorted_tol: float = 0.6,
                 only_describe_cation_polyhedra_connectivity: bool = True,
                 latexify: bool = False):
        """A class to convert condensed structure data into text descriptions.

        Args:
            describe_mineral: Whether to describe the crystal mineral data.
            describe_component_makeup: Whether to describe the component makeup
                of the structure.
            describe_components: Whether to describe the component's sites.
            describe_symmetry_labels: Whether to include the symmetry labels
                when describing atomic sites.
            describe_oxidation_states: Whether to include oxidation states
                in the description.
            describe_bond_lengths: Whether to describe bond lengths.
            bond_length_decimal_places: Number of decimal places to round
                bond lengths.
            distorted_tol: The value under which the site geometry will be
                classified as distorted.
            only_describe_cation_polyhedra_connectivity: Whether to only
                describe cation polyhedra instead of both cation and anion
                polyhedra.
            latexify: Whether to latexify the description.
        """
        self.distorted_tol = distorted_tol
        self.describe_mineral = describe_mineral
        self.describe_component_dimensionality = describe_component_makeup
        self.describe_components = describe_components
        self.describe_symmetry_labels = describe_symmetry_labels
        self.describe_oxidation_state = describe_oxidation_states
        self.describe_bond_lengths = describe_bond_lengths
        self.bond_length_decimal_places = bond_length_decimal_places
        self.cation_polyhedra_only = \
            only_describe_cation_polyhedra_connectivity
        self.latexify = latexify
        self._da: DescriptionAdapter = None

    def describe(self, condensed_structure: Dict[str, Any]) -> str:
        """Convert a condensed structure into a text description.

        Args:
            condensed_structure: The condensed structure data, formatted as
                produced by :meth:`StructureCondenser.condense_structure`.

        Returns:
            A description of the structure.
        """
        self._da = DescriptionAdapter(condensed_structure)

        description = list()

        if self.describe_mineral:
            description.append(self._get_mineral_description())

        if self.describe_component_dimensionality:
            description.append(self._get_component_makeup_summary())

        if self.describe_components:
            description.append(self._get_all_component_descriptions())

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

    def _get_component_makeup_summary(self) -> str:
        """Gets a summary of the makeup of components in a structure.

        Returns:
            A description of the number of components and their dimensionalities
            and orientations.
        """
        desc = "The structure is {}-dimensional".format(
            en.number_to_words(self._da.dimensionality))

        if len(self._da.components) == 1:
            desc += "."

        else:
            desc += " and consists of "

            component_makeup_summaries = []
            for component_group in self._da.get_component_details(
                    group_components=True):
                count = en.number_to_words(component_group.count)
                dimensionality = component_group.dimensionality

                if component_group.molecule_name:
                    shape = en.plural("molecule", count)
                    formula = component_group.molecule_name
                else:
                    shape = en.plural(dimensionality_to_shape[dimensionality],
                                      count)
                    formula = component_group.formula

                comp_desc = "{} {} {}".format(count, formula, shape)

                if component_group.dimensionality in [1, 2]:
                    orientations = set(c.orientation for c in
                                       component_group.components)
                    direction = en.plural("direction", len(orientations))
                    comp_desc += " oriented in the {} {}".format(
                        en.join(orientations), direction)

                component_makeup_summaries.append(comp_desc)

            desc += en.join(component_makeup_summaries) + "."
        return desc

    def _get_all_component_descriptions(self) -> str:
        """Gets the descriptions of all components in the structure.

        Returns:
            A description of all components in the structure.
        """
        if len(self._da.components) == 1:
            return self._get_component_description(
                self._da.get_component_details()[0])

        else:
            component_groups = self._da.get_component_details(
                group_components=True)

            component_descriptions = []
            for group in component_groups:
                for component in group.components:

                    if group.molecule_name:
                        # don't describe known molecules
                        continue

                    formula = group.formula
                    group_count = group.count
                    component_count = component.count
                    shape = dimensionality_to_shape[group.dimensionality]

                    if group_count == component_count:
                        filler = "the" if group_count == 1 else "each"
                    else:
                        filler = "{} of the".format(component_count)

                    desc = "In {} {} {} ".format(filler, formula, shape)
                    desc += self._get_component_description(component.index)

                    component_descriptions.append(desc)

            return " ".join(component_descriptions)

    def _get_component_description(self, component_index: int) -> str:
        """Gets the descriptions of all sites in a component.

        Args:
            component_index: The index of the component

        Returns:
            The description for all sites in the structure.
        """
        site_descriptions = []
        for site_group in self._da.get_site_details(
                component_index, group_sites=True):

            if len(site_group.groups) == 1:
                site_descriptions.append(
                    self._get_site_description(site_groups.groups[0]))

        return " ".join(site_descriptions)

    def _get_site_description(self, site_index: int) -> str:
        """Gets a description of the geometry and bonding of a site.

        If the site likeness (order parameter) is less than ``distorted_tol``,
        "distorted" will be added to the geometry description.

        Args:
            site_index: An inequivalent site index.

        Returns:
            A description of the geometry and bonding of a site.
        """
        site = self._da.sites[site_index]
        element = _get_formatted_el(
            site['elememt'], self._da.sym_labels[site_index],
            use_oxi_state=self.describe_oxidation_state,
            use_sym_label=self.describe_symmetry_labels,
            latexify=self.latexify)

        if site['geometry']['likeness'] < self.distorted_tol:
            geometry_desc = "distorted"
        else:
            geometry_desc = ""
        geometry_desc += site['geometry']['type']

        desc = "{} is bonded in {} geometry to ".format(
            element, en.a(geometry_desc))

        return desc + self._get_nearest_neighbor_description(site_index)

    def _get_nearest_neighbor_description(self, site_index: int) -> str:
        """Gets a description of a sites nearest neighbors.

        Note: This function is intended to be run directly after
        :meth:`get_site_description`, as the output will not form a complete
        sentence on its own.

        Args:
            site_index: An inequivalent site index.

        Returns:
            A description of the nearest neighbors.
        """
        nn_details = self._da.get_nearest_neighbor_details(
            site_index, group_by_element=self.describe_symmetry_labels)

        last_count = 0
        nn_descriptions = []
        bond_length_descriptions = []
        for nn_site in nn_details:
            element = _get_formatted_el(
                nn_site.element, self._da.sym_labels[site_index],
                use_oxi_state=self.describe_oxidation_state,
                use_sym_label=self.describe_symmetry_labels,
                latexify=self.latexify)

            if len(nn_site.sites) == 1 and not self.describe_symmetry_labels:
                equivalent = "equivalent"
            else:
                equivalent = ""

            nn_descriptions.append("{} {} {}".format(
                en.number_to_words(nn_site.count), equivalent, element))
            bond_length_descriptions.append(
                self._get_bond_length_description(site_index, nn_site.sites))

            last_count = nn_site.count

        s_atoms = "atom" if last_count == 1 else "atoms"
        nn_descriptions = "{} {}.".format(en.join(nn_descriptions), s_atoms)

        return " ".join(nn_descriptions) + " ".join(bond_length_descriptions)

    def _get_bond_length_description(self, from_site: int,
                                     to_sites: List[int]) -> str:
        """Gets a description of the bond lengths between sites.

        Args:
            from_site: An inequivalent site index.
            to_sites: A :obj:`list` of site indices. The site indices should
                all be for the same element.

        Returns:
            A description of the bond lengths.
        """
        from_element = _get_formatted_el(
            self._da.elements[from_site],
            self._da.sym_labels[from_site],
            use_oxi_state=False,
            use_sym_label=self.describe_symmetry_labels)
        to_element = _get_formatted_el(
            self._da.elements[to_sites[0]],
            self._da.sym_labels[to_sites[0]],
            use_oxi_state=False,
            use_sym_label=self.describe_symmetry_labels)

        dists = self._da.get_distance_details(from_site, to_sites)

        # if only one bond length
        if len(dists) == 1:
            return "The {}–{} bond length is {}.".format(
                from_element, to_element, _distance_to_string(dists[0]))

        discrete_bond_lengths = _rounded_bond_lengths(dists)

        # if multiple bond lengths but they are all the same
        if len(set(discrete_bond_lengths)) == 1:
            intro = "Both" if len(discrete_bond_lengths) == 2 else "All"
            return "{} {}–{} bond lengths are {}.".format(
                intro, from_element, to_element, _distance_to_string(dists[0]))

        # if two sets of bond lengths
        if len(set(discrete_bond_lengths)) == 2:
            small = min(discrete_bond_lengths)
            small_count = en.number_to_words(discrete_bond_lengths.count(small))
            big = max(discrete_bond_lengths)
            big_count = en.number_to_words(discrete_bond_lengths.count(big))

            # length will only be singular when there is 1 small and 1 big len
            length = en.plural('length', int((small + big) / 2))

            return ("There {} {} shorter ({}) and {} "
                    "longer ({}) {}–{} bond {}.").format(
                en.plural_verb('is', small_count), small_count,
                _distance_to_string(small), big_count,
                _distance_to_string(big), from_element, to_element, length)

        # otherwise just detail the spread of bond lengths
        return ("There are a spread of {}–{} bond distances ranging from "
                "{}.").format(
            from_element, to_element,
            _distance_range_to_string(min(discrete_bond_lengths),
                                      max(discrete_bond_lengths)))


def _get_formatted_el(element: str,
                      sym_label: str,
                      use_oxi_state: bool = True,
                      use_sym_label: bool = True,
                      latexify: bool = False):
    """Formats the element string.

    Performs a variety of functions, including:

    - Changing "Sn+0" to "Sn".
    - Inserting the symmetry label between the element and oxidation state, if
        required.
    - Removing the oxidation state if required.
    - Latexifying the element and oxidation state.

    Args:
        element: The element string (possibly including the oxidation state.
            E.g. "Sn" or "Sn+2".
        sym_label: The symmetry label. E.g. "(1)"
        use_oxi_state: Whether to include the oxidation state, if present.
        use_sym_label: Whether to use the symmetry label.
        latexify: Whether to convert the str for use in latex.

    Returns:
        The formatted element string.
    """
    specie = get_el_sp(element)

    if isinstance(specie, Specie):
        oxi_state = specie.oxi_state
        if oxi_state == 0:
            oxi_state = None
        elif oxi_state % 1 == 0:
            oxi_state = '{:+d}'.format(oxi_state)
        else:
            oxi_state = '{:+.2f}'.format(oxi_state)
    else:
        oxi_state = None

    formatted_element = specie.name

    if use_sym_label:
        formatted_element += sym_label

    if use_oxi_state and oxi_state:
        if latexify:
            oxi_state = "^{{{}}}".format(oxi_state)

        formatted_element += oxi_state

    return formatted_element


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
