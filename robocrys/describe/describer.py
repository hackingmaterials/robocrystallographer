"""This module provides a class for generating descriptions of condensed structure
data.

Todo:
    * Indicate when molecules have been simplified in the mineral description.
    * Handle distortion in connected polyhedra description.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import inflect
from pymatgen.util.string import htmlify, latexify, latexify_spacegroup, unicodeify

from robocrys.describe.adapter import DescriptionAdapter
from robocrys.util import (
    dimensionality_to_shape,
    geometry_to_polyhedra,
    get_el,
    get_formatted_el,
    htmlify_spacegroup,
    polyhedra_plurals,
    unicodeify_spacegroup,
)

if TYPE_CHECKING:
    from typing import Any


en = inflect.engine()


class StructureDescriber:

    def __init__(
        self,
        describe_mineral: bool = True,
        describe_component_makeup: bool = True,
        describe_components: bool = True,
        describe_symmetry_labels: bool = True,
        describe_oxidation_states: bool = True,
        describe_bond_lengths: bool = True,
        bond_length_decimal_places: int = 2,
        distorted_tol: float = 0.6,
        only_describe_cation_polyhedra_connectivity: bool = True,
        only_describe_bonds_once: bool = True,
        fmt: str = "raw",
        return_parts: bool = False,
    ):
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
            only_describe_bonds_once: Whether only describe bond lengths once.
                For example, don't describe the bond lengths from Pb to I and
                also from I to Pb.
            fmt: How to format element strings, formulas and spacegroups.
                Options are:

                - "raw" (default): Don't apply special formatting (e.g. "SnO2").
                - "unicode": Format super/subscripts using unicode characters
                  (e.g. SnO₂).
                - "latex": Use LaTeX markup for formatting (e.g. "SnO$_2$").
                - "html": Use html markup for formatting
                  (e.g. "SnO<sub>2</sub>").
            return_parts: Whether to return the individual parts of the
                description as a :obj:`dict`, or the whole description as a
                :obj:`str`.
        """
        self.distorted_tol = distorted_tol
        self.describe_mineral = describe_mineral
        self.describe_component_makeup = describe_component_makeup
        self.describe_components = describe_components
        self.describe_symmetry_labels = describe_symmetry_labels
        self.describe_oxidation_state = describe_oxidation_states
        self.describe_bond_lengths = describe_bond_lengths
        self.bond_length_decimal_places = bond_length_decimal_places
        self.cation_polyhedra_only = only_describe_cation_polyhedra_connectivity
        self.only_describe_bonds_once = only_describe_bonds_once
        self.fmt = fmt
        self.return_parts = return_parts
        self.angle_decimal_places = 0

        if fmt == "latex":
            self.angstrom = r"$\AA$"
            self.degree = r"$^{\circ}$"
        else:
            self.angstrom = "Å"
            self.degree = "°"

        self._da: DescriptionAdapter
        self._seen_bonds: set = set()

    def describe(self, condensed_structure: dict[str, Any]) -> str | dict[str, str]:
        """Convert a condensed structure into a text description.

        Args:
            condensed_structure: The condensed structure data, formatted as
                produced by :meth:`StructureCondenser.condense_structure`.

        Returns:
            A description of the structure. If
            :attr:`StructureDescriber.return_parts` is ``False``, the
            description will be returned as a :obj:`str`. If it is equal to
            ``True``, the description will be returned as a :obj:`dict` with the
            keys 'mineral', 'component_makeup' and 'components', each containing
            the relevant part of the description.
        """
        self._da = DescriptionAdapter(condensed_structure)
        self._seen_bonds = set()

        description = {}

        if self.describe_mineral:
            description["mineral"] = self.get_mineral_description()

        if self.describe_component_makeup:
            description["component_makeup"] = self.get_component_makeup_summary()

        if self.describe_components:
            description["components"] = self.get_all_component_descriptions()

        if not self.return_parts:
            return " ".join(
                description[part]
                for part in ["mineral", "component_makeup", "components"]
                if part in description and description[part] != ""
            )
        return description

    def get_mineral_description(self) -> str:
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
        spg_symbol = self._da.spg_symbol
        formula = self._da.formula
        if self.fmt == "latex":
            spg_symbol = latexify_spacegroup(self._da.spg_symbol)
            formula = latexify(formula)

        elif self.fmt == "unicode":
            spg_symbol = unicodeify_spacegroup(self._da.spg_symbol)
            formula = unicodeify(formula)

        elif self.fmt == "html":
            spg_symbol = htmlify_spacegroup(self._da.spg_symbol)
            formula = htmlify(formula)

        if mineral_name := get_mineral_name(self._da.mineral):
            # replace latex-like characters with unicode
            latex_reps = {
                "$\\mu$": "\u03BC",
                "\\'{c}": "\u0107",
            }
            if self.fmt in {"html", "unicode"} and any(
                k in mineral_name for k in latex_reps
            ):
                for k, v in latex_reps.items():
                    mineral_name = mineral_name.replace(k, v)

            desc = f"{formula} is {mineral_name} structured and"
        else:
            desc = f"{formula}"

        desc += " crystallizes in the {} {} space group.".format(
            self._da.crystal_system, spg_symbol
        )
        return desc

    def get_component_makeup_summary(self) -> str:
        """Gets a summary of the makeup of components in a structure.

        Returns:
            A description of the number of components and their dimensionalities
            and orientations.
        """
        component_groups = self._da.get_component_groups()

        if (
            len(component_groups) == 1
            and component_groups[0].count == 1
            and component_groups[0].dimensionality == 3
        ):
            desc = ""

        else:
            if self._da.dimensionality == 3:
                desc = "The structure consists of "
            else:
                desc = (
                    f"The structure is {en.number_to_words(self._da.dimensionality)}"  # type: ignore[arg-type]
                    "-dimensional and consists of "
                )

            component_makeup_summaries = []
            nframeworks = len(
                [
                    c
                    for g in component_groups
                    for c in g.components
                    if c.dimensionality == 3
                ]
            )
            for component_group in component_groups:
                if nframeworks == 1 and component_group.dimensionality == 3:
                    s_count: str = "a"
                else:
                    s_count = en.number_to_words(component_group.count)  # type: ignore[assignment]

                dimensionality = component_group.dimensionality

                if component_group.molecule_name:
                    shape = "atom" if component_group.nsites == 1 else "molecule"
                    shape = en.plural(shape, s_count)
                    formula = component_group.molecule_name
                else:
                    shape = en.plural(dimensionality_to_shape[dimensionality], s_count)
                    formula = component_group.formula

                if self.fmt == "latex":
                    formula = latexify(formula)
                elif self.fmt == "unicode":
                    formula = unicodeify(formula)
                elif self.fmt == "html":
                    formula = htmlify(formula)

                comp_desc = f"{s_count} {formula} {shape}"

                if component_group.dimensionality in [1, 2]:
                    orientations = list(
                        {str(c.orientation) for c in component_group.components}
                    )
                    s_direction = en.plural("direction", len(orientations))
                    comp_desc += " oriented in the {} {}".format(
                        en.join(orientations), s_direction
                    )

                component_makeup_summaries.append(comp_desc)

            if nframeworks == 1 and len(component_makeup_summaries) > 1:
                # when there is a single framework, make the description read
                # "... and 8 Sn atoms inside a SnO2 framework" instead of
                # "..., 8 Sn atoms and one SnO2 framework"
                # This works because the component summaries are sorted by
                # dimensionality
                desc += en.join(component_makeup_summaries[:-1])
                desc += f" inside {component_makeup_summaries[-1]}."
            else:
                desc += en.join(component_makeup_summaries) + "."
        return desc

    def get_all_component_descriptions(self) -> str:
        """Gets the descriptions of all components in the structure.

        Returns:
            A description of all components in the structure.
        """
        if len(self._da.components) == 1:
            return self.get_component_description(
                self._da.get_component_groups()[0].components[0].index,
                single_component=True,
            )

        component_groups = self._da.get_component_groups()

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

                if self.fmt == "latex":
                    formula = latexify(formula)
                elif self.fmt == "unicode":
                    formula = unicodeify(formula)
                elif self.fmt == "html":
                    formula = htmlify(formula)

                if group_count == component_count:
                    s_filler = "the" if group_count == 1 else "each"
                else:
                    s_filler = f"{en.number_to_words(component_count)} of the"
                    shape = en.plural(shape)

                desc = f"In {s_filler} {formula} {shape}, "
                desc += self.get_component_description(component.index)

                component_descriptions.append(desc)

        return " ".join(component_descriptions)

    def get_component_description(
        self, component_index: int, single_component: bool = False
    ) -> str:
        """Gets the descriptions of all sites in a component.

        Args:
            component_index: The index of the component
            single_component: Whether the structure contains only a single
                component.

        Returns:
            The description for all sites in the components.
        """
        desc = []
        first_group = True
        for site_group in self._da.get_component_site_groups(component_index):
            if len(site_group.sites) == 1:
                desc.append(self.get_site_description(site_group.sites[0]))

            else:
                element = get_formatted_el(
                    site_group.element,
                    "",
                    use_oxi_state=self.describe_oxidation_state,
                    use_sym_label=False,
                    fmt=self.fmt,
                )

                s_there = "there" if first_group and not single_component else "There"

                s_count = en.number_to_words(len(site_group.sites))  # type: ignore[arg-type]

                desc.append(f"{s_there} are {s_count} inequivalent {element} sites.")

                for i, site in enumerate(site_group.sites):
                    s_ordinal = en.number_to_words(en.ordinal(i + 1))  # type: ignore[arg-type]
                    desc.append(f"In the {s_ordinal} {element} site,")
                    desc.append(self.get_site_description(site))

            first_group = False
        return " ".join(desc)

    def get_site_description(self, site_index: int) -> str:
        """Gets a description of the geometry and bonding of a site.

        If the site likeness (order parameter) is less than ``distorted_tol``,
        "distorted" will be added to the geometry description.

        Args:
            site_index: An inequivalent site index.

        Returns:
            A description of the geometry and bonding of a site.
        """
        site = self._da.sites[site_index]

        if site["poly_formula"] and (
            self.cation_polyhedra_only or "+" in site["element"]
        ):
            desc = self._get_poly_site_description(site_index)
            tilt_desc = self.get_octahedral_tilt_description(site_index)
            if tilt_desc:
                desc += ". " + tilt_desc
        else:
            element = get_formatted_el(
                site["element"],
                self._da.sym_labels[site_index],
                use_oxi_state=self.describe_oxidation_state,
                use_sym_label=self.describe_symmetry_labels,
                fmt=self.fmt,
            )

            if site["geometry"]["likeness"] < self.distorted_tol:
                s_geometry = "distorted "
            else:
                s_geometry = ""
            s_geometry += site["geometry"]["type"]

            desc = f"{element} is bonded in {en.a(s_geometry)} geometry to "
            desc += self._get_nearest_neighbor_description(site_index)

        bond_length_desc = self._get_nearest_neighbor_bond_length_descriptions(
            site_index
        )
        if bond_length_desc:
            desc += ". " + bond_length_desc
        else:
            desc += "."

        return desc

    def _get_poly_site_description(self, site_index: int):
        """Gets a description of a connected polyhedral site.

        If the site likeness (order parameter) is less than ``distorted_tol``,
        "distorted" will be added to the geometry description.

        Args:
            site_index: An inequivalent site index.

        Returns:
            A description the a polyhedral site, including connectivity.
        """
        site = self._da.sites[site_index]
        nnn_details = self._da.get_next_nearest_neighbor_details(
            site_index, group=not self.describe_symmetry_labels
        )

        from_element = get_formatted_el(
            site["element"],
            self._da.sym_labels[site_index],
            use_oxi_state=self.describe_oxidation_state,
            use_sym_label=self.describe_symmetry_labels,
            fmt=self.fmt,
        )

        from_poly_formula = site["poly_formula"]
        if self.fmt == "latex":
            from_poly_formula = latexify(from_poly_formula)
        elif self.fmt == "unicode":
            from_poly_formula = unicodeify(from_poly_formula)
        elif self.fmt == "html":
            from_poly_formula = htmlify(from_poly_formula)

        s_from_poly_formula = get_el(site["element"]) + from_poly_formula

        if site["geometry"]["likeness"] < self.distorted_tol:
            s_distorted = "distorted "
        else:
            s_distorted = ""
        s_polyhedra = geometry_to_polyhedra[site["geometry"]["type"]]
        s_polyhedra = polyhedra_plurals[s_polyhedra]

        nn_desc = self._get_nearest_neighbor_description(site_index)
        desc = f"{from_element} is bonded to {nn_desc} to form "

        # handle the case we were are connected to the same type of polyhedra
        if (
            nnn_details[0].element == site["element"]
            and len(
                {(nnn_site.element, nnn_site.poly_formula) for nnn_site in nnn_details}
            )
        ) == 1:
            connectivities = list({nnn_site.connectivity for nnn_site in nnn_details})
            s_mixture = "a mixture of " if len(connectivities) != 1 else ""
            s_connectivities = en.join(connectivities)

            desc += "{}{}{}-sharing {} {}".format(
                s_mixture,
                s_distorted,
                s_connectivities,
                s_from_poly_formula,
                s_polyhedra,
            )
            return desc

        # otherwise loop through nnn connectivities and describe individually
        desc += "{}{} {} that share ".format(
            s_distorted, s_from_poly_formula, s_polyhedra
        )
        nnn_descriptions = []
        for nnn_site in nnn_details:
            to_element = get_formatted_el(
                nnn_site.element,
                nnn_site.sym_label,
                use_oxi_state=False,
                use_sym_label=self.describe_symmetry_labels,
            )

            to_poly_formula = nnn_site.poly_formula
            if self.fmt == "latex":
                to_poly_formula = latexify(to_poly_formula)
            elif self.fmt == "unicode":
                to_poly_formula = unicodeify(to_poly_formula)
            elif self.fmt == "html":
                to_poly_formula = htmlify(to_poly_formula)

            to_poly_formula = to_element + to_poly_formula
            to_shape = geometry_to_polyhedra[nnn_site.geometry]

            if len(nnn_site.sites) == 1 and nnn_site.count != 1:
                s_equivalent = " equivalent "
            else:
                s_equivalent = " "

            if nnn_site.count == 1:
                s_an = f" {en.an(nnn_site.connectivity)}"
            else:
                s_an = ""
                to_shape = polyhedra_plurals[to_shape]

            nnn_descriptions.append(
                "{}{} with {}{}{} {}".format(
                    s_an,
                    en.plural(nnn_site.connectivity, nnn_site.count),
                    en.number_to_words(nnn_site.count),
                    s_equivalent,
                    to_poly_formula,
                    to_shape,
                )
            )

        return desc + en.join(nnn_descriptions)

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
            site_index, group=not self.describe_symmetry_labels
        )

        last_count = 0
        nn_descriptions = []
        for nn_site in nn_details:
            element = get_formatted_el(
                nn_site.element,
                nn_site.sym_label,
                use_oxi_state=self.describe_oxidation_state,
                use_sym_label=self.describe_symmetry_labels,
                fmt=self.fmt,
            )

            if len(nn_site.sites) == 1 and nn_site.count != 1:
                s_equivalent = " equivalent "
            else:
                s_equivalent = " "

            nn_descriptions.append(
                "{}{}{}".format(
                    en.number_to_words(nn_site.count), s_equivalent, element
                )
            )
            last_count = nn_site.count

        s_atoms = "atom" if last_count == 1 else "atoms"
        return f"{en.join(nn_descriptions)} {s_atoms}"

    def _get_nearest_neighbor_bond_length_descriptions(self, site_index: int) -> str:
        """Gets the descriptions of the bond lengths for nearest neighbor sites.

        Args:
            site_index: An inequivalent site index.

        Returns:
            A description of the nearest neighbor bond lengths.
        """
        if not self.describe_bond_lengths:
            return ""

        nn_details = self._da.get_nearest_neighbor_details(
            site_index, group=not self.describe_symmetry_labels
        )

        bond_descriptions = []
        for nn_site in nn_details:
            bond_descriptions.append(
                self.get_bond_length_description(site_index, nn_site.sites)
            )

        # filter empty bond length description strings
        return " ".join(filter(lambda x: x, bond_descriptions))

    def get_bond_length_description(self, from_site: int, to_sites: list[int]) -> str:
        """Gets a description of the bond lengths between two sets of sites.

        Args:
            from_site: An inequivalent site index.
            to_sites: A :obj:`list` of site indices. The site indices should
                all be for the same element.

        Returns:
            A description of the bond lengths or an empty string if
            :attr:`StructureDescriber.only_describe_bonds_once` is ``True`` and
            all all bond lengths have already been described.
        """
        if self.only_describe_bonds_once:
            to_sites = self._filter_seen_bonds(from_site, to_sites)

            if not to_sites:
                return ""

        from_element = get_formatted_el(
            self._da.elements[from_site],
            self._da.sym_labels[from_site],
            use_oxi_state=False,
            use_sym_label=self.describe_symmetry_labels,
        )
        to_element = get_formatted_el(
            self._da.elements[to_sites[0]],
            self._da.get_sym_label(to_sites),
            use_oxi_state=False,
            use_sym_label=self.describe_symmetry_labels,
        )

        dists = self._da.get_distance_details(from_site, to_sites)

        # if only one bond length
        if len(dists) == 1:
            return "The {}-{} bond length is {}.".format(
                from_element, to_element, self._distance_to_string(dists[0])
            )

        discrete_bond_lengths = self._rounded_bond_lengths(dists)

        # if multiple bond lengths but they are all the same
        if len(set(discrete_bond_lengths)) == 1:
            s_intro = "Both" if len(discrete_bond_lengths) == 2 else "All"
            return "{} {}-{} bond lengths are {}.".format(
                s_intro, from_element, to_element, self._distance_to_string(dists[0])
            )

        # if two sets of bond lengths
        if len(set(discrete_bond_lengths)) == 2:
            small = min(discrete_bond_lengths)
            s_small_count = en.number_to_words(discrete_bond_lengths.count(small))  # type: ignore[arg-type]
            big = max(discrete_bond_lengths)
            s_big_count = en.number_to_words(discrete_bond_lengths.count(big))  # type: ignore[arg-type]

            s_length = en.plural("length", s_big_count)

            return (
                "There {} {} shorter ({}) and {} longer ({}) {}-{} bond {}."
            ).format(
                en.plural_verb("is", s_small_count),
                s_small_count,
                self._distance_to_string(small),
                s_big_count,
                self._distance_to_string(big),
                from_element,
                to_element,
                s_length,
            )

        # otherwise just detail the spread of bond lengths
        return ("There are a spread of {}-{} bond distances ranging from {}.").format(
            from_element,
            to_element,
            self._distance_range_to_string(
                min(discrete_bond_lengths), max(discrete_bond_lengths)
            ),
        )

    def get_octahedral_tilt_description(
        self,
        site_index: int,
    ) -> str:
        """Gets a description of octahedral tilting angles between two sites.

        Currently only implemented for corner-sharing octahedra.
        Will throw an error if the two sites are not next nearest neighbors.

        Args:
            site_index: An inequivalent site index.

        Returns:
            A description of the octahedral tilting angles.
        """
        nnn_details = self._da.get_next_nearest_neighbor_details(
            site_index, group=not self.describe_symmetry_labels
        )
        to_sites = [
            site
            for nnn_site in nnn_details
            for site in nnn_site.sites
            if nnn_site.geometry == "octahedral" and nnn_site.connectivity == "corner"
        ]

        angles = self._da.get_angle_details(site_index, to_sites, "corner")
        discrete_angles = list(set(self._rounded_angles(angles)))
        tilts = [abs(180 - angle) for angle in discrete_angles]

        if not tilts:
            return ""

        # if only one bond length
        if len(tilts) == 1:
            if tilts[0] == 0:
                return "The corner-sharing octahedra are not tilted"

            return "The corner-sharing octahedral tilt angles are {}".format(
                self._angle_to_string(tilts[0])
            )

        # otherwise just detail the spread of bond lengths
        return "The corner-sharing octahedral tilt angles range from {}".format(
            self._angle_range_to_string(min(tilts), max(tilts))
        )

    def _filter_seen_bonds(self, from_site: int, to_sites: list[int]) -> list[int]:
        """Filter the list of to_sites to only include unseen bonds.

        Args:
            from_site: An inequivalent site index.
            to_sites: A :obj:`list` of site indices. The site indices should
                all be for the same element.

        Returns:
            The list of unseen bonds.
        """
        # get a list of tuples of (from_site, to_site) for all to_sites
        bonds = [(f, t) for f, t in zip([from_site] * len(to_sites), to_sites)]

        # use frozen set as it is invariant to the order of the sites
        not_seen = [x for x in bonds if frozenset(x) not in self._seen_bonds]

        # only describe the bonds between unseen site pairs
        not_seen_to_sites = []
        for from_site, to_site in not_seen:
            not_seen_to_sites.append(to_site)
            self._seen_bonds.add(frozenset((from_site, to_site)))

        return not_seen_to_sites

    def _rounded_bond_lengths(self, data: list[float]) -> tuple[float, ...]:
        """Function to round bond lengths to a number of decimal places."""
        return tuple(
            float("{:.{}f}".format(x, self.bond_length_decimal_places)) for x in data
        )

    def _distance_to_string(self, distance: float) -> str:
        """Utility function to round a distance and add an Angstrom symbol."""
        return "{:.{}f} {}".format(
            distance, self.bond_length_decimal_places, self.angstrom
        )

    def _distance_range_to_string(self, dist_a: float, dist_b: float) -> str:
        """Utility function to format a range of distances."""
        return "{:.{}f}-{:.{}f} {}".format(
            dist_a,
            self.bond_length_decimal_places,
            dist_b,
            self.bond_length_decimal_places,
            self.angstrom,
        )

    def _rounded_angles(self, data: list[float]) -> tuple[float, ...]:
        """Function to round angles to a number of decimal places."""
        return tuple(
            float("{:.{}f}".format(x, self.angle_decimal_places)) for x in data
        )

    def _angle_to_string(self, angle: float) -> str:
        """Utility function to round a distance and add an Angstrom symbol."""
        return "{:.{}f}{}".format(angle, self.angle_decimal_places, self.degree)

    def _angle_range_to_string(self, angle_a: float, angle_b: float) -> str:
        """Utility function to format a range of distances."""
        return "{:.{}f}-{:.{}f}{}".format(
            angle_a,
            self.angle_decimal_places,
            angle_b,
            self.angle_decimal_places,
            self.degree,
        )


def get_mineral_name(mineral_dict: dict[str, Any]) -> str | None:
    """Get the mineral name from a mineral dictionary.

    Args:
        mineral_dict: The mineral dictionary from the condensed description.

    Returns:
        If ``mineral_dict["type"]`` is set, the mineral name will be returned as
        a string, else ``None`` will be returned.
    """
    if mineral_dict["type"]:
        if not mineral_dict["n_species_type_match"]:
            suffix = "-derived"
        elif mineral_dict["distance"] >= 0:
            suffix = "-like"
        else:
            suffix = ""

        return "{}{}".format(mineral_dict["type"], suffix)

    return None
