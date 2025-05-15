"""This module implements functions for handling structure components."""

from __future__ import annotations

from copy import deepcopy
from typing import TYPE_CHECKING, Any

import networkx as nx
import numpy as np
from monty.fractions import gcd
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import get_el_sp
from pymatgen.core.structure import PeriodicSite, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.string import formula_double_format

from robocrys.condense.fingerprint import get_structure_fingerprint
from robocrys.util import common_formulas

if TYPE_CHECKING:
    Component = dict[str, Any]


def get_structure_inequiv_components(
    components: list[Component],
    use_structure_graph: bool = False,
    fingerprint_tol: float = 0.01,
) -> list[Component]:
    """Gets and counts the structurally inequivalent components.

    Supports matching through StructureMatcher or by a combined structure graph/
    site fingerprint approach. For the latter method, the component data has to
    have been generated with ``inc_graph=True``.

    Args:
        components: A list of structure components, generated using
            :obj:`pymatgen.analysis.dimensionality.get_structure_components`.
            If ``use_structure_graph=True``, the components should be generated
            with ``inc_graph=True``.
        use_structure_graph: Whether to use the bonding graph and site
            fingerprints to determine if components are structurally equivalent.
            If ``False``,
            :obj:`pymatgen.analysis.structure_matcher.StructureMatcher` will be
            used to compare the components.
        fingerprint_tol: The fingerprint tolerance to determine whether two
            components have a matching fingerprint. This option is ignored if
            ``use_structure_graph=False``.

    Returns:
        A list of the structurally inequivalent components. Any duplicate
        components will only be returned once. The component objects are in the
        same format is given by
        :obj:`pymatgen.analysis.dimensionality.get_structure_components` but
        have an additional field:

        - ``"count"`` (:obj:`int`): The number of times this component appears
          in the structure.
    """
    components = deepcopy(components)

    for component in components:
        component["count"] = 1

    if use_structure_graph:
        # check fingerprints match and components are isomorphic.
        fingerprints = [
            get_structure_fingerprint(c["structure_graph"].structure)
            for c in components
        ]

        seen_components = [components[0]]
        seen_fingers = [fingerprints[0]]
        for component, fingerprint in zip(components[1:], fingerprints[1:]):
            graph_match = [
                components_are_isomorphic(component, c) for c in seen_components
            ]
            finger_match = [
                np.linalg.norm(fingerprint - c) < fingerprint_tol for c in seen_fingers
            ]
            structure_match = [i and f for i, f in zip(graph_match, finger_match)]

            if any(structure_match):
                # there should only ever be a single match so we take index of
                # the first match and increment the component count
                loc = np.where(structure_match)[0][0]
                seen_components[loc]["count"] += 1
            else:
                seen_components.append(component)
                seen_fingers.append(fingerprint)

    else:
        sm = StructureMatcher()
        seen_components = [components[0]]

        for component in components[1:]:
            structure_match = [
                sm.fit(
                    component["structure_graph"].structure,
                    c["structure_graph"].structure,
                )
                for c in seen_components
            ]
            if any(structure_match):
                # there should only ever be a single match so we take index of
                # the first match and increment the component count
                loc = np.where(structure_match)[0][0]
                seen_components[loc]["count"] += 1
            else:
                seen_components.append(component)

    return seen_components


def components_are_isomorphic(
    component_a: Component, component_b: Component, use_weights: bool = False
):
    """Determines whether the graphs of two components are isomorphic.

    Only takes into account graph connectivity and not local geometry (e.g. bond
    angles and distances).

    Args:
        component_a: The first component.
        component_b: The second component.
        use_weights: Whether to use the graph edge weights in comparing graphs.

    Returns:
        Whether the components are isomorphic.
    """

    def node_match(n1, n2):
        return n1["specie"] == n2["specie"]

    def edge_match(e1, e2):
        if use_weights:
            return e1["weight"] == e2["weight"]
        return True

    graph_a = component_a["structure_graph"].graph
    graph_b = component_b["structure_graph"].graph

    species_a = {
        n: {"specie": str(component_a["structure_graph"].structure[n].specie)}
        for n in graph_a
    }
    species_b = {
        n: {"specie": str(component_b["structure_graph"].structure[n].specie)}
        for n in graph_b
    }

    nx.set_node_attributes(graph_a, species_a)
    nx.set_node_attributes(graph_b, species_b)

    return nx.is_isomorphic(
        graph_a, graph_b, node_match=node_match, edge_match=edge_match
    )


def get_sym_inequiv_components(
    components: list[Component], spg_analyzer: SpacegroupAnalyzer
) -> list[Component]:
    """Gets and counts the symmetrically inequivalent components.

    Component data has to have been generated with ``inc_site_ids=True``.

    Args:
        components: A list of structure components, generated using
            :obj:`pymatgen.analysis.dimensionality.get_structure_components`,
            with ``inc_site_ids=True``.
        spg_analyzer: A `pymatgen.symmetry.analyzer.SpacegroupAnalyzer` analyzer
            object for the structure containing the components.

    Returns:
        A list of the symmetrically inequivalent components. Any duplicate
        components will only be returned once. The component objects are in the
        same format is given by
        :obj:`pymatgen.analysis.dimensionality.get_structure_components` but
        the additional property:

        - ``"count"`` (:obj:`int`): The number of times this component appears
          in the structure.
    """
    components = deepcopy(components)

    sym_inequiv_components: dict[frozenset, Component] = {}
    equivalent_atoms = spg_analyzer.get_symmetry_dataset().equivalent_atoms

    for component in components:
        sym_indices = frozenset(equivalent_atoms[x] for x in component["site_ids"])

        # if two components are composed of atoms that are symmetrically
        # equivalent they are the same.
        if sym_indices in sym_inequiv_components:
            sym_inequiv_components[sym_indices]["count"] += 1
            continue

        component["count"] = 1
        sym_inequiv_components[sym_indices] = component

    return list(sym_inequiv_components.values())


def get_formula_inequiv_components(
    components: list[Component],
    use_iupac_formula: bool = True,
    use_common_formulas: bool = True,
) -> list[Component]:
    """Gets and counts the inequivalent components based on their formuula.

    Note that the counting of compounds is different to in
    ``get_sym_inequiv_equivalent``. I.e. the count is not the number of
    components with the same formula. For example, the count of the formula
    "GaAs" in a system with two Ga2As2 components would be 4.

    Args:
        components: A list of structure components, generated using
            :obj:`pymatgen.analysis.dimensionality.get_structure_components`,
            with ``inc_site_ids=True``.
        use_iupac_formula (bool, optional): Whether to order formulas by the
            iupac "electronegativity" series, defined in Table VI of
            "Nomenclature of Inorganic Chemistry (IUPAC Recommendations 2005)".
            This ordering effectively follows the groups and rows of the
            periodic table, except the Lanthanides, Actanides and hydrogen. If
            set to ``False``, the elements will be ordered according to the
            electronegativity values.
        use_common_formulas: Whether to use the database of common formulas.
            The common formula will be used preferentially to the iupac or
            reduced formula.

    Returns:
        A list of the compositionally inequivalent components. Any duplicate
        components will only be returned once. The component objects are in the
        same format is given by
        :obj:`pymatgen.analysis.dimensionality.get_structure_components` but
        have two additional fields:

        - ``"count"`` (:obj:`int`): The number of formula units of this
          component. Note, this is not the number of components with the same
          formula. For example, the count of the formula "GaAs" in a system
          with two Ga2As2 components would be 4.
        - ``"formula"`` (``list[int]``): The reduced formula of the component.
    """
    components = deepcopy(components)
    inequiv_components: dict[str, Component] = {}

    for component in components:
        formula, factor = get_component_formula_and_factor(
            component,
            use_iupac_formula=use_iupac_formula,
            use_common_formulas=use_common_formulas,
        )

        # if two components have the same composition we treat them as the same
        if formula in inequiv_components:
            inequiv_components[formula]["count"] += factor
            continue

        component["count"] = factor
        component["formula"] = formula

        inequiv_components[formula] = component

    return list(inequiv_components.values())


def filter_molecular_components(
    components: list[Component],
) -> tuple[list[Component], list[Component]]:
    """Separate list of components into molecular and non-molecular components.

    Args:
        components: A list of structure components, generated using
            :obj:`pymatgen.analysis.dimensionality.get_structure_components`.

    Returns:
        The filtered components as a tuple of ``(molecular_components,
        other_components)``.
    """
    molecular_components = [c for c in components if c["dimensionality"] == 0]
    other_components = [c for c in components if c["dimensionality"] != 0]

    return molecular_components, other_components


def get_reconstructed_structure(
    components: list[Component], simplify_molecules: bool = True
) -> Structure:
    """Reconstructs a structure from a list of components.

    Has the option to simplify molecular components into a single site
    positioned at the centre of mass of the molecular. If using this option,
    the components must have been generated with ``inc_molecule_graph=True``.

    Args:
        components: A list of structure components, generated using
            :obj:`pymatgen.analysis.dimensionality.get_structure_components`,
            with ``inc_molecule_graph=True``.
        simplify_molecules: Whether to simplify the molecular components into
            a single site positioned at the centre of mass of the molecule.

    Returns:
        The reconstructed structure.
    """
    mol_sites = []
    other_sites = []

    if simplify_molecules:
        mol_components, components = filter_molecular_components(components)

        if mol_components:
            lattice = mol_components[0]["structure_graph"].structure.lattice

            mol_sites = [
                PeriodicSite(
                    c["structure_graph"].structure[0].specie,
                    c["molecule_graph"].molecule.center_of_mass,
                    lattice,
                    coords_are_cartesian=True,
                )
                for c in mol_components
            ]

    if components:
        other_sites = [
            site for c in components for site in c["structure_graph"].structure
        ]

    return Structure.from_sites(other_sites + mol_sites)


def get_component_formula_and_factor(
    component: Component,
    use_iupac_formula: bool = True,
    use_common_formulas: bool = True,
) -> tuple[str, int]:
    """Gets the reduced formula and factor of a single component.

    Args:
        component: A structure component, generated using
            :obj:`pymatgen.analysis.dimensionality.get_structure_components`.
        use_iupac_formula (bool, optional): Whether to order formulas by the
            iupac "electronegativity" series, defined in Table VI of
            "Nomenclature of Inorganic Chemistry (IUPAC Recommendations 2005)".
            This ordering effectively follows the groups and rows of the
            periodic table, except the Lanthanides, Actanides and hydrogen. If
            set to ``False``, the elements will be ordered according to the
            electronegativity values.
        use_common_formulas: Whether to use the database of common formulas.
            The common formula will be used preferentially to the iupac or
            reduced formula.

    Returns:
        The formula and factor of the component.
    """
    formula, factor = component[
        "structure_graph"
    ].structure.composition.get_reduced_formula_and_factor(
        iupac_ordering=use_iupac_formula
    )

    reduced_formula = component["structure_graph"].structure.composition.reduced_formula
    if use_common_formulas and reduced_formula in common_formulas:
        formula = common_formulas[reduced_formula]
    return formula, factor


def get_component_formula(
    component: Component,
    use_iupac_formula: bool = True,
    use_common_formulas: bool = True,
) -> str:
    """Gets the reduced formula of a single component.

    Args:
        component: A structure component, generated using
            :obj:`pymatgen.analysis.dimensionality.get_structure_components`.
        use_iupac_formula (bool, optional): Whether to order formulas by the
            iupac "electronegativity" series, defined in Table VI of
            "Nomenclature of Inorganic Chemistry (IUPAC Recommendations 2005)".
            This ordering effectively follows the groups and rows of the
            periodic table, except the Lanthanides, Actanides and hydrogen. If
            set to ``False``, the elements will be ordered according to the
            electronegativity values.
        use_common_formulas: Whether to use the database of common formulas.
            The common formula will be used preferentially to the iupac or
            reduced formula.

    Returns:
        The formula and factor of the component.
    """
    return get_component_formula_and_factor(
        component,
        use_iupac_formula=use_iupac_formula,
        use_common_formulas=use_common_formulas,
    )[0]


def get_formula_from_components(
    components: list[Component],
    molecules_first: bool = False,
    use_iupac_formula: bool = True,
    use_common_formulas: bool = True,
) -> str:
    """Reconstructs a chemical formula from structure components.

    The chemical formulas for the individual components will be grouped
    together. If two components share the same composition, they will be
    treated as equivalent.

    Args:
        components: A list of structure components, generated using
            :obj:`pymatgen.analysis.dimensionality.get_structure_components`.
        molecules_first: Whether to put any molecules (zero-dimensional
            components) at the beginning of the formula.
        use_iupac_formula (bool, optional): Whether to order formulas by the
            iupac "electronegativity" series, defined in Table VI of
            "Nomenclature of Inorganic Chemistry (IUPAC Recommendations 2005)".
            This ordering effectively follows the groups and rows of the
            periodic table, except the Lanthanides, Actanides and hydrogen. If
            set to ``False``, the elements will be ordered according to the
            electronegativity values.
        use_common_formulas: Whether to use the database of common formulas.
            The common formula will be used preferentially to the iupac or
            reduced formula.

    Returns:
        The chemical formula.
    """

    def order(comp_formula):
        composition = Composition(comp_formula)
        if use_iupac_formula:
            return sum(
                [get_el_sp(s).iupac_ordering for s in composition.elements]
            ) / len(composition.elements)
        return composition.average_electroneg

    components = get_formula_inequiv_components(
        components,
        use_iupac_formula=use_iupac_formula,
        use_common_formulas=use_common_formulas,
    )

    if molecules_first:
        mol_comps, other_comps = filter_molecular_components(components)
    else:
        mol_comps = []
        other_comps = components

    formulas = sorted([c["formula"] for c in mol_comps], key=order) + sorted(
        [c["formula"] for c in other_comps], key=order
    )

    # if components include special formulas, then the count can be 0.5
    # therefore if any non integer amounts we can just use a factor of 2
    all_int = all(v["count"] % 1 == 0 for v in components)
    prefactor = 1 if all_int else 2
    form_count_dict = {c["formula"]: int(c["count"] * prefactor) for c in components}

    # the following is based on ``pymatgen.core.composition.reduce_formula``
    num_comps = len(formulas)
    factor = abs(gcd(*(form_count_dict.values())))

    reduced_form = []
    for i in range(0, num_comps):
        formula = formulas[i]
        norm_amt = form_count_dict[formula] * 1.0 / factor
        formatted_formula = formula if norm_amt == 1 else f"({formula})"
        reduced_form.append(formatted_formula)
        reduced_form.append(str(formula_double_format(norm_amt)))

    return "".join(reduced_form)


def components_are_vdw_heterostructure(components: list[Component]) -> bool:
    """Whether a list of components form a van der Waals heterostructure.

    A heterostructure is defined here as a structure with more than one
    formula inequivalent 2D component.

    Args:
        components: A list of structure components, generated using
            :obj:`pymatgen.analysis.dimensionality.get_structure_components`.

    Returns:
        Whether the list of components from a heterostructure.
    """
    components = get_formula_inequiv_components(components)

    return bool(len([c for c in components if c["dimensionality"] == 2]))


def get_vdw_heterostructure_information(
    components: list[Component],
    use_iupac_formula: bool = True,
    use_common_formulas: bool = True,
    inc_ordered_components: bool = False,
    inc_intercalants: bool = False,
) -> dict[str, Any]:
    """Gets information about ordering of components in a vdw heterostructure.

    Args:
        components: A list of structure components, generated using
            :obj:`pymatgen.analysis.dimensionality.get_structure_components`
            with ``inc_orientation=True``.
        use_iupac_formula (bool, optional): Whether to order formulas by the
            iupac "electronegativity" series, defined in Table VI of
            "Nomenclature of Inorganic Chemistry (IUPAC Recommendations 2005)".
            This ordering effectively follows the groups and rows of the
            periodic table, except the Lanthanides, Actanides and hydrogen. If
            set to ``False``, the elements will be ordered according to the
            electronegativity values.
        use_common_formulas: Whether to use the database of common formulas.
            The common formula will be used preferentially to the iupac or
            reduced formula.
        inc_ordered_components: Whether to return a list of the ordered
            components. If False, just the component formulas will be returned.
        inc_intercalants: Whether to return a list of the intercalants. If
            False, just the intercalant formulas will be returned.

    Returns:
        Information on the heterostructure, as an :obj:`dict` with they keys:

        - ``"repeating_unit"`` (``list[str]``): A :obj:`List` of formulas of the
          smallest repeating series of components. For example. if the
          structure consists of A and B components ordered as "A B A B A B",
          the repeating unit is "A B".
        - ``"num_repetitions"`` (:obj:`int`): The number of repetitions of the
          repeating unit that forms the overall structure. For example. if
          the structure consists of A and B components ordered as
          "A B A B A B", the number of repetitions is 3.
        - ``"intercalant_formulas"`` (:obj:`list[str]`): The formulas of the
          intercalated compounds.
        - ``"ordered_components"`` (``list[component]``): If
          ``inc_ordered_components``, a :obj:`List` of components, ordered as
          they appear in the heteostructure stacking direction.
        - ``"intercalants"`` (``list[component]``: If ``inc_intercalants``, a
          :obj:`List` of intercalated components.
    """
    if not components_are_vdw_heterostructure(components):
        raise ValueError("Components do not form a heterostructure.")

    try:
        millers = {c["orientation"] for c in components if c["dimensionality"] == 2}
    except KeyError as exc:
        if "orientation" in str(exc):
            raise KeyError(
                "Components not generated with inc_orientation=True"
            ) from exc
        raise exc

    if len(millers) != 1:
        raise ValueError("2D components don't all have the same orientation.")

    cart_miller = (
        components[0]["structure_graph"]
        .structure.lattice.get_cartesian_coords(millers.pop())
        .tolist()
    )

    # plane is used to find the distances of all components along a certain axis
    # should use normal vector of plane to get exact distances but we just care
    # about relative distances.
    def distances_to_plane(points):
        return [np.dot(cart_miller, pp) for pp in points]

    min_distances = [
        min(distances_to_plane(c["structure_graph"].structure.cart_coords))
        for c in components
    ]

    # sort the components by distance to plane
    ordering = np.argsort(min_distances)
    ordered_components = [components[x] for x in ordering]

    # only consider the layered components formulae
    ordered_layers = [c for c in ordered_components if c["dimensionality"] == 2]
    ordered_layers_formula = [
        get_component_formula(
            c,
            use_iupac_formula=use_iupac_formula,
            use_common_formulas=use_common_formulas,
        )
        for c in ordered_layers
    ]
    num_layer_formulas = len(set(ordered_layers_formula))

    repeating_formula = ordered_layers_formula
    num_repetitions = 1

    # depending on the number of inequivalent formulae, there is a maximum
    # number of repetitions that can occur. To avoid unnecessary work we start
    # from this number of repetitions and move to 1 repetition (e.g. no
    # repetition)
    max_repetitions = int(np.floor(len(ordered_layers) / num_layer_formulas))
    for n in range(max_repetitions, 0, -1):
        if (
            all(
                len(set(ordered_layers_formula[i::num_layer_formulas])) == 1
                for i in range(n)
            )
            and len(ordered_layers) % n == 0
        ):
            repeating_formula = ordered_layers_formula[: int(len(ordered_layers) / n)]
            num_repetitions = n
            break

    intercalants = [c for c in components if c["dimensionality"] < 2]
    intercalant_formulas = [get_component_formula(c) for c in intercalants]

    data = {
        "repeating_unit": repeating_formula,
        "num_repetitions": num_repetitions,
        "intercalant_formulas": intercalant_formulas,
    }

    if inc_intercalants:
        data["intercalants"] = intercalants

    if inc_ordered_components:
        data["ordered_components"] = ordered_components

    return data
