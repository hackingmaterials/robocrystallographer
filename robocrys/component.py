"""
This module implements functions for handling structure components.
"""
from copy import deepcopy

import numpy as np

from monty.fractions import gcd
from typing import List, Dict, Text, Any, Tuple

from pymatgen.core.structure import Structure, PeriodicSite
from pymatgen.core.periodic_table import get_el_sp
from pymatgen.core.composition import Composition, iupac_ordering_dict
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.string import formula_double_format

from robocrys import formula_mapping

Component = Dict[Text, Any]


def get_sym_inequiv_components(components: List[Component],
                               spg_analyzer: SpacegroupAnalyzer
                               ) -> List[Component]:
    """Gets and counts the symmetrically inequivalent components.

    Component data has to have been generated with ``inc_site_ids=True``.

    Args:
        components: A list of structure components, generated using
            :obj:`pymatgen.analysis.dimensionality.get_structure_components`,
            with ``inc_site_ids=True``.
        spg_analyzer: A spacegroup analyzer object for the structure containing
            the components.

    Returns:
        A list of the symmetrically inequivalent components. Any duplicate
        components will only be returned once. The component objects are in the
        same format is given by
        :obj:`pymatgen.analysis.dimensionality.get_structure_components` but
        have two additional fields:

        - ``"count"`` (:obj:`int`): The number of times this component appears
            in the structure.
        - ``"inequivalent_ids"`` (``list[int]``): The site indices of the
            symmetry inequivalent atoms in the component.
    """
    components = deepcopy(components)
    sym_inequiv_components = {}
    equivalent_atoms = spg_analyzer.get_symmetry_dataset()['equivalent_atoms']

    for component in components:
        sym_ids = frozenset(
            equivalent_atoms[x] for x in component['site_ids'])

        # if two components are composed of atoms that are symmetrically
        # equivalent they are the same.
        if sym_ids in sym_inequiv_components:
            sym_inequiv_components[sym_ids]['count'] += 1
            continue

        component['count'] = 1
        component['inequivalent_ids'] = tuple(sym_ids)
        sym_inequiv_components[sym_ids] = component

    return list(sym_inequiv_components.values())


def get_formula_inequiv_components(components: List[Component],
                                   use_iupac_formula: bool=True,
                                   ) -> List[Component]:
    """Gets and counts the inequivalent components based on their formuula.

    Note that the counting of compounds is different to in
    ``get_sym_inequiv_equivalent``. I.e. the count is not the number of
    components with the same formula. For example, the count of the formula
    "GaAs" in a system with two Ga2As2 components would be 4.

    Args:
        components: A list of structure components, generated using
            :obj:`pymatgen.analysis.dimensionality.get_structure_components`,
            with ``inc_site_ids=True``.
        use_iupac_formula (bool, optional): Whether to order the
            formula by the iupac "electronegativity" series, defined in
            Table VI of "Nomenclature of Inorganic Chemistry (IUPAC
            Recommendations 2005)". This ordering effectively follows
            the groups and rows of the periodic table, except the
            Lanthanides, Actanides and hydrogen. If set to None, the elements
            will be ordered according to the electronegativity values.

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
    inequiv_components = {}

    for component in components:
        formula, factor = component['structure'].composition. \
            get_reduced_formula_and_factor(iupac_ordering=use_iupac_formula)

        # if the formula is known from the list of 100,000 known formulae we
        # preferentially use that
        reduced_formula = component['structure'].composition.reduced_formula
        if reduced_formula in formula_mapping:
            formula = formula_mapping[reduced_formula]

        # if two components have the same composition we treat them as the same
        if formula in inequiv_components:
            inequiv_components[formula]['count'] += factor
            continue

        component['count'] = factor
        component['formula'] = formula

        inequiv_components[formula] = component

    return list(inequiv_components.values())


def filter_molecular_components(components: List[Component]
                                ) -> Tuple[List[Component], List[Component]]:
    """Separate a list of components into molecular and non-molecular components.

    Args:
        components: A list of structure components, generated using
            :obj:`pymatgen.analysis.dimensionality.get_structure_components`.

    Returns:
        The filtered components as a tuple of ``(molecular_components,
        other_components)``.
    """
    molecular_components = [c for c in components if c['dimensionality'] == 0]
    other_components = [c for c in components if c['dimensionality'] != 0]

    return molecular_components, other_components


def get_reconstructed_structure(components: List[Component],
                                simplify_molecules: bool=True
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
            lattice = mol_components[0]['structure'].lattice

            mol_sites = [
                PeriodicSite(c['structure'][0].specie,
                             c['molecule_graph'].molecule.center_of_mass,
                             lattice, coords_are_cartesian=True)
                for c in mol_components]

    if components:
        other_sites = [site for c in components for site in c['structure']]

    return Structure.from_sites(other_sites + mol_sites)


def get_formula_from_components(components: List[Component],
                                molecules_first: bool=False,
                                use_iupac_formula: bool=True) -> Text:
    """Reconstructs a chemical formula from structure components.

    The chemical formulas for the individual components will be grouped
    together. If two components share the same composition, they will be
    treated as equivalent.

    Args:
        components: A list of structure components, generated using
            :obj:`pymatgen.analysis.dimensionality.get_structure_components`.
        molecules_first: Whether to put any molecules (zero-dimensional
            components) at the beginning of the formula.
        use_iupac_formula (bool, optional): Whether to order the
            formula by the iupac "electronegativity" series, defined in
            Table VI of "Nomenclature of Inorganic Chemistry (IUPAC
            Recommendations 2005)". This ordering effectively follows
            the groups and rows of the periodic table, except the
            Lanthanides, Actanides and hydrogen. If set to None, the elements
            will be ordered according to the electronegativity values.

    Returns:
        The chemical formula.
    """
    def order(comp_formula):
        composition = Composition(comp_formula)
        if use_iupac_formula:
            return (sum([iupac_ordering_dict[get_el_sp(s)]
                         for s in composition.elements]) /
                    len(composition.elements))
        else:
            return composition.average_electroneg

    components = get_formula_inequiv_components(
        components, use_iupac_formula=use_iupac_formula)

    if molecules_first:
        mol_comps, other_comps = filter_molecular_components(components)
    else:
        mol_comps = []
        other_comps = components

    formulas = (sorted([c['formula'] for c in mol_comps], key=order) +
                sorted([c['formula'] for c in other_comps], key=order))

    # if components include special formulas, then the count can be 0.5
    # therefore if any non integer amounts we can just use a factor of 2
    all_int = all(v['count'] % 1 == 0 for v in components)
    prefactor = 1 if all_int else 2
    form_count_dict = {c['formula']: int(c['count'] * prefactor)
                       for c in components}

    # the following is based on ``pymatgen.core.composition.reduce_formula``
    num_comps = len(formulas)
    factor = abs(gcd(*(form_count_dict.values())))

    reduced_form = []
    for i in range(0, num_comps):
        formula = formulas[i]
        normamt = form_count_dict[formula] * 1.0 / factor
        formatted_formula = formula if normamt == 1 else "({})".format(formula)
        reduced_form.append(formatted_formula)
        reduced_form.append(formula_double_format(normamt))

    reduced_form = "".join(reduced_form)
    return reduced_form


def components_are_vdw_heterostructure(components: List[Component]
                                       ) -> bool:
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

    if len([c for c in components if c['dimensionality'] == 2]):
        return True
    else:
        return False


def get_vdw_heterostructure_information(components: List[Component],
                                        use_iupac_formula: bool=True,
                                        ) -> Dict[Text, Any]:
    """Gets information about ordering of components in a vdw heterostructure.

    Args:
        components: A list of structure components, generated using
            :obj:`pymatgen.analysis.dimensionality.get_structure_components`
            with ``inc_orientation=True``.
        use_iupac_formula (bool, optional): Whether to order the
            formula by the iupac "electronegativity" series, defined in
            Table VI of "Nomenclature of Inorganic Chemistry (IUPAC
            Recommendations 2005)". This ordering effectively follows
            the groups and rows of the periodic table, except the
            Lanthanides, Actanides and hydrogen. If set to None, the elements
            will be ordered according to the electronegativity values.

    Returns:
        Information on the heterostructure, as an :obj:`dict` with they keys:

        - ``"ordered_components"`` (``list[component]``): A :obj:`List` of
            components, ordered as they appear in the heteostructure stacking
            direction.
        - ``"repeating_unit"`` (``list[str]``: A :obj:`List` of formulas of the
            smallest repeating series of components. For example. if the
            structure consists of A and B components ordered as "A B A B A B",
            the repeating unit is "A B".
        - ``"num_repetitions"`` (:obj:`int`): The number of repetitions of the
            repeating unit that forms the overall structure. For example. if
            the structure consists of A and B components ordered as
            "A B A B A B", the number of repetitions is 3.
        - ``"intercalants"`` (``list[component]``: A :obj:`List` of intercalated
            components.
    """
    if not components_are_vdw_heterostructure(components):
        raise ValueError("Components do not form a heterostructure.")

    try:
        millers = set([c['orientation'] for c in components
                       if c['dimensionality'] == 2])
    except KeyError as e:
        if "orientation" in str(e):
            raise KeyError("Components not generated with inc_orientation=True")
        else:
            raise e

    if len(millers) != 1:
        raise ValueError("2D components don't all have the same orientation.")

    cart_miller = components[0]['structure'].lattice.get_cartesian_coords(
        millers.pop()).tolist()

    # plane is used to find the distances of all components along a certain axis
    # should use normal vector of plane to get exact distances but we just care
    # about relative distances.
    def distances_to_plane(points):
        return [np.dot(cart_miller, pp) for pp in points]

    min_distances = [min(distances_to_plane(c['structure'].cart_coords))
                     for c in components]

    # sort the components by distance to plane
    ordering = np.argsort(min_distances)
    ordered_components = [components[x] for x in ordering]

    # only consider the layered components formulae
    ordered_layers = [c for c in ordered_components if c['dimensionality'] == 2]
    ordered_layers_formula = [
        c['structure'].composition.get_reduced_formula_and_factor(
            use_iupac_formula)[0] for c in ordered_layers]
    num_layer_formulas = len(set(ordered_layers_formula))

    repeating_formula = ordered_layers_formula
    num_repetitions = 1

    # depending on the number of inequivalent formulae, there is a maximum
    # number of repetitions that can occur. To avoid unnecessary work we start
    # from this number of repetitions and move to 1 repetition (e.g. no
    # repetition)
    max_repetitions = int(np.floor(len(ordered_layers)/num_layer_formulas))
    for n in range(max_repetitions, 0, -1):
        if (all([len(set(ordered_layers_formula[i::num_layer_formulas])) == 1
                 for i in range(n)]) and len(ordered_layers) % n == 0):
            repeating_formula = ordered_layers_formula[
                                :int(len(ordered_layers)/n)]
            num_repetitions = n
            break

    intercalants = [c for c in components if c['dimensionality'] < 2]
    data = {"ordered_components": ordered_components,
            'intercalants': intercalants,
            'repeating_unit': repeating_formula,
            'num_repetitions': num_repetitions}
    return data
