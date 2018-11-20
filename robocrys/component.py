"""
This module implements functions for handling structure components.
"""

from typing import List, Dict, Text, Any, Tuple

from pymatgen import Structure, PeriodicSite
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

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
        A list of the symmetrically inequivalent spacegroups. Any duplicate
        components will only be returned once. The component objects are in the
        same format is given by
        :obj:`pymatgen.analysis.dimensionality.get_structure_components` but
        have two additional fields:

        - "count" (int): The number of times this component appears in the
          structure.
        - "inequivalent_ids" (list[int]): The site indices of the symmetry
          inequivalent atoms in the component.
    """
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
