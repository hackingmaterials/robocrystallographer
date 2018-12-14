"""
This module defines a class to turn a structure into a dict representation.
"""

from typing import Optional, Dict, Any

from pymatgen.analysis.dimensionality import get_structure_components
from pymatgen.analysis.local_env import NearNeighbors, CrystalNN
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from robocrys import common_formulas
from robocrys.condense.component import (get_sym_inequiv_components,
                                         get_reconstructed_structure,
                                         get_formula_from_components,
                                         get_component_formula,
                                         components_are_vdw_heterostructure,
                                         get_vdw_heterostructure_information,
                                         get_structure_inequiv_components)
from robocrys.condense.mineral import MineralMatcher
from robocrys.condense.molecule import MoleculeNamer
from robocrys.condense.site import SiteAnalyzer


class StructureCondenser(object):
    """Class to transform a structure into an intermediate dict representation.

    Args:
        force_conventional_cell: Whether to always use the convention cell
            representation of the structure.
        near_neighbors: A ``NearNeighbors`` instance used to calculate the
            bonding in the structure. For example, one of
            :class:`pymatgen.analysis.local_env.CrystalNN`,
            :class:`pymatgen.analysis.local_env.VoronoiNN`, etc. Defaults to
            ``None``, in which case
            :class:`pymatgen.analysis.local_env.CrystalNN` will be used.
        mineral_matcher: A ``MineralMatcher`` instance. Defaults to ``None``
            in which case the default ``MineralMatcher`` settings will be used.
        use_symmetry: Whether to use symmetry to determine if structure
            compoents and sites are equivalent. If ``False``, the site geometry
            and bonding graph information will be used.
        symprec: The tolerance used when determining the symmetry of
            the structure. The symmetry is used to determine if multiple
            sites are symmetrically equivalent. If ``use_symmetry=False`` this
            option will be ignored.
        use_iupac_formula (bool, optional): Whether to order formulas
            by the iupac "electronegativity" series, defined in
            Table VI of "Nomenclature of Inorganic Chemistry (IUPAC
            Recommendations 2005)". This ordering effectively follows
            the groups and rows of the periodic table, except the
            Lanthanides, Actanides and hydrogen. If set to ``False``, the
            elements will be ordered according to the electronegativity values.
        use_common_formulas: Whether to use the database of common formulas.
            The common formula will be used preferentially to the iupac or
            reduced formula.
        group_bond_lengths: Whether to group together sites with different
            nearest neighbor bond lengths but otherwise comparable geometries
            and (next) nearest neighbor information. In some cases this will
            significantly simplify the condensed represenation.
        store_oxidation_states: Whether to store the oxidation states of the
            sites.
    """

    def __init__(self,
                 force_conventional_cell: bool = False,
                 near_neighbors: Optional[NearNeighbors] = None,
                 mineral_matcher: Optional[MineralMatcher] = None,
                 use_symmetry: bool = False,
                 symprec: float = 0.01,
                 simplify_molecules: bool = True,
                 use_iupac_formula: bool = True,
                 use_common_formulas: bool = True,
                 group_bond_lengths: bool = True,
                 store_oxidation_states: bool = True):
        if not near_neighbors:
            near_neighbors = CrystalNN()

        if not mineral_matcher:
            mineral_matcher = MineralMatcher()

        self.force_conventional_cell = force_conventional_cell
        self.near_neighbors = near_neighbors
        self.mineral_matcher = mineral_matcher
        self.use_symmetry = use_symmetry
        self.symprec = symprec
        self.simplify_molecules = simplify_molecules
        self.use_common_formulas = use_common_formulas
        self.use_iupac_formula = use_iupac_formula
        self.group_bond_lengths = group_bond_lengths
        self.store_oxidation_states = store_oxidation_states

    def condense_structure(self, structure: Structure) -> Dict[str, Any]:
        """Condenses the structure into a dict representation.

        Args:
            structure: A pymatgen structure object.

        Returns:
            WIP
        """
        # sort so we can give proper symmetry labels
        structure = structure.get_sorted_structure()

        # wrap all site coords into unit cell
        structure.translate_sites(range(structure.num_sites), [1, 1, 1])

        sga = SpacegroupAnalyzer(structure, symprec=self.symprec)
        if self.force_conventional_cell:
            structure = sga.get_conventional_standard_structure()
        else:
            structure = sga.get_symmetrized_structure()

        bonded_structure = self.near_neighbors.get_bonded_structure(structure)

        components = get_structure_components(
            bonded_structure, inc_orientation=True, inc_site_ids=True,
            inc_molecule_graph=True)

        dimensionality = max(c['dimensionality'] for c in components)
        mineral = self._condense_mineral(structure, components)
        formula = self._condense_formula(structure, components)

        structure_data = {
            'formula': formula,
            'spg_symbol': sga.get_space_group_symbol(),
            'crystal_system': sga.get_crystal_system(),
            'mineral': mineral,
            'dimensionality': dimensionality,
        }

        site_analyzer = SiteAnalyzer(
            bonded_structure, use_symmetry_equivalent_sites=self.use_symmetry,
            symprec=self.symprec)
        structure_data['sites'] = site_analyzer.get_all_site_summaries()
        structure_data['distances'] = site_analyzer. \
            get_all_bond_distance_summaries()
        structure_data['angles'] = site_analyzer. \
            get_all_connectivity_angle_summaries()

        component_summary, component_makeup = self._condense_components(
            components, sga, site_analyzer)
        structure_data['components'] = component_summary
        structure_data['component_makeup'] = component_makeup

        if components_are_vdw_heterostructure(components):
            hs_info = get_vdw_heterostructure_information(
                components, use_iupac_formula=self.use_iupac_formula,
                use_common_formulas=self.use_common_formulas)
        else:
            hs_info = False

        structure_data['vdw_heterostructure_info'] = hs_info

        return structure_data

    def _condense_mineral(self, structure, components):
        """Utility function to condense the mineral data."""

        mineral = self.mineral_matcher.get_best_mineral_name(structure)

        if not mineral['type']:
            mineral_structure = get_reconstructed_structure(
                components, simplify_molecules=self.simplify_molecules)
            mineral = self.mineral_matcher.get_best_mineral_name(
                mineral_structure)
            mineral['simplified'] = True
        else:
            mineral['simplified'] = False

        return mineral

    def _condense_formula(self, structure, components):
        """Utility function to condense the formula data."""
        # if the formula is known from the list of 100,000 known formulae we
        # preferentially use that, else we reconstruct it from the components.
        reduced_formula = structure.composition.reduced_formula
        if reduced_formula in common_formulas:
            formula = common_formulas[reduced_formula]
        else:
            formula = get_formula_from_components(
                components, use_common_formulas=self.use_common_formulas)
        return formula

    def _condense_components(self, components, sga, site_analyzer):
        """Utility function to condense the component data."""
        if self.use_symmetry:
            inequiv_components = get_sym_inequiv_components(components, sga)
        else:
            inequiv_components = get_structure_inequiv_components(components)

        molecule_namer = MoleculeNamer()

        components = {}
        component_makeup = []

        for i, component in inequiv_components:
            formula = get_component_formula(
                component, use_iupac_formula=self.use_iupac_formula,
                use_common_formulas=self.use_common_formulas)

            sites = site_analyzer.get_inequivalent_site_ids(
                component['site_ids'])

            if component['dimensionality'] == 0:
                molecule_name = molecule_namer.get_name_from_molecule_graph(
                    component['molecule_graph'])
            else:
                molecule_name = None

            components[i] = {'formula': formula,
                             'dimensionality': component['dimensionality'],
                             'orientation': component['orientation'],
                             'molecule_name': molecule_name,
                             'sites': sites}
            component_makeup.extend([i] * component['count'])

        return components, component_makeup
