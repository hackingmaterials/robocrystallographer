"""
This module defines a function to turn a structure into a dict representation.
"""
from collections import defaultdict
from typing import Optional, Dict, Text, Any

from pymatgen.analysis.dimensionality import get_structure_components
from pymatgen.analysis.local_env import NearNeighbors, CrystalNN
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from robocrys import common_formulas
from robocrys.component import (get_sym_inequiv_components,
                                get_reconstructed_structure,
                                get_formula_from_components,
                                get_component_formula,
                                components_are_vdw_heterostructure,
                                get_vdw_heterostructure_information)
from robocrys.mineral import MineralMatcher
from robocrys.site import SiteAnalyzer


class StructureCondenser(object):
    """Class to transform a structure into an intermediate dict representation.

    Args:
        near_neighbors: A ``NearNeighbors`` instance used to calculate the
            bonding in the structure. For example, one of
            :class:`pymatgen.analysis.local_env.CrystalNN`,
            :class:`pymatgen.analysis.local_env.VoronoiNN`, etc. Defaults to
            ``None``, in which case
            :class:`pymatgen.analysis.local_env.CrystalNN` will be used.
        mineral_matcher: A ``MineralMatcher`` instance. Defaults to ``None``
            in which case the default ``MineralMatcher`` settings will be used.
        symprec: The tolerance used when determining the symmetry of
            the structure. The symmetry is used to determine if multiple
            sites are symmetrically equivalent.
        use_iupac_formula (bool, optional): Whether to order the
            formula by the iupac "electronegativity" series, defined in
            Table VI of "Nomenclature of Inorganic Chemistry (IUPAC
            Recommendations 2005)". This ordering effectively follows
            the groups and rows of the periodic table, except the
            Lanthanides, Actanides and hydrogen. If set to None, the elements
            will be ordered according to the electronegativity values.
        use_common_formulas: Whether to use the database of common formulas.
            The common formula will be used preferentially to the iupac or
            reduced formula.
    """

    def __init__(self,
                 near_neighbors: Optional[NearNeighbors] = None,
                 mineral_matcher: Optional[MineralMatcher] = None,
                 symprec: float = 0.01,
                 simplify_molecules: bool = True,
                 use_iupac_formula: bool = True,
                 use_common_formulas: bool = True):
        if not near_neighbors:
            near_neighbors = CrystalNN()

        if not mineral_matcher:
            mineral_matcher = MineralMatcher()

        self.near_neighbors = near_neighbors
        self.mineral_matcher = mineral_matcher
        self.symprec = symprec
        self.simplify_molecules = simplify_molecules
        self.use_common_formulas = use_common_formulas
        self.use_iupac_formula = use_iupac_formula

    def condense_structure(self, structure: Structure) -> Dict[Text, Any]:
        """Condenses the structure into a dict representation.

        Args:
            structure: A pymatgen structure object.

        Returns:
            WIP
        """
        # wrap all site coords into unit cell
        structure.translate_sites(range(structure.num_sites), [1, 1, 1])

        sga = SpacegroupAnalyzer(structure, symprec=self.symprec)
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

        hs_info = None
        is_vdw_heterostructure = components_are_vdw_heterostructure(components)
        if is_vdw_heterostructure:
            hs_info = get_vdw_heterostructure_information(
                components, use_iupac_formula=self.use_iupac_formula,
                use_common_formulas=self.use_common_formulas)

        structure_data['is_vdw_heterostructure'] = is_vdw_heterostructure
        structure_data['vdw_heterostructure_info'] = hs_info

        components_data, n_components = self._condense_components(
            components, bonded_structure, sga)
        structure_data['components'] = components_data
        structure_data['n_components'] = n_components

        return structure_data

    def _condense_mineral(self, structure, components):
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
        # if the formula is known from the list of 100,000 known formulae we
        # preferentially use that, else we reconstruct it from the components.
        reduced_formula = structure.composition.reduced_formula
        if reduced_formula in common_formulas:
            formula = common_formulas[reduced_formula]
        else:
            formula = get_formula_from_components(
                components, use_common_formulas=self.use_common_formulas)
        return formula

    def _condense_components(self, components, bonded_structure, sga):
        sym_inequiv_components = get_sym_inequiv_components(components, sga)
        site_analyzer = SiteAnalyzer(bonded_structure, self.symprec)

        # defaultdict of defaultdicts
        cc = defaultdict(lambda: defaultdict(dict))
        total_components = 0
        for component in sym_inequiv_components:
            formula = get_component_formula(
                component, use_iupac_formula=self.use_iupac_formula,
                use_common_formulas=self.use_common_formulas)
            dimen = component['dimensionality']
            count = component['count']

            component_data = {
                'orientation': component['orientation'],
                'count': count,
                'sites': []
            }

            # for each inequivalent site in the component get a site description
            for site_id in component['inequivalent_ids']:
                geometry = site_analyzer.get_site_geometry(site_id)
                nn_info = site_analyzer.get_nearest_neighbor_summary(site_id)
                nnn_info = site_analyzer.get_next_nearest_neighbor_summary(
                    site_id)
                component_data['sites'].append({
                    'element': bonded_structure.structure[site_id].specie.name,
                    'geometry': geometry,
                    'nn_info': nn_info,
                    'next_nn_info': nnn_info
                })

            if dimen in cc and formula in cc[dimen]:
                cc[dimen][formula]['sym_inequiv_components'].append(
                    component_data)
            else:
                cc[dimen][formula]['sym_inequiv_components'] = [component_data]

            if 'count' in cc[dimen][formula]:
                cc[dimen][formula]['count'] += count
            else:
                cc[dimen][formula]['count'] = count

            total_components += count
        return cc, total_components
