"""
This module defines a function to turn a structure into a dict representation.
"""

from typing import Optional, Dict, Text, Any

from pymatgen.analysis.dimensionality import get_structure_component_info
from pymatgen.analysis.local_env import NearNeighbors, CrystalNN
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from robocrys import MineralMatcher, SiteAnalyzer
from robocrys.component import get_sym_inequiv_components


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
    """

    def __init__(self,
                 near_neighbors: Optional[NearNeighbors]=None,
                 mineral_matcher: Optional[MineralMatcher]=None,
                 symprec: float=0.01):
        if not near_neighbors:
            near_neighbors = CrystalNN()

        if not mineral_matcher:
            mineral_matcher = MineralMatcher()

        self.near_neighbors = near_neighbors
        self.mineral_matcher = mineral_matcher
        self.symprec = symprec

    def condense_structure(self, structure: Structure) -> Dict[Text, Any]:
        """Condenses the structure into a dict representation.

        Args:
            structure: A pymatgen structure object.

        Returns:

        """
        # wrap all site coords into unit cell
        structure.translate_sites(range(structure.num_sites), [1, 1, 1])

        sga = SpacegroupAnalyzer(structure, symprec=self.symprec)
        structure = sga.get_symmetrized_structure()

        bonded_structure = self.near_neighbors.get_bonded_structure(structure)

        mineral = self.mineral_matcher.get_best_mineral_name(
            bonded_structure.structure)

        components = get_structure_component_info(
            bonded_structure, inc_orientation=True, inc_site_ids=True)

        dimensionality = max(c['dimensionality'] for c in components)

        structure_data = {
            'formula': structure.composition.reduced_formula,
            'spg': sga.get_space_group_symbol(),
            'mineral': mineral,
            'dimensionality': dimensionality,
            'components': []
        }

        sym_inequiv_components = get_sym_inequiv_components(
            sga.get_symmetry_dataset()['equivalent_atoms'], components)

        site_analyzer = SiteAnalyzer(bonded_structure, self.symprec)

        for component in sym_inequiv_components:
            component_data = {
                'dimensionality': component['dimensionality'],
                'orientation': component['orientation'],
                'count': component['count'],
                'formula': component['structure'].composition.reduced_formula,
                'sites': []
            }

            # for each inequivalent site in the component get a site description
            for site_id in component['inequivalent_ids']:
                geometry = site_analyzer.get_site_geometry(site_id)
                nn_info = site_analyzer.get_nearest_neighbor_summary(site_id)
                nnn_info = site_analyzer.get_next_nearest_neighbor_summary(
                    site_id)
                component_data['sites'].append({
                    'element': structure[site_id].specie.name,
                    'geometry': geometry,
                    'nn_info': nn_info,
                    'next_nn_info': nnn_info
                })

            structure_data['components'].append(component_data)

        return structure_data
