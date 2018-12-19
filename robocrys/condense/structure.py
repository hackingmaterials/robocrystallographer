"""
This module defines a class for condensing structures into dict representations.
"""

from typing import Optional, Dict, Any, List, Tuple

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
                                         get_structure_inequiv_components,
                                         Component)
from robocrys.condense.mineral import MineralMatcher
from robocrys.condense.molecule import MoleculeNamer
from robocrys.condense.site import SiteAnalyzer


class StructureCondenser(object):
    """Class to transform a structure into an intermediate dict representation.

    Args:
        use_conventional_cell: Whether to always use the convention cell
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
    """

    def __init__(self,
                 use_conventional_cell: bool = True,
                 near_neighbors: Optional[NearNeighbors] = None,
                 mineral_matcher: Optional[MineralMatcher] = None,
                 use_symmetry: bool = False,
                 symprec: float = 0.01,
                 simplify_molecules: bool = True,
                 use_iupac_formula: bool = True,
                 use_common_formulas: bool = True):
        if not near_neighbors:
            near_neighbors = CrystalNN()

        if not mineral_matcher:
            mineral_matcher = MineralMatcher()

        self.use_conventional_cell = use_conventional_cell
        self.near_neighbors = near_neighbors
        self.mineral_matcher = mineral_matcher
        self.use_symmetry = use_symmetry
        self.symprec = symprec
        self.simplify_molecules = simplify_molecules
        self.use_common_formulas = use_common_formulas
        self.use_iupac_formula = use_iupac_formula

    def condense_structure(self, structure: Structure) -> Dict[str, Any]:
        """Condenses the structure into an intermediate dict representation.

        Args:
            structure: A pymatgen structure object.

        Returns:
            The condensed structure information. The data is formatted as a
            :obj:`dict` with a fixed set of keys. For example, the condensed
            representation of MoS2 is shown below::

                {'mineral': {'type': 'Molybdenite',
                             'distance': -1,
                             'n_species_type_match': True,
                             'simplified': False},
                 'formula': 'MoS2',
                 'spg_symbol': 'P6_3/mmc',
                 'crystal_system': 'hexagonal',
                 'dimensionality': 2,
                 'sites': {0: {'element': 'Mo4+',
                               'geometry': {'likeness': 0.5566,
                                            'type': 'pentagonal pyramidal'},
                               'nn': [2, 2, 2, 2, 2, 2],
                               'nnn': {'edge': [0, 0, 0, 0, 0, 0]},
                               'poly_formula': 'S6',
                               'sym_labels': (1,)},
                           2: {'element': 'S2-',
                               'geometry': {'likeness': 0.8888,
                                   'type': '3-coordinate'},
                               'nn': [0, 0, 0],
                               'nnn': {'corner': [2, 2, 2, ...],
                                       'face': [2]},
                               'poly_formula': None,
                               'sym_labels': (1,)}},
                 'distances': {0: {2: [2.42, 2.42, 2.42, ...]},
                               2: {0: [2.42, 2.42, 2.42]}},
                 'angles': {0: {0: {'edge': [82.60, 82.60, 82.60, ...]}},
                            2: {2: {'corner': [135.20, 82.60, 135.20, ...],
                                    'face': [80.70, 80.70, 80.70]}}},
                 'components': {0: {'formula':
                                    'MoS2',
                                    'sites': [0, 2]
                                    'dimensionality': 2,
                                    'molecule_name': None,
                                    'orientation': (0, 0, 1)}},
                 'component_makeup': [0, 0],
                 'vdw_heterostructure_info': {'intercalant_formulas': [],
                                              'num_repetitions': 2,
                                              'repeating_unit': ['MoS2']}}

            For more information on the ``mineral`` key, see
            :meth:`StructureCondenser._condense_mineral`. For more information
            on the ``site`` key see :meth:`SiteAnalyzer.get_all_site_summaries`.
            For more information on the ``distance`` and ``angles`` keys, see
            :meth:`SiteAnalyzer.get_all_bond_distance_summaries` and
            :meth:`SiteAnalyzer.get_all_connectivity_angle_summaries`. For more
            information on the ``components`` and ``components_makeup`` keys,
            see :meth:`StructureCondenser._condense_components`.
        """
        # sort so we can give proper symmetry labels
        structure = structure.get_sorted_structure()

        # wrap all site coords into unit cell
        structure.translate_sites(range(structure.num_sites), [1, 1, 1])

        sga = SpacegroupAnalyzer(structure, symprec=self.symprec)
        if self.use_conventional_cell:
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
            hs_info = None

        structure_data['vdw_heterostructure_info'] = hs_info

        return structure_data

    def _condense_mineral(self,
                          structure: Structure,
                          components: List[Component]
                          ) -> Dict[str, Any]:
        """Condenses the mineral data.

        Initially the original structure will be matched against a library
        of known minerals. If no match is found, the structure will be
        reconstructed, with all zero-dimensional components replaced by single
        atoms at their centre of mass, and the matching performed again.
        This allows for matching of systems with molecular substituents in
        known mineral prototypes.

        Args:
            structure: A pymatgen structure object.
            components: A list of structure components, generated using
                :obj:`pymatgen.analysis.dimensionality.get_structure_components`
                with ``inc_molecule_graph=True``.

        Returns:
            The mineral information as a :obj:`dict`, formatted as::

                {
                    'type': mineral (str),
                    'distance': distance (float),
                    'n_species_type_match': match (bool),
                    'simplified': simplified (bool
                }

            The  ``type``, ``distance``, and ``n_species_types_match`` keys
            correspond to the mineral name (e.g. Perovskite), the fingerprint
            distance between the prototype and known mineral, and whether the
            number of species types in the structure matches the number in the
            known prototype. If no mineral match is determined, the mineral type
            will be ``None``. If an exact mineral match is found the distance
            will be set to ``-1``.
        """

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

    def _condense_formula(self,
                          structure: Structure,
                          components: List[Component]
                          ) -> str:
        """Condenses the structure formula.

        If :attr:`StructureCondenser.use_common_formulas` is ``True`` and the
        formula is found in a database of 100,000 known formulas we
        preferentially use the known formula, otherwise we reconstruct the
        formula from the structure components.

        Args:
            structure: A pymatgen structure object.
            components: A list of structure components, generated using
                :obj:`pymatgen.analysis.dimensionality.get_structure_components`

        Returns:
            The chemical formula.
        """
        reduced_formula = structure.composition.reduced_formula

        if self.use_common_formulas and reduced_formula in common_formulas:
            formula = common_formulas[reduced_formula]
        else:
            formula = get_formula_from_components(
                components, use_common_formulas=self.use_common_formulas)
        return formula

    def _condense_components(self,
                             components: List[Component],
                             spacegroup_analyzer: SpacegroupAnalyzer,
                             site_analyzer: SiteAnalyzer
                             ) -> Tuple[Dict[int, Any], List[int]]:
        """Condenses the component data.

        Args:
            components: A list of structure components, generated using
                :obj:`pymatgen.analysis.dimensionality.get_structure_components`
            site_analyzer: A site analyzer object for the structure containing
                the components.
            spacegroup_analyzer: A space group analyzer object for the structure
                containing the components.

        Returns:
            The condensed component data and the component makeup of the
            structure. The condensed components have the form::

                {
                    0: {
                        'formula': 'MoS2',
                        'sites': [0, 2]
                        'dimensionality': 2,
                        'molecule_name': None,
                        'orientation': (0, 0, 1)
                    }
                }

            Where the ``0`` key is the component identifier. The ``sites``
            key gives a :obj:`list` of the indexes of the inequivalent sites
            in the structure. If the component is zero-dimensional and
            is a known molecule, the ``molecule_name`` key will be a :obj:`str`
            with the molecule name. The ``orientation`` key is the miller
            index (for two-dimensional components) or direction of propagation
            (for one-dimensional components).
        """
        if self.use_symmetry:
            inequiv_components = get_sym_inequiv_components(
                components, spacegroup_analyzer)
        else:
            inequiv_components = get_structure_inequiv_components(components)

        molecule_namer = MoleculeNamer()

        components = {}
        component_makeup = []
        for i, component in enumerate(inequiv_components):
            formula = get_component_formula(
                component, use_iupac_formula=self.use_iupac_formula,
                use_common_formulas=self.use_common_formulas)

            sites = site_analyzer.get_inequivalent_site_indices(
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
