"""
This module defines a function to turn a structure into a dict representation.
"""
import copy
from collections import defaultdict
from typing import Optional, Dict, Text, Any, List

from pymatgen.analysis.dimensionality import get_structure_components
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import NearNeighbors, CrystalNN
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from robocrys import common_formulas
from robocrys.component import (get_sym_inequiv_components,
                                get_reconstructed_structure,
                                get_formula_from_components,
                                get_component_formula,
                                components_are_vdw_heterostructure,
                                get_vdw_heterostructure_information,
                                get_structure_inequiv_components, Component)
from robocrys.mineral import MineralMatcher
from robocrys.molecule import MoleculeNamer
from robocrys.site import SiteAnalyzer, geometries_match, nn_summaries_match, \
    nnn_summaries_match


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
        group_bond_lengths: Whether to group together sites with different
            nearest neighbor bond lengths but otherwise comparable geometries
            and (next) nearest neighbor information. In some cases this will
            significantly simplify the condensed represenation.
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
                 group_bond_lengths: bool =True):
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

    def _condense_components(self, components, bonded_structure, sga):
        """Utility function to condense the component data."""
        if self.use_symmetry:
            inequiv_components = get_sym_inequiv_components(components, sga)
        else:
            inequiv_components = get_structure_inequiv_components(components)

        site_analyzer = SiteAnalyzer(
            bonded_structure, use_symmetry=self.use_symmetry,
            symprec=self.symprec)

        molecule_namer = MoleculeNamer()

        # defaultdict of defaultdicts
        cc = defaultdict(lambda: defaultdict(dict))
        total_components = 0
        for component in inequiv_components:

            formula = get_component_formula(
                component, use_iupac_formula=self.use_iupac_formula,
                use_common_formulas=self.use_common_formulas)
            dimen = component['dimensionality']
            count = component['count']

            component_data = {
                'orientation': component['orientation'],
                'count': count,
            }

            sites = self._condense_component_site_data(
                bonded_structure, site_analyzer, component)
            component_data['sites'] = sites

            if dimen in cc and formula in cc[dimen]:
                cc[dimen][formula]['inequiv_components'].append(
                    component_data)
            else:
                cc[dimen][formula]['inequiv_components'] = [component_data]

            if 'count' in cc[dimen][formula]:
                cc[dimen][formula]['count'] += count
            else:
                cc[dimen][formula]['count'] = count

            if dimen == 0:
                molecule_name = molecule_namer.get_name_from_molecule_graph(
                    component['molecule_graph'])

                if molecule_name:
                    cc[dimen][formula]['molecule_name'] = molecule_name

            total_components += count
        return cc, total_components

    def _condense_component_site_data(self,
                                      bonded_structure: StructureGraph,
                                      site_analyzer: SiteAnalyzer,
                                      component: Component
                                      ) -> List[Dict[str, Any]]:
        """Utility function to condense component site data."""

        def site_order(index):
            el = bonded_structure.structure[index].specie.element
            return el.iupac_ordering if self.use_iupac_formula else el.X

        inequivalent_ids = site_analyzer.get_inequivalent_site_ids(
            component['site_ids'])
        inequivalent_ids = sorted(inequivalent_ids, key=site_order)

        sites = []
        # for each inequivalent site in the component get a site description
        for site_id in inequivalent_ids:
            geometry = site_analyzer.get_site_geometry(site_id)
            nn_data = site_analyzer.get_nearest_neighbor_summary(site_id)
            nnn_data = site_analyzer.get_next_nearest_neighbor_summary(
                site_id)
            sites.append({
                'element': bonded_structure.structure[site_id].specie.name,
                'geometry': geometry,
                'nn_data': nn_data,
                'nnn_data': nnn_data
            })

        if self.group_bond_lengths:
            # merge sites with same geometry, NN info and NNN info but with
            # different NN bond lengths and NNN angles.
            sites = merge_similar_sites(sites)

        return sites


def merge_similar_sites(sites: List[Dict[Text, Any]]):
    """Merges sites with the same properties except bond angles and distances.

    Args:
        sites: A list of sites. Each site is formatted as a :ob:`dict` with the
            keys:

            - ``'element'`` (``str``): The element of the site.
            - ``'geometry'`` (``dict``): The geometry, as output by
                :meth:`SiteAnalyzer.get_site_geometry`.
            - ``'nn_data'`` (``dict``): The nearest neighbor data, as output by
                :meth:`SiteAnalyzer.get_nearest_neighbor_summary`.
            - ``'nnn_data'`` (``dict``): The next nearest neighbor data, as
                given by :meth:`SiteAnalyzer.get_next_nearest_neighbor_summary`.

    Returns:
        A list of merged sites with the same format as above. Merged sites
        have a different ``nn_data`` format than unmerged sites. For example,
        ``nn_data`` in unmerged sites is formatted as::

            {
                'Sn': {
                    'n_sites': 6,
                    'inequiv_groups': [
                        {
                            'n_sites': 4,
                            'inequiv_id': 0,
                            'dists': [1, 1, 2, 2]
                        },
                        {
                            'n_sites': 2,
                            'inequiv_id': 1,
                            'dists': [3, 3]
                        }
                    ]
                }
            }

        Merged sites do not contain an ``inequiv_groups`` key and are instead
        formatted as::

            {
                'n_sites': 6
                'dists': [1, 1, 1, 2, 2, 2, 2, 3, 3]
                )
            }

        Note that there are now more distances than there are number of sites.
        This is because n_sites gives the number of bonds to a specific site,
        whereas the distances are for the complete set of distances for all
        similar (merged) sites. Similarly, merged next nearest neighbor
        data can contain more angles than number of sites, however, the
        general format of the ``nnn_data`` dict is unaltered.
    """
    sites = copy.deepcopy(sites)
    new_sites = []

    for site in sites:

        matched = False
        for new_site in new_sites:
            elem_match = site['element'] == new_site['element']
            geom_match = geometries_match(
                site['geometry'], new_site['geometry'], likeness_tol=1)
            nn_match = nn_summaries_match(
                site['nn_data'], new_site['nn_data'],
                match_bond_dists=False)
            nnn_match = nnn_summaries_match(
                site['nnn_data'], new_site['nnn_data'], match_bond_angles=False)

            if elem_match and geom_match and nn_match and nnn_match:
                new_site['nn_data'] = _merge_nn_data(site['nn_data'],
                                                     new_site['nn_data'])
                new_site['nnn_data'] = _merge_nnn_data(site['nnn_data'],
                                                       new_site['nnn_data'])
                matched = True
                break

        if not matched:
            # no matches therefore store original site id
            new_sites.append(site)

    return new_sites


def _merge_nn_data(site_nn_data, new_site_nn_data):
    """Utility function to merge nearest neighbor data.

    See the ``merge_similar_sites`` docstring for information on the format of
    the merged data.

    Note an error will be thrown if this function is called on two sites that do
    not have matching nearest neighbor summaries (ignoring bond distances).
    """

    for el in site_nn_data:
        site_dists = [dist for group in
                      site_nn_data[el]['inequiv_groups']
                      for dist in group['dists']]

        if 'inequiv_groups' in new_site_nn_data[el]:
            # remove inequiv_groups key and group all distances
            # together
            groups = new_site_nn_data[el].pop('inequiv_groups')
            dists = [dist for dist_set in groups
                     for dist in dist_set['dists']]
            new_site_nn_data[el]['dists'] = dists + site_dists
        else:
            new_site_nn_data[el]['dists'] += site_dists

    return new_site_nn_data


def _merge_nnn_data(site_nnn_data, new_site_nnn_data):
    """Utility function to merge next nearest neighbor data.

    See the ``merge_similar_sites`` docstring for information on the format of
    the merged data.

    Note an error will be thrown if this function is called on two sites that do
    not have matching next nearest neighbor summaries (ignoring bond angles).
    """
    for el in site_nnn_data:
        for geometry in site_nnn_data[el]:
            for connectivity in site_nnn_data[el][geometry]:
                new_site_nnn_data[el][geometry][connectivity]['angles'].extend(
                    site_nnn_data[el][geometry][connectivity]['angles'])

    return new_site_nnn_data
