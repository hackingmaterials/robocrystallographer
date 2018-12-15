"""
This module implements a class to resolve the symbolic references in condensed
structure data.
"""

from typing import Dict, Any, List, Union


class DescriptionAdapter(object):
    """Class to facilitate pulling data from the condensed structure dictionary.

    Args:
        condensed_structure: The condensed structure data, formatted as produced
            by :meth:`robocrys.condense.StructureCondenser.condense_structure`.
    """

    def __init__(self, condensed_structure: Dict[str, Any]):
        self._condensed_structure = condensed_structure

    def get_nearest_neighbor_details(self, site_index: int) -> Dict[str, Any]:
        """Gets a summary of all the nearest neighbors to a site.

        Args:
            site_index: The site index (zero based).

        Returns:
            A summary of the nearest neighbor information formatted as::

                {
                    'Sn': {
                        'count': 6,
                        'groups': [
                            {
                                'count': 4,
                                'sym_label': '(0)',
                                'dists': [1, 1, 2, 2]
                            },
                            {
                                'count': 2,
                                'sym_label': '(1)',
                                'dists': [3, 3]
                            }
                        ]
                    }
                }
        """
        nn_sites = self.sites[site_index]['nn']
        nn_elements = [self.sites[nn_site]['element']
                       for nn_site in set(nn_sites)]
        nn_details = {el: {'groups': [], 'count': 0} for el in nn_elements}

        for nn_site in set(nn_sites):
            count = nn_sites.count(nn_site)
            element = self.sites[nn_site]['element']

            # convert the sym_labels tuple into a str. E.g. (2, ) -> "(2)"
            label = "({})".format(
                ",".join(map(str, self.sites[nn_site]['sym_labels'])))

            group = {'count': count,
                     'sym_label': label,
                     'distances': self.distances[site_index][nn_site]}
            nn_details[element]['groups'].append(group)
            nn_details[element]['count'] += count

        return nn_details

    @property
    def mineral(self) -> Dict[str, Union[str, int, bool]]:
        """The mineral data.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure['mineral']

    @property
    def formula(self) -> str:
        """The structure formula.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure['formula']

    @property
    def spg_symbol(self) -> str:
        """The space group symbol.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure['spg_symbol']

    @property
    def crystal_system(self) -> str:
        """The crystal system.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure['crystal_system']

    @property
    def dimensionality(self) -> int:
        """The overall dimensionality.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure['dimensionality']

    @property
    def sites(self) -> Dict[int, Dict[str, Any]]:
        """The site data.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure['sites']

    @property
    def distances(self) -> Dict[int, Dict[int, List[float]]]:
        """The distance data.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure['distances']

    @property
    def angles(self) -> Dict[int, Dict[int, Dict[str, List[float]]]]:
        """The angle data.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure['angles']

    @property
    def components(self) -> Dict[int, Dict[str, Any]]:
        """The component data.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure['components']

    @property
    def component_makeup(self) -> List[int]:
        """The component makeup of the structure.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure['component_makeup']


# def get_next_nearest_neighbor_data(self, site_index: int
#                                    ) -> Dict[str, Any]:
#     """Gets a summary of the next nearest neighbor connectivity.
#
#     Args:
#         site_index: The site index (zero based).
#
#     Returns:
#         A summary of the next nearest neighbor information as a dict.
#         Formatted as::
#
#             {
#                 'Sn': {
#                     'octahedral': {
#                         'corner-sharing': {
#                             'n_sites': 8,
#                             'angles': [180, 180, 180, ...]
#                         }
#                     }
#                 }
#             }
#
#     """
#     nnn_info = self.get_next_nearest_neighbors(site_index)
#
#     # group next nearest neighbors by element, connectivity and geometry
#     # e.g. grouped_nnn looks like {el: {connectivity: {geometry: [sites]}}}
#     grouped_nnn = defaultdict(
#         lambda: defaultdict(lambda: defaultdict(list)))
#
#     for site in nnn_info:
#         grouped_nnn[site['element']][
#             site['geometry']['type']][site['connectivity']].append(site)
#
#     nnn_data = {}
#     for element, geom_data in grouped_nnn.items():
#         nnn_el_data = {}
#         for geometry, con_data in geom_data.items():
#             nnn_geom_data = {}
#             for connectivity, sites in con_data.items():
#                 nnn_geom_data[connectivity] = {
#                     'n_sites': len(sites),
#                     'angles': [angle for site in sites
#                                for angle in site['angles']]}
#             nnn_el_data[geometry] = nnn_geom_data
#         nnn_data[element] = nnn_el_data
#     return nnn_data

# def merge_similar_sites(sites: List[Dict[str, Any]]):
#     """Merges sites with the same properties except bond angles and distances.
#
#     Args:
#         sites: A list of sites. Each site is formatted as a :ob:`dict` with
#           the keys:
#
#             - ``'element'`` (``str``): The element of the site.
#             - ``'geometry'`` (``dict``): The geometry, as output by
#                 :meth:`SiteAnalyzer.get_site_geometry`.
#             - ``'nn_data'`` (``dict``): The nearest neighbor data, as output
#             by :meth:`SiteAnalyzer.get_nearest_neighbor_data`.
#             - ``'nnn_data'`` (``dict``): The next nearest neighbor data, as
#                 given by :meth:`SiteAnalyzer.get_next_nearest_neighbor_data`.
#
#     Returns:
#         A list of merged sites with the same format as above. Merged sites
#         have a different ``nn_data`` format than unmerged sites. For example,
#         ``nn_data`` in unmerged sites is formatted as::
#
#             {
#                 'Sn': {
#                     'n_sites': 6,
#                     'inequiv_groups': [
#                         {
#                             'n_sites': 4,
#                             'inequiv_id': 0,
#                             'dists': [1, 1, 2, 2]
#                         },
#                         {
#                             'n_sites': 2,
#                             'inequiv_id': 1,
#                             'dists': [3, 3]
#                         }
#                     ]
#                 }
#             }
#
#         Merged sites do not contain an ``inequiv_groups`` key and are instead
#         formatted as::
#
#             {
#                 'n_sites': 6
#                 'dists': [1, 1, 1, 2, 2, 2, 2, 3, 3]
#                 )
#             }
#
#         Note that there are now more distances than there are number of sites.
#         This is because n_sites gives the number of bonds to a specific site,
#         whereas the distances are for the complete set of distances for all
#         similar (merged) sites. Similarly, merged next nearest neighbor
#         data can contain more angles than number of sites, however, the
#         general format of the ``nnn_data`` dict is unaltered.
#     """
#     sites = copy.deepcopy(sites)
#     new_sites = []
#
#     for site in sites:
#
#         matched = False
#         for new_site in new_sites:
#             elem_match = site['element'] == new_site['element']
#             geom_match = geometries_match(
#                 site['geometry'], new_site['geometry'], likeness_tol=1)
#             nn_match = nn_summaries_match(
#                 site['nn_data'], new_site['nn_data'],
#                 match_bond_dists=False)
#             nnn_match = nnn_summaries_match(
#                 site['nnn_data'], new_site['nnn_data'],
#                 match_bond_angles=False)
#
#             if elem_match and geom_match and nn_match and nnn_match:
#                 new_site['nn_data'] = _merge_nn_data(site['nn_data'],
#                                                      new_site['nn_data'])
#                 new_site['nnn_data'] = _merge_nnn_data(site['nnn_data'],
#                                                        new_site['nnn_data'])
#                 matched = True
#                 break
#
#         if not matched:
#             # no matches therefore store original site id
#             new_sites.append(site)
#
#     return new_sites
#
#
# def _merge_nn_data(site_nn_data, new_site_nn_data):
#     """Utility function to merge nearest neighbor data.
#
#     See the ``merge_similar_sites`` docstring for information on the format of
#     the merged data.
#
#     Note an error will be thrown if this function is called on two sites that
#     not have matching nearest neighbor summaries (ignoring bond distances).
#     """
#
#     for el in site_nn_data:
#         site_dists = [dist for group in
#                       site_nn_data[el]['inequiv_groups']
#                       for dist in group['dists']]
#
#         if 'inequiv_groups' in new_site_nn_data[el]:
#             # remove inequiv_groups key and group all distances
#             # together
#             groups = new_site_nn_data[el].pop('inequiv_groups')
#             dists = [dist for dist_set in groups
#                      for dist in dist_set['dists']]
#             new_site_nn_data[el]['dists'] = dists + site_dists
#         else:
#             new_site_nn_data[el]['dists'] += site_dists
#
#     return new_site_nn_data
#
#
# def _merge_nnn_data(site_nnn_data, new_site_nnn_data):
#     """Utility function to merge next nearest neighbor data.
#
#     See the ``merge_similar_sites`` docstring for information on the format of
#     the merged data.
#
#     Note an error will be thrown if this function is called on two sites that
#     not have matching next nearest neighbor summaries (ignoring bond angles).
#     """
#     for el in site_nnn_data:
#         for geometry in site_nnn_data[el]:
#             for connectivity in site_nnn_data[el][geometry]:
#                 new_site_nnn_data[el][geometry][connectivity]['angles'].extend(
#                     site_nnn_data[el][geometry][connectivity]['angles'])
#
#     return new_site_nnn_data
