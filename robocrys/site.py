"""
This module provides functions for extracting information about site geometry.

TODO: make functions for
 - distortion of geometry e.g. elongated along an axis
 - connectivity e.g. edge-sharing
 - maybe have custom descriptions for octahedrons, tetrahedron etc
 - handle the case where no geometry type is given, just the CN
"""

from typing import Tuple, Dict, Any, Text, List
from collections import defaultdict

from pymatgen.analysis.graphs import StructureGraph
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from robocrys.fingerprint import get_site_fingerprints


class SiteAnalyzer(object):
    """Class to extract information on site geometry and bonding.

    Args:
        bonded_structure: A bonded structure with nearest neighbor data
            included. For example generated using
            :class:`pymatgen.analysis.local_env.CrystalNN` or
            :class:`pymatgen.analysis.local_env.VoronoiNN`.
        symprec: The tolerance used when determining the symmetry of
            the structure. The symmetry is used to determine if multiple
            sites are symmetrically equivalent.
    """

    def __init__(self, bonded_structure: StructureGraph, symprec: float=0.01):
        self.bonded_structure = bonded_structure
        self.site_fingerprints = get_site_fingerprints(
            bonded_structure.structure)

        sga = SpacegroupAnalyzer(bonded_structure.structure, symprec=symprec)
        self.equivalent_sites = sga.get_symmetry_dataset()['equivalent_atoms']

        self.nearest_neighbor_summary = {}
        self.next_nearest_neighbor_summary = {}

    def get_site_geometry(self, site_index: int) -> Dict[Text, Any]:
        """Gets the bonding geometry of a site.

        For example, "octahedral" or "square-planar".

        Args:
            site_index: The site index (zero based).

        Returns:
            The site geometry information as a :obj:`dict` with keys "type"
            and "likeness", corresponding to the geometry type (e.g. octahedral)
            and order parameter, respectively.
        """
        # get fingerprint as a list of tuples, e.g. [("op name", val), ...]
        site_fingerprint = list(self.site_fingerprints[site_index].items())

        # get coordination number with largest weight, ignore op names with
        # just the coordination number weight (e.g. containing "wt")
        parameter = max(site_fingerprint,
                        key=lambda x: x[1] if "wt" not in x[0] else 0)

        # return the geometry type without the CN at the end, e.g.
        # "square co-planar CN_4" -> "square co-planar"
        geometry = " ".join(parameter[0].split()[:-1])

        return {'type': geometry, 'likeness': parameter[1]}

    def get_nearest_neighbor_summary(self, site_index: int) -> Dict[Text, Any]:
        """Gets a summary of all the nearest neighbors to a site.

        Args:
            site_index: The site index (zero based).

        Returns:
            A summary of the nearest neighbor information as a dict. Formatted
            as::

                {
                    'Sn': {
                        'n_sites': 6
                        'sym_groups': (
                            {
                                'n_sites': 4
                                'sym_id': 0
                                'dists': [1, 1, 1, 2, 2, 2]
                            },
                            {
                                'n_sites': 2
                                'sym_id': 1
                                'dists': [3, 3]
                            }
                        )
                    }
                }
        """
        if site_index in self.nearest_neighbor_summary:
            return self.nearest_neighbor_summary[site_index]

        nn_info = self._get_nearest_neighbor_info(site_index)

        # first group nearest neighbors by element and sym_id
        # e.g. grouped_nn looks like {'el': {'sym_id': [sites]}}
        grouped_nn = defaultdict(lambda: defaultdict(list))
        for site in nn_info:
            grouped_nn[site['element']][site['sym_id']].append(site)

        data = {}
        for element, sym_data in grouped_nn.items():
            n_sites = sum([len(sites) for sites in sym_data.values()])
            sym_groups = tuple({
                'n_sites': len(sites),
                'sym_id': sym_id,
                'dists': [x['dist'] for x in sites]
            } for sym_id, sites in sym_data.items())
            data[element] = {'n_sites': n_sites, 'sym_groups': sym_groups}

        return data

    def get_next_nearest_neighbor_summary(self,
                                          site_index: int) -> Dict[Text, Any]:
        """Gets a summary of the next nearest neighbor connectivity.

        Args:
            site_index: The site index (zero based).

        Returns:
            A summary of the next nearest neighbor information as a dict.
            Formatted as::

                {
                    'Sn': {
                        'corner-sharing': {
                            'n_sites': 8,
                            'geometries': ('octahedral')
                        }
                    }
                }
        """
        if site_index in self.next_nearest_neighbor_summary:
            return self.next_nearest_neighbor_summary[site_index]

        nnn_info = self._get_next_nearest_neighbor_info(site_index)

        # first group next nearest neighbors by element and connectivity
        # e.g. grouped_nnn looks like {'el': {'connectivity': [sites]}}
        grouped_nnn = defaultdict(lambda: defaultdict(list))
        for site in nnn_info:
            grouped_nnn[site['element']][site['connectivity']].append(site)

        data = {}
        for element, con_data in grouped_nnn.items():
            el_data = {'n_sites': sum([len(sites) for sites in
                                       con_data.values()])}
            con_groups = {
                connectivity: {
                    'n_sites': len(sites),
                    'geometries': tuple(set([site['geometry']['type']
                                             for site in sites]))
                } for connectivity, sites in con_data.items()}
            el_data.update(con_groups)
            data[element] = el_data

        return data

    def _get_nearest_neighbor_info(
            self, site_index: int) -> Tuple[Dict[Text, Any]]:
        """Gets information about the bonded nearest neighbors.

        Args:

            site_index: The site index (zero based).

        Returns:
            For each site bonded to ``site_index``, returns a :obj:`dict`
            with the format::

                {'element': el, 'sym_id': i, 'dist': distance}

            The ``sym_id`` property is the symmetry index for the site. E.g. if
            two sites are symmetrically  equivalent then they will have the
            same ``sym_id``.
        """

        nn_sites = self.bonded_structure.get_connected_sites(site_index)
        return tuple({'element': site.site.specie.name,
                      'sym_id': self.equivalent_sites[site.index],
                      'dist': site.dist} for site in nn_sites)

    def _get_next_nearest_neighbor_info(
            self, site_index: int) -> List[Dict[Text, Any]]:
        """Gets information about the bonded next nearest neighbors.

        Args:
            site_index: The site index (zero based).

        Returns:
            A list of the next nearest neighbor information. For each next
            nearest neighbor site, returns a :obj:`dict` with the format::

                {'element': el, 'connectivity': con, 'geometry': geom}

            The ``connectivity`` property is the connectivity type to the
            next nearest neighbor, e.g. "face-sharing", "corner-sharing", or
            "edge-sharing". The ``geometry`` property gives the geometry of the
            next nearest neighbor site. See the ``get_site_geometry`` method for
            the format of this data.
        """

        nn_sites = self.bonded_structure.get_connected_sites(site_index)
        next_nn_sites = [site for nn_site in nn_sites for site in
                         self.bonded_structure.get_connected_sites(
                             nn_site.index, jimage=nn_site.jimage)]

        nn_sites_set = set((site.index, site.jimage) for site in nn_sites)

        next_nn_summary = []
        for nnn_site in next_nn_sites:
            if nnn_site.index == site_index and nnn_site.jimage == (0, 0, 0):
                # skip the nnn site if it is the original atom of interest
                continue

            sites = set((site.index, site.jimage) for site in
                        self.bonded_structure.get_connected_sites(
                            nnn_site.index, jimage=nnn_site.jimage))
            n_shared_atoms = len(nn_sites_set.intersection(sites))

            if n_shared_atoms == 1:
                connectivity = 'corner-sharing'
            elif n_shared_atoms == 2:
                connectivity = 'edge-sharing'
            else:
                connectivity = 'face-sharing'

            geometry = self.get_site_geometry(nnn_site.index)

            next_nn_summary.append({
                'element': nnn_site.site.specie.name,
                'connectivity': connectivity,
                'geometry': geometry})

        return next_nn_summary
