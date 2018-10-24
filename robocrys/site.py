"""
This module provides functions for extracting information about site geometry.

TODO: make functions for
 - distortion of geometry e.g. elongated along an axis
 - connectivity e.g. edge-sharing
 - maybe have custom descriptions for octahedrons, tetrahedron etc
 - handle the case where no geometry type is given, just the CN
"""

from typing import Optional, Iterable, Dict
from collections import defaultdict

from pymatgen.analysis.local_env import CrystalNN, NearNeighbors
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import Structure

from robocrys.fingerprint import get_site_fingerprints


class SiteAnalyzer(object):
    """Class to extract information on site geometry and bonding.

    Args:
        structure: A structure.
        bonded_structure: A structure with nearest neighbour data included.
            One should use a subclass of
            :class:`pymatgen.analysis.local_env.NearNeighbors` such as
            :class:`pymatgen.analysis.local_env.CrystalNN` or
            :class:`pymatgen.analysis.local_env.VoronoiNN`. If no bonded
            structure is provided, the bonding will be calculated using
            :class:`pymatgen.analysis.local_env.CrystalNN`.
        symprec: The tolerance used when determining the symmetry of
            the structure. The symmetry is used to determine if multiple
            sites are symmetrically equivalent.
    """

    def __init__(self,
                 structure: Structure,
                 bonded_structure: Optional[NearNeighbors]=None,
                 symprec: float=0.01):
        self.structure = structure
        self.site_fingerprints = get_site_fingerprints(structure)

        if not bonded_structure:
            cnn = CrystalNN()
            bonded_structure = cnn.get_bonded_structure(structure)

        sga = SpacegroupAnalyzer(structure, symprec=symprec)
        equivalent_sites = sga.get_symmetry_dataset()['equivalent_atoms']

        self.nearest_neighbour_data = []
        for site_index in range(structure.num_sites):
            con_sites = bonded_structure.get_connected_sites(site_index)
            data = tuple({'element': x.site.specie.name,
                          'sym_id': equivalent_sites[x.index],
                          'dist': x.dist} for x in con_sites)
            self.nearest_neighbour_data.append(data)

        self.nearest_neighbour_summary = {}

    def get_site_geometry(self, site_index: int) -> dict:
        """Gets the bonding geometry of a site.

        For example, "octahedral" or "square-planar".

        Args:
            site_index: The site index (zero based).

        Returns:
            The site geometry information as a :obj:`dict` with keys "geometry"
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

        return {'geometry': geometry, 'likeness': parameter[1]}

    def get_nearest_neighbour_info(self, site_index: int) -> Iterable[Dict]:
        """Gets information about the bonded nearest neighbours.

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
        return self.nearest_neighbour_data[site_index]

    def get_nearest_neighbour_summary(self, site_index: int) -> dict:
        """Gets a summary of all the nearest neighbours to a site.

        Args:
            site_index: The site index (zero based).

        Returns:
            A summary of the nearest neighbour information as a dict. Formatted
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
        if site_index in self.nearest_neighbour_summary:
            return self.nearest_neighbour_summary[site_index]

        nn_info = self.get_nearest_neighbour_info(site_index)

        # first group nearest neighbours by element and sym_id
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
