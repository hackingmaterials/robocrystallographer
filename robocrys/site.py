"""
This module provides functions for extracting information about site geometry.

TODO: distortion of geometry e.g. elongated along an axis
TODO: maybe have custom descriptions for octahedron, tetrahedron etc
TODO: handle the case where no geometry type is given, just the CN
"""
import copy
from collections import defaultdict
from typing import Dict, Any, List

import numpy as np

from pymatgen import Composition
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.core.periodic_table import get_el_sp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.coord import get_angle
from robocrys.fingerprint import get_site_fingerprints
from robocrys.util import connected_geometries, get_el


class SiteAnalyzer(object):
    """Class to extract information on site geometry and bonding.

    Args:
        bonded_structure: A bonded structure with nearest neighbor data
            included. For example generated using
            :class:`pymatgen.analysis.local_env.CrystalNN` or
            :class:`pymatgen.analysis.local_env.VoronoiNN`.
        use_symmetry: Whether to use symmetry to determine if sites are
            inequivalent. If ``False``, the site geometry and (next) nearest
            neighbor information will be used.
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
    """

    def __init__(self,
                 bonded_structure: StructureGraph,
                 use_symmetry: bool = False,
                 symprec: float = 0.01,
                 use_iupac_formula: bool = True):
        self.bonded_structure = bonded_structure
        self.site_fingerprints = get_site_fingerprints(
            bonded_structure.structure)

        if use_symmetry:
            sga = SpacegroupAnalyzer(bonded_structure.structure,
                                     symprec=symprec)
            sym_data = sga.get_symmetry_dataset()
            self.equivalent_sites = sym_data['equivalent_atoms']

        else:
            self.equivalent_sites = self._calculate_equivalent_sites()

        self.use_iupac_formula = use_iupac_formula

    def get_site_summary(self, site_index: int) -> Dict[str, Any]:
        element = str(self.bonded_structure.structure[site_index].specie)
        geometry = self.get_site_geometry(site_index)
        nn_data = self.get_nearest_neighbor_data(site_index)
        nnn_data = self.get_next_nearest_neighbor_data(site_index)

        site = {'element': element, 'geometry': geometry, 'nn_data': nn_data,
                'nnn_data': nnn_data}

        if site_is_connected_polyhedra(site):
            nn_els = "".join(["{}{}".format(get_el(el), el_data['n_sites'])
                              for el, el_data in nn_data.items()])
            comp = Composition(get_el(element) + nn_els)

            site['polyhedra_formula'] = comp.get_reduced_formula_and_factor(
                iupac_ordering=self.use_iupac_formula)[0]

        return site

    def get_site_geometry(self, site_index: int) -> Dict[str, Any]:
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

    def get_nearest_neighbor_data(self, site_index: int,
                                  split_into_groups: bool = True
                                  ) -> Dict[str, Any]:
        """Gets a summary of all the nearest neighbors to a site.

        Args:
            site_index: The site index (zero based).
            split_into_groups: Whether to split the nearest neighbors into
                groups of inequivalent atoms.

        Returns:
            A summary of the nearest neighbor information as a dict. If
            ``split_into_groups=True``, the data is formatted as::

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

            If ``split_into_groups=False``, the data is formatted as::

                {
                    'Sn': {
                        'n_sites': 6,
                        'dists': [1, 1, 2, 2, 3, 3]
                    }
                }

        """
        nn_info = self._get_nearest_neighbor_info(
            site_index, inc_inequiv_id=split_into_groups)

        if split_into_groups:
            # first group nearest neighbors by element and inequiv_id
            # e.g. grouped_nn looks like {'el': {'sym_id': [sites]}}
            grouped_nn = defaultdict(lambda: defaultdict(list))
            for site in nn_info:
                grouped_nn[site['element']][site['inequiv_id']].append(site)

            data = {}
            for element, sym_data in grouped_nn.items():
                n_sites = sum([len(sites) for sites in sym_data.values()])
                sym_groups = [
                    {'n_sites': len(sites),
                     'inequiv_id': sym_id,
                     'dists': [x['dist'] for x in sites]
                     } for sym_id, sites in sym_data.items()]
                data[element] = {'n_sites': n_sites,
                                 'inequiv_groups': sym_groups}

        else:
            # first group nearest neighbors by element
            # e.g. grouped_nn looks like {'el': [sites]}
            grouped_nn = defaultdict(list)
            for site in nn_info:
                grouped_nn[site['element']].append(site)

            data = {element: {'n_sites': len(sites),
                              'dists': [x['dist'] for x in sites]}
                    for element, sites in grouped_nn.items()}

        return data

    def get_next_nearest_neighbor_data(self, site_index: int
                                       ) -> Dict[str, Any]:
        """Gets a summary of the next nearest neighbor connectivity.

        Args:
            site_index: The site index (zero based).

        Returns:
            A summary of the next nearest neighbor information as a dict.
            Formatted as::

                {
                    'Sn': {
                        'octahedral': {
                            'corner-sharing': {
                                'n_sites': 8,
                                'angles': [180, 180, 180, ...]
                            }
                        }
                    }
                }

        """
        nnn_info = self._get_next_nearest_neighbor_info(site_index)

        # group next nearest neighbors by element, connectivity and geometry
        # e.g. grouped_nnn looks like {el: {connectivity: {geometry: [sites]}}}
        grouped_nnn = defaultdict(
            lambda: defaultdict(lambda: defaultdict(list)))

        for site in nnn_info:
            grouped_nnn[site['element']][
                site['geometry']['type']][site['connectivity']].append(site)

        nnn_data = {}
        for element, geom_data in grouped_nnn.items():
            nnn_el_data = {}
            for geometry, con_data in geom_data.items():
                nnn_geom_data = {}
                for connectivity, sites in con_data.items():
                    nnn_geom_data[connectivity] = {
                        'n_sites': len(sites),
                        'angles': [angle for site in sites
                                   for angle in site['angles']]}
                nnn_el_data[geometry] = nnn_geom_data
            nnn_data[element] = nnn_el_data
        return nnn_data

    def get_inequivalent_site_ids(self, site_ids: List[int]) -> List[int]:
        """Gets the inequivalent sites from a list of site indices.

        Args:
            site_ids: The site indices.

        Returns:
            The inequivalent site indices. For example, if a structure has 4
            sites where the first two are equivalent and the last two are
            inequivalent. If  ``site_ids=[0, 1, 2, 3]`` the output will be::

                [0, 2, 3]

        """
        return list(set(self.equivalent_sites[i] for i in site_ids))

    def _get_nearest_neighbor_info(self, site_index: int,
                                   inc_inequiv_id: bool = True
                                   ) -> List[Dict[str, Any]]:
        """Gets information about the bonded nearest neighbors.

        Args:

            site_index: The site index (zero based).
            inc_inequiv_id: Whether to include the inequivalent site ids.

        Returns:
            For each site bonded to ``site_index``, returns a :obj:`dict`
            with the format::

                {'element': el, 'dist': distance}

            If ``inc_inequiv_id=True``, the data will have an additional key
            ``'inequiv_id'`` corresponding to the inequivalent site index. E.g.
            if two sites are structurally/symmetrically equivalent (depending
            on the value of ``self.use_symmetry`` then they will have the same
            ``inequiv_id``.
        """

        nn_sites = self.bonded_structure.get_connected_sites(site_index)

        if inc_inequiv_id:
            return [{'element': str(site.site.specie),
                     'inequiv_id': self.equivalent_sites[site.index],
                     'dist': site.dist} for site in nn_sites]
        else:
            return [{'element': str(site.site.specie),
                     'dist': site.dist} for site in nn_sites]

    def _get_next_nearest_neighbor_info(
            self, site_index: int) -> List[Dict[str, Any]]:
        """Gets information about the bonded next nearest neighbors.

        Args:
            site_index: The site index (zero based).

        Returns:
            A list of the next nearest neighbor information. For each next
            nearest neighbor site, returns a :obj:`dict` with the format::

                {'element': el, 'connectivity': con, 'geometry': geom,
                 'angles': angles}

            The ``connectivity`` property is the connectivity type to the
            next nearest neighbor, e.g. "face-sharing", "corner-sharing", or
            "edge-sharing". The ``geometry`` property gives the geometry of the
            next nearest neighbor site. See the ``get_site_geometry`` method for
            the format of this data. The ``angles`` property gives the bond
            angles between the site and the next nearest neighbour. Returned as
            a :obj:`list` of :obj:`int`. Multiple bond angles are given when
            the two sites share more than nearest neighbor (e.g. if they are
            face-sharing or edge-sharing).
        """

        def get_coords(a_site_index, a_site_image):
            return np.asarray(
                self.bonded_structure.structure.lattice.get_cartesian_coords(
                    self.bonded_structure.structure.frac_coords[a_site_index] +
                    a_site_image))

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
            shared_sites = nn_sites_set.intersection(sites)
            n_shared_atoms = len(shared_sites)

            if n_shared_atoms == 1:
                connectivity = 'corner-sharing'
            elif n_shared_atoms == 2:
                connectivity = 'edge-sharing'
            else:
                connectivity = 'face-sharing'

            site_coords = get_coords(site_index, (0, 0, 0))
            nnn_site_coords = get_coords(nnn_site.index, nnn_site.jimage)
            nn_site_coords = [get_coords(nn_site_index, nn_site_image)
                              for nn_site_index, nn_site_image in shared_sites]

            # can't just use Structure.get_angles to calculate angles as it
            # doesn't take into account the site image
            angles = [get_angle(site_coords - x, nnn_site_coords - x)
                      for x in nn_site_coords]

            geometry = self.get_site_geometry(nnn_site.index)

            next_nn_summary.append({
                'element': str(nnn_site.site.specie),
                'connectivity': connectivity,
                'geometry': geometry,
                'angles': angles})

        return next_nn_summary

    def _calculate_equivalent_sites(self,
                                    geometry_likeness_tol: float = 0.1,
                                    bond_dist_tol: float = 0.1,
                                    bond_angle_tol: float = 0.1
                                    ) -> List[int]:
        """Determines the indices of the structurally inequivalent sites.

        Args:
            geometry_likeness_tol: The tolerance to use when determining whether
                two geometry likeness parameters are the same.
            bond_dist_tol: The tolerance used to determine if two bond lengths
                are the same.
            bond_angle_tol: The tolerance used to determine if two bond angles
                are the same.

        Two sites are considered equivalent if they are the same element,
        geometry, and (next) nearest neighbor summaries.

        Returns:
            A :obj:`list` of indices mapping each site in the structure to a
            structurally equivalent site. For example, if the first two sites
            are equivalent and the last two are both inequivalent, the data will
            be formatted as::

                [0, 0, 2, 3]

        """
        # TODO: Consider using site fingerprint rather than geometry type.
        inequiv_sites = {}
        equivalent_sites = []

        for site_id, site in enumerate(self.bonded_structure.structure):
            element = get_el_sp(site.specie)
            geometry = self.get_site_geometry(site_id)
            nn_data = self.get_nearest_neighbor_data(
                site_id, split_into_groups=False)
            nnn_data = self.get_next_nearest_neighbor_data(site_id)

            matched = False
            for inequiv_id, inequiv_site in inequiv_sites.items():
                elem_match = element == inequiv_site['element']
                geom_match = geometries_match(
                    geometry, inequiv_site['geometry'],
                    likeness_tol=geometry_likeness_tol)
                nn_match = nn_summaries_match(
                    nn_data, inequiv_site['nn_data'],
                    bond_dist_tol=bond_dist_tol)
                nnn_match = nnn_summaries_match(
                    nnn_data, inequiv_site['nnn_data'],
                    bond_angle_tol=bond_angle_tol)

                if elem_match and geom_match and nn_match and nnn_match:
                    equivalent_sites.append(inequiv_id)
                    matched = True
                    break

            if not matched:
                # no matches therefore store original site id
                equivalent_sites.append(site_id)
                site_data = {'element': element,
                             'geometry': geometry,
                             'nn_data': nn_data,
                             'nnn_data': nnn_data}
                inequiv_sites[site_id] = site_data

        return equivalent_sites


def geometries_match(geometry_a: Dict[str, Any], geometry_b: Dict[str, Any],
                     likeness_tol: float = 0.01) -> bool:
    """Determine whether two site geometries match.

    Geometry data should be formatted the same as produced by
    :meth:`robocrys.site.SiteAnalyzer.get_site_geometry`.

    Args:
        geometry_a: The first set of geometry data.
        geometry_b: The second set of geometry data.
        likeness_tol: The tolerance to use when determining whether two
            geometry likeness parameters are the same.

    Returns:
        Whether the two geometries are the same.
    """
    # TODO: Handle case when geometry types are both None
    return (geometry_a['type'] == geometry_b['type'] and
            abs(geometry_a['likeness'] - geometry_b['likeness']) < likeness_tol)


def nn_summaries_match(nn_summary_a: Dict[str, Any],
                       nn_summary_b: Dict[str, Any],
                       bond_dist_tol: float = 0.1,
                       match_bond_dists: bool = True) -> bool:
    """Determine whether two nearest neighbors summaries match.

    Nearest neighbor data should be formatted the same as produced by
    :meth:`robocrys.site.SiteAnalyzer.get_nearest_neighbor_data` with
    ``split_into_groups=False``.

    Args:
        nn_summary_a: The first set of nearest neighbor data.
        nn_summary_b: The second set of nearest neighbor data.
        bond_dist_tol: The tolerance used to determine if two bond lengths
            are the same.
        match_bond_dists: Whether to consider bond distances when matching.

    Returns:
        Whether the two nearest neighbor summaries match.
    """

    nn_summary_a = copy.deepcopy(nn_summary_a)
    nn_summary_b = copy.deepcopy(nn_summary_b)

    for elem_a, data_a in nn_summary_a.items():
        if elem_a not in nn_summary_b:
            return False

        data_b = nn_summary_b[elem_a]
        if data_a['n_sites'] != data_b['n_sites']:
            return False

        if data_a['n_sites'] > 0 and match_bond_dists:
            diff_dists = (np.array(sorted(data_a['dists'])) -
                          np.array(sorted(data_b['dists'])))
            if not all(np.abs(diff_dists) < bond_dist_tol):
                return False

        nn_summary_b.pop(elem_a)

    if nn_summary_b:
        # if summary_b still contains data then summary_b has additional
        # neighbors not in summary_a
        return False
    else:
        return True


def nnn_summaries_match(nnn_summary_a: Dict[str, Any],
                        nnn_summary_b: Dict[str, Any],
                        bond_angle_tol: float = 0.1,
                        match_bond_angles: bool = True):
    """Determine whether two next nearest neighbors summaries match.

    Next nearest neighbor data should be formatted the same as produced by
    :meth:`robocrys.site.SiteAnalyzer.get_next_nearest_neighbor_data`.

    Args:
        nnn_summary_a: The first set of next nearest neighbor data.
        nnn_summary_b: The second set of next nearest neighbor data.
        bond_angle_tol: The tolerance used to determine if two bond angles
            are the same.
        match_bond_angles: Whether to consider bond angles when matching.

    Returns:
        Whether the two next nearest neighbor summaries match.
    """
    nnn_summary_a = copy.deepcopy(nnn_summary_a)
    nnn_summary_b = copy.deepcopy(nnn_summary_b)

    for elem_a, geom_data_a in nnn_summary_a.items():
        if elem_a not in nnn_summary_b:
            return False

        geom_data_b = nnn_summary_b[elem_a]
        for geometry_a, con_data_a in geom_data_a.items():
            if geometry_a not in geom_data_b:
                return False

            con_data_b = geom_data_b[geometry_a]
            for con_type_a, data_a in con_data_a.items():
                if con_type_a not in con_data_b:
                    return False

                data_b = con_data_b[con_type_a]

                if data_a['n_sites'] != data_b['n_sites']:
                    return False

                if match_bond_angles:
                    diff_angles = (np.array(sorted(data_a['angles'])) -
                                   np.array(sorted(data_b['angles'])))
                    if not all(np.abs(diff_angles) < bond_angle_tol):
                        return False

                con_data_b.pop(con_type_a)

            if con_data_b:
                # if geom_data_b still contains data then it has additional
                # data relative to geom_data_a
                return False

            geom_data_b.pop(geometry_a)

        if geom_data_b:
            # if con_data_b still contains data then it has additional
            # data relative to con_data_a
            return False

        nnn_summary_b.pop(elem_a)

    if nnn_summary_b:
        # if summary_b still contains data then summary_b has additional
        # next nearest neighbors not in summary_a
        return False
    else:
        return True


def site_is_connected_polyhedra(site: Dict[str, Any]) -> bool:
    """Determines whether a site is corner, edge, or face sharing polyhedra.

    Essentially, if two octahedra are bonded, two tetrahedra are bonded or
    a mixture of tetrahedra and octahedra are bonded, then you have some sort
    of edge, corner or face-sharing connectivity. This description does not
    indicate how persistent the connectivity is throughout the whole cell.

    Args:
        site: A site, formatted as a :ob:`dict` with the keys:

            - ``'element'`` (``str``): The element of the site.
            - ``'geometry'`` (``dict``): The geometry, as output by
                :meth:`SiteAnalyzer.get_site_geometry`.
            - ``'nn_data'`` (``dict``): The nearest neighbor data, as output by
                :meth:`SiteAnalyzer.get_nearest_neighbor_data`.
            - ``'nnn_data'`` (``dict``): The next nearest neighbor data, as
                given by :meth:`SiteAnalyzer.get_next_nearest_neighbor_data`.

    Returns:
        Whether the site is a edge, corner, or face sharing octahedral or
        tetrahedral site.
    """
    if (site['geometry']['type'] in connected_geometries and
            any([nnn_site_geometry in connected_geometries
                 for nnn_geometry_data in site['nnn_data'].values()
                 for nnn_site_geometry in nnn_geometry_data.keys()])):
        return True
    else:
        return False


def merge_similar_sites(sites: List[Dict[str, Any]]):
    """Merges sites with the same properties except bond angles and distances.

    Args:
        sites: A list of sites. Each site is formatted as a :ob:`dict` with the
            keys:

            - ``'element'`` (``str``): The element of the site.
            - ``'geometry'`` (``dict``): The geometry, as output by
                :meth:`SiteAnalyzer.get_site_geometry`.
            - ``'nn_data'`` (``dict``): The nearest neighbor data, as output by
                :meth:`SiteAnalyzer.get_nearest_neighbor_data`.
            - ``'nnn_data'`` (``dict``): The next nearest neighbor data, as
                given by :meth:`SiteAnalyzer.get_next_nearest_neighbor_data`.

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
