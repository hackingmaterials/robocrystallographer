"""This module provides a class to extract geometry and neighbor information.

Todo:
    * distortion of geometry e.g. elongated along an axis
"""

from __future__ import annotations

from collections import defaultdict
from typing import TYPE_CHECKING

import numpy as np
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import get_el_sp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.coord import get_angle
from pymatgen.util.string import formula_double_format

from robocrys.condense.fingerprint import get_site_fingerprints
from robocrys.util import connected_geometries, defaultdict_to_dict, get_el

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any


class SiteAnalyzer:
    """Class to extract information on site geometry and bonding.

    Attributes:
        symmetry_labels: A :obj:`dict` mapping the site indices to the symmetry
            label for that site. If two sites are symmetrically equivalent they
            share the same symmetry label. The numbering begins at 1 for each
            element in the structure.
        equivalent_sites: A :obj:`list` of indices mapping each site in
            the structure to a symmetrically or structurally equivalent site,
            depending on the value of ``use_symmetry_equivalent_sites``.

    Args:
        bonded_structure: A bonded structure with nearest neighbor data
            included. For example generated using
            :class:`pymatgen.analysis.local_env.CrystalNN` or
            :class:`pymatgen.analysis.local_env.VoronoiNN`.
        use_symmetry_equivalent_sites: Whether to use symmetry to determine if
            sites are inequivalent. If ``False``, the site geometry and (next)
            nearest neighbor information will be used.
        symprec: The tolerance used when determining the symmetry of
            the structure. The symmetry can used both to determine if multiple
            sites are symmetrically equivalent and to obtain the symmetry labels
            for each site.
        minimum_geometry_op: The minimum geometrical order parameter for a
            geometry match to be returned.
        use_iupac_formula (bool, optional): Whether to order formulas
            by the iupac "electronegativity" series, defined in
            Table VI of "Nomenclature of Inorganic Chemistry (IUPAC
            Recommendations 2005)". This ordering effectively follows
            the groups and rows of the periodic table, except the
            Lanthanides, Actanides and hydrogen. If set to ``False``, the
            elements will be ordered according to the electronegativity values.
    """

    def __init__(
        self,
        bonded_structure: StructureGraph,
        use_symmetry_equivalent_sites: bool = False,
        symprec: float = 0.01,
        minimum_geometry_op: float = 0.4,
        use_iupac_formula: bool = True,
    ):
        self.bonded_structure = bonded_structure
        self.use_iupac_formula = use_iupac_formula
        self.minimum_geometry_op = minimum_geometry_op

        self.site_fingerprints = get_site_fingerprints(bonded_structure.structure)

        sga = SpacegroupAnalyzer(bonded_structure.structure, symprec=symprec)
        equivalent_sites = sga.get_symmetry_dataset().equivalent_atoms

        if use_symmetry_equivalent_sites:
            self.equivalent_sites = list(equivalent_sites)
        else:
            self.equivalent_sites = self._calculate_equivalent_sites()

        self.symmetry_labels = self._calculate_symmetry_labels(
            equivalent_sites.tolist()
        )

    def get_site_geometry(self, site_index: int) -> dict[str, str | float]:
        """Gets the bonding geometry of a site.

        For example, "octahedral" or "square-planar".

        Args:
            site_index: The site index (zero based).

        Returns:
            The site geometry information formatted at as::

                {'type': geometry_type, 'likeness': order_parameter}

            Where ``geometry_type`` is a :obj:`str` corresponding to the
            geometry type (e.g. octahedral) and ``order_parameter`` is a
            :obj:`float` indicating whether how close the geometry is to the
            perfect geometry. If the largest geometrical order parameter falls
            beneath :attr:`robocrys.site.SiteAnalyzer.minimum_geometry_op`, the
            geometry type will be returned as "X-coordinate", where X is the
            coordination number.
        """
        # get fingerprint as a list of tuples, e.g. [("op name", val), ...]
        site_fingerprint: list[tuple[str, int]] = list(
            self.site_fingerprints[site_index].items()
        )

        # get coordination number with largest weight, ignore op names with
        # just the coordination number weight (e.g. containing "wt")
        parameter = max(site_fingerprint, key=lambda x: x[1] if "wt" not in x[0] else 0)

        if parameter[1] < self.minimum_geometry_op:
            # the largest filtered weight is less than the tolerance for determining the
            # coordination geometry (i.e., we can no longer say the site is octahedral
            # or tetrahedral etc). Now we don't care about the actual geometry, we just
            # want the most likely coordination number; take this from the largest
            # of all weights
            parameter = max(site_fingerprint, key=lambda x: x[1])
            cn = parameter[0].split()[-1].split("_")[-1]
            geometry = f"{cn}-coordinate"
            likeness = 1.0

        else:
            # return the geometry type without the CN at the end, e.g.
            # "square co-planar CN_4" -> "square co-planar"
            geometry = " ".join(parameter[0].split()[:-1])
            geometry = "single-bond" if geometry == "sgl_bd" else geometry
            likeness = parameter[1]

        return {"type": geometry, "likeness": likeness}

    def get_nearest_neighbors(
        self, site_index: int, inc_inequivalent_site_index: bool = True
    ) -> list[dict[str, Any]]:
        """Gets information about the bonded nearest neighbors.

        Args:
            site_index: The site index (zero based).
            inc_inequivalent_site_index: Whether to include the inequivalent
                site indices in the nearest neighbor information.

        Returns:
            For each site bonded to ``site_index``, returns a :obj:`dict`
            with the format::

                {'element': el, 'dist': distance}

            If ``inc_inequivalent_site_index=True``, the data will have an
            additional key ``'inequiv_index'`` corresponding to the inequivalent
            site index. E.g. if two sites are structurally/symmetrically
            equivalent (depending on the value of ``self.use_symmetry_equivalent_sites``
            then they will have the same ``inequiv_index``.
        """
        nn_sites = self.bonded_structure.get_connected_sites(site_index)

        if inc_inequivalent_site_index:
            return [
                {
                    "element": str(site.site.specie),
                    "inequiv_index": self.equivalent_sites[site.index],
                    "dist": site.dist,
                }
                for site in nn_sites
            ]
        return [
            {"element": str(site.site.specie), "dist": site.dist} for site in nn_sites
        ]

    def get_next_nearest_neighbors(
        self, site_index: int, inc_inequivalent_site_index: bool = True
    ) -> list[dict[str, Any]]:
        """Gets information about the bonded next nearest neighbors.

        Args:
            site_index: The site index (zero based).
            inc_inequivalent_site_index: Whether to include the inequivalent
                site indices.

        Returns:
            A list of the next nearest neighbor information. For each next
            nearest neighbor site, returns a :obj:`dict` with the format::

                {'element': el, 'connectivity': con, 'geometry': geom,
                 'angles': angles, 'distance': distance}

            The ``connectivity`` property is the connectivity type to the
            next nearest neighbor, e.g. "face", "corner", or
            "edge". The ``geometry`` property gives the geometry of the
            next nearest neighbor site. See the ``get_site_geometry`` method for
            the format of this data. The ``angles`` property gives the bond
            angles between the site and the next nearest neighbour. Returned as
            a :obj:`list` of :obj:`int`. Multiple bond angles are given when
            the two sites share more than nearest neighbor (e.g. if they are
            face-sharing or edge-sharing). The ``distance`` property gives the
            distance between the site and the next nearest neighbor.
            If ``inc_inequivalent_site_index=True``, the data will have an
            additional key ``'inequiv_index'`` corresponding to the inequivalent
            site index. E.g. if two sites are structurally/symmetrically
            equivalent (depending on the value of ``self.use_symmetry_equivalent_sites``
            then they will have the same ``inequiv_index``.
        """

        def get_coords(a_site_index, a_site_image):
            return np.asarray(
                self.bonded_structure.structure.lattice.get_cartesian_coords(
                    self.bonded_structure.structure.frac_coords[a_site_index]
                    + a_site_image
                )
            )

        nn_sites = self.bonded_structure.get_connected_sites(site_index)
        next_nn_sites = [
            site
            for nn_site in nn_sites
            for site in self.bonded_structure.get_connected_sites(
                nn_site.index, jimage=nn_site.jimage
            )
        ]

        nn_sites_set = {(site.index, site.jimage) for site in nn_sites}

        seen_nnn_sites = set()
        next_nn_summary = []
        for nnn_site in next_nn_sites:
            if (
                nnn_site.index == site_index
                and nnn_site.jimage == (0, 0, 0)
                or (nnn_site.index, nnn_site.jimage) in seen_nnn_sites
            ):
                # skip the nnn site if it is the original atom of interest
                continue

            seen_nnn_sites.add((nnn_site.index, nnn_site.jimage))

            sites = {
                (site.index, site.jimage)
                for site in self.bonded_structure.get_connected_sites(
                    nnn_site.index, jimage=nnn_site.jimage
                )
            }
            shared_sites = nn_sites_set.intersection(sites)
            n_shared_atoms = len(shared_sites)

            if n_shared_atoms == 1:
                connectivity = "corner"
            elif n_shared_atoms == 2:
                connectivity = "edge"
            else:
                connectivity = "face"

            site_coords = get_coords(site_index, (0, 0, 0))
            nnn_site_coords = get_coords(nnn_site.index, nnn_site.jimage)
            nn_site_coords = [
                get_coords(nn_site_index, nn_site_image)
                for nn_site_index, nn_site_image in shared_sites
            ]

            # can't just use Structure.get_angles to calculate angles as it
            # doesn't take into account the site image
            angles = [
                get_angle(site_coords - x, nnn_site_coords - x) for x in nn_site_coords
            ]

            distance = np.linalg.norm(site_coords - nnn_site_coords)

            geometry = self.get_site_geometry(nnn_site.index)

            summary = {
                "element": str(nnn_site.site.specie),
                "connectivity": connectivity,
                "geometry": geometry,
                "angles": angles,
                "distance": distance,
            }

            if inc_inequivalent_site_index:
                summary["inequiv_index"] = self.equivalent_sites[nnn_site.index]

            next_nn_summary.append(summary)

        return next_nn_summary

    def get_site_summary(self, site_index: int) -> dict[str, Any]:
        """Gets a summary of the site information.

        Args:
            site_index: The site index (zero based).

        Returns:
            A summary of the site information, formatted as::

                {
                    'element': 'Mo4+',
                    'geometry': {
                        'likesness': 0.5544,
                        'type': 'pentagonal pyramidal'
                    },
                    'nn': [2, 2, 2, 2, 2, 2],
                    'nnn': {'edge': [0, 0, 0, 0, 0, 0]},
                    'poly_formula': 'S6',
                    'sym_labels': (1,)
                }

            Where ``element`` is the species string (if the species has
            oxidation states, these will be included in the string). The
            ``geometry`` key is the geometry information as produced by
            :meth:`SiteAnalyzer.get_site_geometry`. The `nn` key lists
            the site indices of the nearest neighbor bonding sites. Note the
            inequivalent site index is given for each site. The `nnn` key gives
            the next nearest neighbor information, broken up by the connectivity
            to that neighbor. The ``poly_formula`` key gives the formula of the
            bonded nearest neighbors. ``poly_formula`` will be ``None`` if the
            site geometry is not in :data:`robocrys.util.connected_geometries`.
            The ``sym_labels`` key gives the symmetry labels of the site. If
            two sites are symmetrically equivalent they share the same symmetry
            label. The numbering begins at 1 for each element in the structure.
            If :attr:`SiteAnalyzer.use_symmetry_inequivalnt_sites` is ``False``,
            each site may have more than one symmetry label, as structural
            features have instead been used to determine the site equivalences,
            i.e. two sites are symmetrically distinct but share the same
            geometry, nearest neighbor and next nearest neighbor properties.
        """
        element = str(self.bonded_structure.structure[site_index].specie)
        geometry = self.get_site_geometry(site_index)

        nn_sites = self.get_nearest_neighbors(
            site_index, inc_inequivalent_site_index=True
        )
        nn_indices = [nn_site["inequiv_index"] for nn_site in nn_sites]

        nnn_sites = self.get_next_nearest_neighbors(
            site_index, inc_inequivalent_site_index=True
        )
        nnn = defaultdict(list)
        for nnn_site in nnn_sites:
            nnn[nnn_site["connectivity"]].append(nnn_site["inequiv_index"])

        equiv_sites = [
            i
            for i in range(len(self.equivalent_sites))
            if self.equivalent_sites[i] == self.equivalent_sites[site_index]
        ]
        sym_labels = tuple({self.symmetry_labels[x] for x in equiv_sites})

        poly_formula = self._get_poly_formula(geometry, nn_sites, nnn_sites)

        return {
            "element": element,
            "geometry": geometry,
            "nn": nn_indices,
            "nnn": dict(nnn),
            "poly_formula": poly_formula,
            "sym_labels": sym_labels,
        }

    def get_bond_distance_summary(self, site_index: int) -> dict[int, list[float]]:
        """Gets the bond distance summary for a site.

        Args:
            site_index: The site index (zero based).

        Returns:
            The bonding data for the site, formatted as::

                {to_site: [dist_1, dist_2, dist_3, ...]}

            Where ``to_site`` is the index of a nearest neighbor site
            and ``dist_1`` etc are the bond distances as :obj:`float`.
        """
        bonds = defaultdict(list)
        for nn_site in self.get_nearest_neighbors(site_index):
            to_site = nn_site["inequiv_index"]
            bonds[to_site].append(nn_site["dist"])

        return defaultdict_to_dict(bonds)

    def get_connectivity_angle_summary(
        self, site_index: int
    ) -> dict[int, dict[str, list[float]]]:
        """Gets the connectivity angle summary for a site.

        The connectivity angles are the angles between a site and its
        next nearest neighbors.

        Args:
            site_index: The site index (zero based).

        Returns:
            The connectivity angle data for the site, formatted as::

                {
                    to_site: {
                        connectivity_a: [angle_1, angle_2, ...]
                        connectivity_b: [angle_1, angle_2, ...]
                    }
                }

            Where ``to_site`` is the index of a next nearest neighbor site,
            ``connectivity_a`` etc are the bonding connectivity type, e.g.
            ``'edge'`` or ``'corner'`` (for edge-sharing and corner-sharing
            connectivity), and ``angle_1`` etc are the bond angles as
            :obj:`float`.
        """
        connectivities: defaultdict[int, defaultdict[str, list[float]]] = defaultdict(
            lambda: defaultdict(list)
        )

        for nnn_site in self.get_next_nearest_neighbors(
            site_index, inc_inequivalent_site_index=True
        ):
            to_site = nnn_site["inequiv_index"]
            connectivity = nnn_site["connectivity"]
            connectivities[to_site][connectivity].extend(nnn_site["angles"])

        return defaultdict_to_dict(connectivities)

    def get_nnn_distance_summary(
        self, site_index: int
    ) -> dict[int, dict[str, list[float]]]:
        """Gets the next nearest neighbor distance summary for a site.

        Args:
            site_index: The site index (zero based).

        Returns:
            The connectivity distance data for the site, formatted as::

                {
                    to_site: {
                        connectivity_a: [distance_1, distance_2, ...]
                        connectivity_b: [distance_1, distance_2, ...]
                    }
                }

            Where ``to_site`` is the index of a next nearest neighbor site,
            ``connectivity_a`` etc are the bonding connectivity type, e.g.
            ``'edge'`` or ``'corner'`` (for edge-sharing and corner-sharing
            connectivity), and ``distance_1`` etc are the bond angles as
            :obj:`float`.
        """
        connectivities: defaultdict[int, defaultdict[str, list[float]]] = defaultdict(
            lambda: defaultdict(list)
        )

        for nnn_site in self.get_next_nearest_neighbors(
            site_index, inc_inequivalent_site_index=True
        ):
            to_site = nnn_site["inequiv_index"]
            connectivity = nnn_site["connectivity"]
            connectivities[to_site][connectivity].append(nnn_site["distance"])

        return defaultdict_to_dict(connectivities)

    def get_all_site_summaries(self):
        """Gets the site summaries for all sites.

        Returns:
            The site summaries for all sites, formatted as::

                {
                    site_index: site_summary
                }

            Where ``site_summary`` has the same format as produced by
            :meth:`SiteAnalyzer.get_site_summary`.
        """
        return {
            site: self.get_site_summary(site) for site in set(self.equivalent_sites)
        }

    def get_all_bond_distance_summaries(self) -> dict[int, dict[int, list[float]]]:
        """Gets the bond distance summaries for all sites.

        Returns:
            The bond distance summaries for all sites, formatted as::

                {
                    from_site: {
                        to_site: distances
                    }
                }

            Where ``from_site`` and ``to_site`` are site indices and
            ``distances`` is a :obj:`list` of :obj:`float` of bond distances.
        """
        return {
            from_site: self.get_bond_distance_summary(from_site)
            for from_site in set(self.equivalent_sites)
        }

    def get_all_connectivity_angle_summaries(
        self,
    ) -> dict[int, dict[int, dict[str, list[float]]]]:
        """Gets the connectivity angle summaries for all sites.

        The connectivity angles are the angles between a site and its
        next nearest neighbors.

        Returns:
            The connectivity angle summaries for all sites, formatted as::

                {
                    from_site: {
                        to_site: {
                            connectivity: angles
                        }
                    }
                }

            Where ``from_site`` and ``to_site`` are the site indices of
            two sites, ``connectivity`` is the connectivity type (e.g.
            ``'edge'`` or ``'face'``) and ``angles`` is a :obj:`list` of
            :obj:`float` of connectivity angles.
        """
        return {
            from_site: self.get_connectivity_angle_summary(from_site)
            for from_site in set(self.equivalent_sites)
        }

    def get_all_nnn_distance_summaries(
        self,
    ) -> dict[int, dict[int, dict[str, list[float]]]]:
        """Gets the next nearest neighbor distance summaries for all sites.

        Returns:
            The next nearest neighbor distance summaries for all sites,
            formatted as::

                {
                    from_site: {
                        to_site: {
                            connectivity: distances
                        }
                    }
                }

            Where ``from_site`` and ``to_site`` are the site indices of
            two sites, ``connectivity`` is the connectivity type (e.g.
            ``'edge'`` or ``'face'``) and ``distances`` is a :obj:`list` of
            :obj:`float` of distances.
        """
        return {
            from_site: self.get_nnn_distance_summary(from_site)
            for from_site in set(self.equivalent_sites)
        }

    def get_inequivalent_site_indices(self, site_indices: list[int]) -> list[int]:
        """Gets the inequivalent site indices from a list of site indices.

        Args:
            site_indices: The site indices.

        Returns:
            The inequivalent site indices. For example, if a structure has 4
            sites where the first two are equivalent and the last two are
            inequivalent. If  ``site_indices=[0, 1, 2, 3]`` the output will be::

                [0, 0, 2, 3]

        """
        return [self.equivalent_sites[i] for i in site_indices]

    def _calculate_equivalent_sites(
        self,
        likeness_tol: float = 0.001,
        bond_dist_tol: float = 0.01,
        bond_angle_tol: float = 0.1,
    ) -> list[int]:
        """Determines the indices of the structurally inequivalent sites.

        Args:
            likeness_tol: The tolerance used to determine if two likeness
                parameters are the same.
            bond_dist_tol: The tolerance used to determine if two bond lengths
                are the same.
            bond_angle_tol: The tolerance used to determine if two bond angles
                are the same.

        Two sites are considered equivalent if they are the same element, and
        have the same geometry and (next) nearest neighbors.

        Returns:
            A :obj:`list` of indices mapping each site in the structure to a
            structurally equivalent site. For example, if the first two sites
            are equivalent and the last two are both inequivalent, the data will
            be formatted as::

                [0, 0, 2, 3]

        """
        # TODO: Use site fingerprint rather than geometry type.
        inequiv_sites: dict[int, dict[str, Any]] = {}
        equivalent_sites: list[int] = []

        for site_index, site in enumerate(self.bonded_structure.structure):
            element = get_el_sp(site.specie)
            geometry = self.get_site_geometry(site_index)
            nn_sites = self.get_nearest_neighbors(
                site_index, inc_inequivalent_site_index=False
            )
            nnn_sites = self.get_next_nearest_neighbors(
                site_index, inc_inequivalent_site_index=False
            )

            matched = False
            for inequiv_index, inequiv_site in inequiv_sites.items():
                elem_match = element == inequiv_site["element"]
                geom_match = geometries_match(
                    geometry, inequiv_site["geometry"], likeness_tol=likeness_tol
                )
                nn_match = nn_summaries_match(
                    nn_sites, inequiv_site["nn_sites"], bond_dist_tol=bond_dist_tol
                )
                nnn_match = nnn_summaries_match(
                    nnn_sites, inequiv_site["nnn_sites"], bond_angle_tol=bond_angle_tol
                )

                if elem_match and geom_match and nn_match and nnn_match:
                    equivalent_sites.append(inequiv_index)
                    matched = True
                    break

            if not matched:
                # no matches therefore store original site index
                equivalent_sites.append(site_index)
                site_data = {
                    "element": element,
                    "geometry": geometry,
                    "nn_sites": nn_sites,
                    "nnn_sites": nnn_sites,
                }
                inequiv_sites[site_index] = site_data

        return equivalent_sites

    def _calculate_symmetry_labels(
        self, sym_equivalent_atoms: Sequence[int]
    ) -> list[int]:
        """Calculates the symmetry labels for all sites in the structure.

        The symmetry labels number the sites in the structure. If two sites
        are symmetrically equivalent they share the same symmetry label. The
        numbering begins at 1 for each element in the structure.

        Args:
            sym_equivalent_atoms: A :obj:`list` of indices mapping each site in
                the structure to a symmetrically equivalent site. The data
                should be formatted as given by the ``equivalent_atoms`` key in
                :meth`SpacegroupAnalyzer.get_symmetry_dataset()`.

        Returns:
            A mapping between the site index and symmetry label for that site.
        """
        symmetry_labels = dict()

        # this way is a little long winded but works if the sites aren't
        # grouped together by element
        for specie in self.bonded_structure.structure.species:
            el_indices = self.bonded_structure.structure.indices_from_symbol(
                get_el(specie)
            )
            equiv_indices = [sym_equivalent_atoms[x] for x in el_indices]

            count = 1
            equiv_index_to_sym_label: dict[int, int] = {}
            for el_index, equiv_index in zip(el_indices, equiv_indices):
                if equiv_index in equiv_index_to_sym_label:
                    symmetry_labels[el_index] = equiv_index_to_sym_label[equiv_index]
                else:
                    equiv_index_to_sym_label[equiv_index] = count
                    symmetry_labels[el_index] = count
                    count += 1

        return [symmetry_labels[i] for i in sorted(symmetry_labels.keys())]

    def _get_poly_formula(
        self,
        geometry: dict[str, Any],
        nn_sites: list[dict[str, Any]],
        nnn_sites: list[dict[str, Any]],
    ) -> str | None:
        """Gets the polyhedra formula of the nearest neighbor atoms.

        The polyhedral formula is effectively the sorted nearest neighbor
        atoms in a reduced format. For example, if the nearest neighbors are
        3 I atoms, 2 Br atoms and 1 Cl atom, the polyhedral formula will be
        "I3Br2Cl". The polyhedral formula will be ``None`` if the site geometry
        is not in :data:`robocrys.util.connected_geometries`.

        Args:
            geometry: The site geometry as produced by
                :meth:`SiteAnalyzer.get_site_geometry`.
            nn_sites: The nearest neighbor sites as produced by
                :meth:`SiteAnalyzer.get_nearest_neighbors`.
            nnn_sites: The next nearest neighbor sites as produced by
                :meth:`SiteAnalyzer.get_next_nearest_neighbors`.

        Returns:
            The polyhedral formula if the site geometry is in
            :data:`robocrys.util.connected_geometries` else ``None``.
        """

        def order_elements(el):
            if self.use_iupac_formula:
                return [get_el_sp(el).X, el]
            return [get_el_sp(el).iupac_ordering, el]

        nnn_geometries = [nnn_site["geometry"] for nnn_site in nnn_sites]

        poly_formula = ""
        if geometry["type"] in connected_geometries and any(
            nnn_geometry["type"] in connected_geometries
            for nnn_geometry in nnn_geometries
        ):
            nn_els = [get_el(nn_site["element"]) for nn_site in nn_sites]
            comp = Composition("".join(nn_els))
            el_amt_dict = comp.get_el_amt_dict()

            for e in sorted(el_amt_dict.keys(), key=order_elements):
                poly_formula += e
                poly_formula += str(formula_double_format(el_amt_dict[e]))

        return poly_formula or None


def geometries_match(
    geometry_a: dict[str, Any], geometry_b: dict[str, Any], likeness_tol: float = 0.001
) -> bool:
    """Determine whether two site geometries match.

    Geometry data should be formatted the same as produced by
    :meth:`robocrys.site.SiteAnalyzer.get_site_geometry`.

    Args:
        geometry_a: The first set of geometry data.
        geometry_b: The second set of geometry data.
        likeness_tol: The tolerance used to determine if two likeness parameters
            are the same.

    Returns:
        Whether the two geometries are the same.
    """
    return (
        geometry_a["type"] == geometry_b["type"]
        and abs(geometry_a["likeness"] - geometry_b["likeness"]) < likeness_tol
    )


def nn_summaries_match(
    nn_sites_a: list[dict[str, int | str]],
    nn_sites_b: list[dict[str, int | str]],
    bond_dist_tol: float = 0.01,
    match_bond_dists: bool = True,
) -> bool:
    """Determine whether two sets of nearest neighbors match.

    Nearest neighbor data should be formatted the same as produced by
    :meth:`robocrys.site.SiteAnalyzer.get_nearest_neighbors`.

    Args:
        nn_sites_a: The first set of nearest neighbors.
        nn_sites_b: The second set of nearest neighbors.
        bond_dist_tol: The tolerance used to determine if two bond lengths
            are the same.
        match_bond_dists: Whether to consider bond distances when matching.

    Returns:
        Whether the two sets of nearest neighbors match.
    """

    def nn_sites_order(nn_site):
        return [nn_site["element"], nn_site["dist"]]

    if len(nn_sites_a) != len(nn_sites_b):
        return False

    nn_sites_a = sorted(nn_sites_a, key=nn_sites_order)
    nn_sites_b = sorted(nn_sites_b, key=nn_sites_order)

    dists_match = [
        (
            abs(site_a["dist"] - site_b["dist"]) < bond_dist_tol  # type: ignore[operator]
            if match_bond_dists
            else True
        )
        for site_a, site_b in zip(nn_sites_a, nn_sites_b)
    ]
    elements_match = [
        site_a["element"] == site_b["element"]
        for site_a, site_b in zip(nn_sites_a, nn_sites_b)
    ]

    return all(d and e for d, e in zip(dists_match, elements_match))


def nnn_summaries_match(
    nnn_sites_a: list[dict[str, Any]],
    nnn_sites_b: list[dict[str, Any]],
    likeness_tol: float = 0.001,
    bond_angle_tol: float = 0.1,
    match_bond_angles: bool = True,
):
    """Determine whether two sets of next nearest neighbors match.

    Next nearest neighbor data should be formatted the same as produced by
    :meth:`robocrys.site.SiteAnalyzer.get_next_nearest_neighbors`.

    Args:
        nnn_sites_a: The first set of next nearest neighbors.
        nnn_sites_b: The second set of next nearest neighbors.
        likeness_tol: The tolerance used to determine if two likeness parameters
            are the same.
        bond_angle_tol: The tolerance used to determine if two bond angles
            are the same.
        match_bond_angles: Whether to consider bond angles when matching.

    Returns:
        Whether the two sets of next nearest neighbors match.
    """

    def nnn_sites_order(nnn_site):
        return [
            nnn_site["element"],
            nnn_site["geometry"]["type"],
            nnn_site["connectivity"],
            sorted(nnn_site["angles"]),
        ]

    if len(nnn_sites_a) != len(nnn_sites_b):
        return False

    nnn_sites_a = sorted(nnn_sites_a, key=nnn_sites_order)
    nnn_sites_b = sorted(nnn_sites_b, key=nnn_sites_order)

    elements_match = [
        site_a["element"] == site_b["element"]
        for site_a, site_b in zip(nnn_sites_a, nnn_sites_b)
    ]
    cons_match = [
        site_a["connectivity"] == site_b["connectivity"]
        for site_a, site_b in zip(nnn_sites_a, nnn_sites_b)
    ]
    geoms_match = [
        geometries_match(
            site_a["geometry"], site_b["geometry"], likeness_tol=likeness_tol
        )
        for site_a, site_b in zip(nnn_sites_a, nnn_sites_b)
    ]
    angles_match = [
        (
            all(
                abs(a_a - a_b) < bond_angle_tol
                for a_a, a_b in zip(sorted(site_a["angles"]), sorted(site_b["angles"]))
            )
            if match_bond_angles
            else True
        )
        for site_a, site_b in zip(nnn_sites_a, nnn_sites_b)
    ]

    return all(
        e and c and g and a
        for e, c, g, a in zip(elements_match, cons_match, geoms_match, angles_match)
    )
