"""This module contains a class to obtain robocrystallographer ML features."""

from __future__ import annotations

from itertools import product
from typing import TYPE_CHECKING

from matminer.featurizers.base import BaseFeaturizer  # type: ignore[import-untyped]
from numpy import mean
from pymatgen.analysis.local_env import cn_opt_params
from pymatgen.core.structure import Structure

from robocrys import StructureCondenser
from robocrys.featurize.adapter import FeaturizerAdapter
from robocrys.util import connected_geometries

if TYPE_CHECKING:
    pass

_geometries = [geometry for cn in cn_opt_params.values() for geometry in cn]
_dimensionalities = (3, 2, 1, 0)
_dimensionality_sets = (
    (3, 2, 1, 0),
    (3, 2, 1),
    (3, 2, 0),
    (3, 1, 0),
    (3, 2),
    (3, 1),
    (3, 0),
    (2, 1, 0),
    (2, 1),
    (2, 0),
    (1, 0),
)
_molecules = ["water", "oxygen", "ammonia", "methane"]
_connectivities = ["corner", "edge", "face"]
_cns = range(1, 13)


class RobocrysFeaturizer(BaseFeaturizer):
    """Class to generate structure features from robocrystallographer output.

    Args:
        condenser_kwargs: Keyword arguments that will be passed to
            :obj:`robocrys.condense.StructureCondenser`.
        distorted_tol: The value under which the site geometry will be
            classified as distorted.
    """

    def __init__(
        self, condenser_kwargs: dict | None = None, distorted_tol: float = 0.6
    ):
        condenser_kwargs = condenser_kwargs if condenser_kwargs else {}
        self._sc = StructureCondenser(**condenser_kwargs)
        self._distorted_tol = distorted_tol

    def featurize(self, s: Structure) -> list[float | bool | str | None]:
        """Featurizes a structure using robocrystallographer.

        Args:
            s: A structure.

        Returns:
            The robocrystallographer features.
        """
        fa = FeaturizerAdapter(
            self._sc.condense_structure(s), distorted_tol=self._distorted_tol
        )

        # add general structure features
        features: list[float | bool | str | None] = [
            fa.mineral["type"],
            fa.spg_symbol,
            fa.crystal_system,
            fa.dimensionality,
            fa.is_vdw_heterostructure,
            fa.is_interpenetrated,
            fa.is_intercalated,
        ]

        # add dimensionality features
        features += [fa.is_dimensionality(d) for d in _dimensionalities]
        features += [fa.is_dimensionality(d) for d in _dimensionality_sets]
        features += [d in fa.component_dimensionalities for d in _dimensionalities]

        # add molecule features
        features += [fa.contains_named_molecule]
        features += [fa.contains_molecule(m) for m in _molecules]

        # add geometry features
        features += [fa.contains_geometry_type(g) for g in _geometries]
        features += [fa.contains_geometry_type(g, distorted=True) for g in _geometries]
        features += [
            fa.average_coordination_number,
            fa.average_cation_coordination_number,
            fa.average_anion_coordination_number,
        ]

        # add polyhedral features
        features += [
            fa.contains_polyhedra,
            fa.contains_corner_sharing_polyhedra,
            fa.contains_edge_sharing_polyhedra,
            fa.contains_face_sharing_polyhedra,
        ]

        # add connectivity features
        features += [
            fa.contains_connected_geometry(c, g)
            for c, g in product(_connectivities, connected_geometries)
        ]
        features += [fa.average_corner_sharing_octahedral_tilt_angle]

        # add fractional features
        features += [fa.frac_sites_polyhedra]
        features += [fa.frac_site_geometry(g) for g in _geometries]
        features += [fa.frac_sites_n_coordinate(n) for n in _cns]

        all_distances = fa.all_bond_lengths()
        # add bond length features
        features += [max(all_distances), min(all_distances), mean(all_distances)]

        return features

    def feature_labels(self):
        # general features
        labels = [
            "mineral_prototype",
            "spg_symbol",
            "crystal_system",
            "dimensionality",
            "is_vdw_heterostructure",
            "is_interpenetrated",
            "is_intercalated",
        ]

        # dimensionality features
        labels += [f"is_only_{d}d" for d in _dimensionalities]
        labels += [
            "is_{}".format("_".join([f"{d}d" for d in ds]))
            for ds in _dimensionality_sets
        ]
        labels += [f"contains_{d}d_component" for d in _dimensionalities]

        # molecule features
        labels += ["contains_named_molecule"]
        labels += [f"contains_{m}" for m in _molecules]

        # geometry features
        labels += [f"contains_{g}" for g in _geometries]
        labels += [f"contains_distorted_{g}" for g in _geometries]
        labels += ["average_site_cn", "average_cation_cn", "average_anion_cn"]

        # polyhedral features
        labels += [
            "contains_polyhedra",
            "contains_corner_sharing_polyhedra",
            "contains_edge_sharing_polyhedra",
            "contains_face_sharing_polyhedra",
        ]

        # connectivity features
        labels += [
            f"contains_{c}_{g}"
            for c, g in product(_connectivities, connected_geometries)
        ]
        labels += ["corner_sharing_octahedral_tilt_angle"]

        # fractional features
        labels += ["frac_site_polyhedra"]
        labels += [f"frac_sites_{g}" for g in _geometries]
        labels += [f"frac_sites_{n}_coordinate" for n in _cns]

        # bond length features
        labels += ["max_bond_length", "min_bond_length", "average_bond_length"]

        return labels

    def citations(self):
        return ["in prep."]

    def implementors(self):
        return ["Alex Ganose"]
