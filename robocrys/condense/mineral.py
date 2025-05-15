"""This module provides tools for matching structures to known mineral class."""

from __future__ import annotations

from importlib.resources import files as import_resource_file
from itertools import islice
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from pymatgen.analysis.prototypes import AflowPrototypeMatcher
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure
from robocrys.condense.fingerprint import (
    get_fingerprint_distance,
    get_structure_fingerprint,
)


if TYPE_CHECKING:
    from typing import Any
    import pandas as pd

    from pymatgen.core.structure import Structure


_mineral_db_file = str(import_resource_file("robocrys.condense") / "mineral_db.json.gz")


def load_mineral_db_file(min_db_path: str | Path = _mineral_db_file) -> pd.DataFrame:
    """Define protocol to load mineral database.

    Args:
        min_db_path : str or Path
            Path to the Mineral DB
    """
    mineral_db = pd.read_json(min_db_path, compression="gzip", orient="split")
    mineral_db.structure = mineral_db.structure.apply(Structure.from_dict)
    return mineral_db


class MineralMatcher:
    """Class to match a structure to a mineral name.

    Uses a precomputed database of minerals and their fingerprints, extracted
    from the AFLOW prototype database. For more information on this database
    see reference [aflow]_:

    .. [aflow] Mehl, M. J., Hicks, D., Toher, C., Levy, O., Hanson, R. M., Hart,
               G., & Curtarolo, S. (2017), The AFLOW library of crystallographic
               prototypes: part 1. Computational Materials Science, 136,
               S1-S828. doi: 10.1016/j.commatsci.2017.01.017

    Args:
        initial_ltol: The fractional length tolerance used in the AFLOW
            structure matching.
        initial_stol : The site coordinate tolerance used in the AFLOW
            structure matching.
        initial_angle_tol: The angle tolerance used in the AFLOW structure
            matching.
        use_fingerprint_matching: Whether to use the fingerprint distance to
            match minerals.
        fingerprint_distance_cutoff: Cutoff to determine how similar a match
            must be to be returned. The distance is measured between the
            structural fingerprints in euclidean space.
        mineral_db : Optional path or pandas .DataFrame object containing the
            mineral fingerprint database.
    """

    def __init__(
        self,
        initial_ltol: float = 0.2,
        initial_stol: float = 0.3,
        initial_angle_tol: float = 5.0,
        use_fingerprint_matching: bool = True,
        fingerprint_distance_cutoff: float = 0.4,
        mineral_db: str | Path | pd.DataFrame | None = None,
    ):
        mineral_db = mineral_db or _mineral_db_file
        if isinstance(mineral_db, (str, Path)):
            self.mineral_db = load_mineral_db_file(mineral_db)
        else:
            self.mineral_db = mineral_db

        self.initial_ltol = initial_ltol
        self.initial_stol = initial_stol
        self.initial_angle_tol = initial_angle_tol
        self.fingerprint_distance_cutoff = fingerprint_distance_cutoff
        self.use_fingerprint_matching = use_fingerprint_matching
        self._structure: Structure | None = None
        self._mineral_db: pd.DataFrame | None = None

    def get_best_mineral_name(self, structure: Structure) -> dict[str, Any]:
        """Gets the "best" mineral name for a structure.

        Uses a combination of AFLOW prototype matching and fingerprinting to
        get the best mineral name.

        The AFLOW structure prototypes are detailed in reference [aflow]_.

        The algorithm works as follows:

        1. Check for AFLOW match. If single match return mineral name.
        2. If multiple matches, return the one with the smallest fingerprint
           distance.
        3. If no AFLOW match, get fingerprints within tolerance. If there are
           any matches, take the one with the smallest distance.
        4. If no fingerprints within tolerance, check get fingerprints without
           constraining the number of species types. If any matches, take the
           best one.

        Args:
            structure (Structure): A pymatgen ``Structure`` object to match.

        Return:
            (dict): The mineral name information. Stored as a dict with the keys
            "type", "distance", "n_species_types_match", corresponding to the
            mineral name, the fingerprint distance between the prototype and
            known mineral, and whether the number of species types in the
            structure matches the number in the known prototype, respectively.
            If no mineral match is determined, the mineral type will be
            ``None``. If an AFLOW match is found, the distance will be set to
            -1.
        """
        self._set_distance_matrix(structure)  # pre-calculate distance matrix

        aflow_matches = self.get_aflow_matches(structure)

        fingerprint_matches = self.get_fingerprint_matches(structure)
        fingerprint_derived = self.get_fingerprint_matches(structure, match_n_sp=False)

        distance = -1
        n_species_types_match = True
        if aflow_matches:
            # mineral db sorted by fingerprint distance so first result always
            # has a smaller distance
            mineral = aflow_matches[0]["type"]

        elif fingerprint_matches and self.use_fingerprint_matching:
            mineral = fingerprint_matches[0]["type"]
            distance = fingerprint_matches[0]["distance"]

        elif fingerprint_derived and self.use_fingerprint_matching:
            mineral = fingerprint_derived[0]["type"]
            distance = fingerprint_derived[0]["distance"]
            n_species_types_match = False

        else:
            mineral = None

        return {
            "type": mineral,
            "distance": distance,
            "n_species_type_match": n_species_types_match,
        }

    def get_aflow_matches(
        self,
        structure: Structure,
    ) -> list[dict[str, Any]] | None:
        """Gets minerals for a structure by matching to AFLOW prototypes.

        Overrides
        :class:`pymatgen.analysis.prototypes.AflowPrototypeMatcher` to
        only return matches to prototypes with known mineral names.

        The AFLOW tolerance parameters (defined in the init method) are passed
        to a :class:`pymatgen.analysis.structure_matcher.StructureMatcher`
        object. The tolerances are gradually decreased until only a single match
        is found (if possible).

        The AFLOW structure prototypes are detailed in reference [aflow]_.

        Args:
            structure: A pymatgen structure to match.

        Returns:
            A :obj:`list` of :obj:`dict`, sorted by how close the match is, with
            the keys 'type', 'distance', 'structure'. Distance is the
            euclidean distance between the structure and prototype fingerprints.
            If no match was found within the tolerances, ``None`` will be
            returned.
        """
        self._set_distance_matrix(structure)

        matcher = ReducedAflowPrototypeMatcher(
            initial_ltol=self.initial_ltol,
            initial_stol=self.initial_stol,
            initial_angle_tol=self.initial_angle_tol,
            mineral_db=self._mineral_db,
        )

        return matcher.get_prototypes(structure)

    def get_fingerprint_matches(
        self,
        structure: Structure,
        max_n_matches: int | None = None,
        match_n_sp: bool = True,
        mineral_name_constraint: str | None = None,
    ) -> list[dict[str, Any]] | None:
        """Gets minerals for a structure by matching to AFLOW fingerprints.

        Only AFLOW prototypes with mineral names are considered. The AFLOW
        structure prototypes are detailed in reference [aflow]_.

        Args:
            structure: A structure to match.
            max_n_matches: Maximum number of matches to return. Set to ``None``
                to return all matches within the cutoff.
            match_n_sp: Whether the structure and mineral must have the same
                number of species. Defaults to True.
            mineral_name_constraint: Whether to limit the matching to a specific
                mineral name.

        Returns:
            A :obj:`list` of :obj:`dict`, sorted by how close the match is, with
            the keys 'type', 'distance', 'structure'. Distance is the
            euclidean distance between the structure and prototype fingerprints.
            If no match was found within the tolerances, ``None`` will be
            returned.
        """
        self._set_distance_matrix(structure)

        mineral_db: pd.DataFrame = self._mineral_db  # type: ignore[assignment]

        if mineral_name_constraint:
            mineral_db = mineral_db[
                mineral_db["mineral"].str.lower() == mineral_name_constraint
            ]

        if match_n_sp:
            mineral_db = mineral_db[mineral_db["n_elems"] == structure.n_elems]

        num_rows = mineral_db.shape[0]
        max_n_matches = max_n_matches or num_rows
        max_n_matches = min(max_n_matches, num_rows)

        minerals = [
            _get_row_data(row)
            for i, row in islice(mineral_db.iterrows(), max_n_matches)
            if row["distance"] < self.fingerprint_distance_cutoff
        ]

        return minerals if minerals else None

    def _set_distance_matrix(self, structure: Structure) -> None:
        """Utility func to calculate distance between structure and minerals.

        First checks to see if the distances have already been calculated for
        the structure. If not, the distances are stored in a class variable
        for use by other class methods.

        Args:
            structure: A structure.
        """
        if self._structure == structure and self._mineral_db is not None:
            return

        data = self.mineral_db.copy()
        fingerprint = get_structure_fingerprint(structure)

        if np.linalg.norm(fingerprint) < 0.4:
            # fingerprint is too small for a reasonable match, indicates very
            # little bonding or small order parameter matches
            fingerprint = get_structure_fingerprint(structure, prototype_match=False)

        data["distance"] = data["fingerprint"].apply(
            lambda x: get_fingerprint_distance(x, fingerprint)
        )

        self._mineral_db = data.sort_values(by="distance")
        self._structure = structure


def _get_row_data(row: pd.Series) -> dict[str, Any]:
    """Utility function to extract mineral data from pandas `DataFrame` row."""
    return {
        "type": row["mineral"],
        "distance": row["distance"],
        "structure": row["structure"],
    }


class ReducedAflowPrototypeMatcher(AflowPrototypeMatcher):
    """Match structures to predefined library of reduced, named AFLOW prototypes.

    Based on `pymatgen.analysis.prototypes.AflowPrototypeMatcher`,
    which uses the methods of:
        Mehl, M. J., Hicks, D., Toher, C., Levy, O., Hanson, R. M.,
            Hart, G., & Curtarolo, S. (2017).
        The AFLOW library of crystallographic prototypes: part 1.
        Computational Materials Science, 136, S1-S828.
        https://doi.org/10.1016/j.commatsci.2017.01.017
    """

    def __init__(
        self,
        initial_ltol: float = 0.2,
        initial_stol: float = 0.3,
        initial_angle_tol: float = 5,
        mineral_db: pd.DataFrame | None = None,
    ) -> None:
        """
        Tolerances as defined in StructureMatcher. Tolerances will be
        gradually decreased until only a single match is found (if possible).

        Args:
            initial_ltol (float): fractional length tolerance.
            initial_stol (float): site tolerance.
            initial_angle_tol (float): angle tolerance.
            mineral_db (pd.DataFrame or None) : precomputed database of prototypes
                and their fingerprints.
        """
        self.initial_ltol = initial_ltol
        self.initial_stol = initial_stol
        self.initial_angle_tol = initial_angle_tol
        if mineral_db is None:
            self._mineral_db = load_mineral_db_file()
        else:
            self._mineral_db = mineral_db

    def _match_prototype(
        self, structure_matcher: StructureMatcher, s: Structure
    ) -> list[dict[str, Any]]:
        """Match a structure to the predefined table of reduced structures.

        Args:
            structure_matcher (StructureMatcher)
            s (Structure)
        """
        tags = []
        for _, row in self._mineral_db.iterrows():
            p = row["structure"]
            m = structure_matcher.fit_anonymous(p, s)
            if m:
                tags.append(_get_row_data(row))
        return tags
