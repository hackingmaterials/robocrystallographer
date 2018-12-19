"""
This module provides tools for matching structures to known mineral class.
"""

from itertools import islice
from typing import List, Optional, Dict, Text, Any

import numpy as np
from matminer.utils.io import load_dataframe_from_json
from pkg_resources import resource_filename

from pymatgen.analysis.aflow_prototypes import AflowPrototypeMatcher
from pymatgen.core.structure import IStructure
from robocrys.condense.fingerprint import (get_structure_fingerprint,
                                           get_fingerprint_distance)


class MineralMatcher(object):
    """Class to match a structure to a mineral name.

    Uses a precomputed database of minerals and their fingerprints, extracted
    from the AFLOW prototype database. For more information on this database
    see reference [aflow]_:

    .. [aflow] Mehl, M. J., Hicks, D., Toher, C., Levy, O., Hanson, R. M., Hart,
               G., & Curtarolo, S. (2017), The AFLOW library of crystallographic
               prototypes: part 1. Computational Materials Science, 136,
               S1-S828. doi: 10.1016/j.commatsci.2017.01.017

    Args:
        aflow_initial_ltol: The fractional length tolerance used in the AFLOW
            structure matching.
        aflow_initial_stol : The site coordinate tolerance used in the AFLOW
            structure matching.
        aflow_initial_angle_tol: The angle tolerance used in the AFLOW structure
            matching.
        use_fingerprint_matching: Whether to use the fingerprint distance to
            match minerals.
        fingerprint_distance_cutoff: Cutoff to determine how similar a match
            must be to be returned. The distance is measured between the
            structural fingerprints in euclidean space.
    """

    def __init__(self,
                 aflow_initial_ltol: float = 0.2,
                 aflow_initial_stol: float = 0.3,
                 aflow_initial_angle_tol: float = 5.,
                 use_fingerprint_matching: bool = True,
                 fingerprint_distance_cutoff: float = 0.4):
        db_file = resource_filename('robocrys.condense', 'mineral_db.json.gz')
        self.mineral_db = load_dataframe_from_json(db_file)
        self.aflow_initial_ltol = aflow_initial_ltol
        self.aflow_initial_stol = aflow_initial_stol
        self.aflow_initial_angle_tol = aflow_initial_angle_tol
        self.fingerprint_distance_cutoff = fingerprint_distance_cutoff
        self.use_fingerprint_matching = use_fingerprint_matching
        self._structure = None
        self._mineral_db = None

    def get_best_mineral_name(self, structure: IStructure) -> Dict[Text, Any]:
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
        fingerprint_derived = self.get_fingerprint_matches(
            structure, match_n_sp=False)

        distance = -1
        n_species_types_match = True
        if aflow_matches:
            # mineral db sorted by fingerprint distance so first result always
            # has a smaller distance
            mineral = aflow_matches[0]['type']

        elif fingerprint_matches and self.use_fingerprint_matching:
            mineral = fingerprint_matches[0]['type']
            distance = fingerprint_matches[0]['distance']

        elif fingerprint_derived and self.use_fingerprint_matching:
            mineral = fingerprint_derived[0]['type']
            distance = fingerprint_derived[0]['distance']
            n_species_types_match = False

        else:
            mineral = None

        return {'type': mineral, 'distance': distance,
                'n_species_type_match': n_species_types_match}

    def get_aflow_matches(self, structure: IStructure,
                          ) -> Optional[List[Dict[Text, Any]]]:
        """Gets minerals for a structure by matching to AFLOW prototypes.

        Overrides
        :class:`pymatgen.analysis.aflow_prototypes.AflowPrototypeMatcher` to
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

        # redefine AflowPrototypeMatcher._match_prototype function to run over
        # our custom pandas DataFrame of AFLOW prototypes. This DataFrame only
        # contains entries from the AFLOW database with mineral names. We
        # have also pre-calculated the fingerprints and distances to make this
        # quicker.
        def _match_prototype(structure_matcher, s):
            tags = []
            for _, row in self._mineral_db.iterrows():
                p = row['structure']
                m = structure_matcher.fit_anonymous(p, s)
                if m:
                    tags.append(_get_row_data(row))
            return tags

        matcher = AflowPrototypeMatcher(
            initial_ltol=self.aflow_initial_ltol,
            initial_stol=self.aflow_initial_stol,
            initial_angle_tol=self.aflow_initial_angle_tol)
        matcher._match_prototype = _match_prototype

        return matcher.get_prototypes(structure)

    def get_fingerprint_matches(self,
                                structure: IStructure,
                                max_n_matches: Optional[int] = None,
                                match_n_sp: bool = True
                                ) -> Optional[List[Dict[Text, Any]]]:
        """Gets minerals for a structure by matching to AFLOW fingerprints.

        Only AFLOW prototypes with mineral names are considered. The AFLOW
        structure prototypes are detailed in reference [aflow]_.

        Args:
            structure: A structure to match.
            max_n_matches: Maximum number of matches to return. Set to ``None``
                to return all matches within the cutoff.
            match_n_sp: Whether the structure and mineral must have the same
                number of species. Defaults to True.

        Returns:
            A :obj:`list` of :obj:`dict`, sorted by how close the match is, with
            the keys 'type', 'distance', 'structure'. Distance is the
            euclidean distance between the structure and prototype fingerprints.
            If no match was found within the tolerances, ``None`` will be
            returned.
        """
        self._set_distance_matrix(structure)

        mineral_db = self._mineral_db

        if match_n_sp:
            ntypesp = structure.ntypesp
            mineral_db = mineral_db[mineral_db['ntypesp'] == ntypesp]

        num_rows = mineral_db.shape[0]
        max_n_matches = max_n_matches if max_n_matches else num_rows
        max_n_matches = num_rows if max_n_matches > num_rows else max_n_matches

        minerals = [_get_row_data(row)
                    for i, row in islice(mineral_db.iterrows(), max_n_matches)
                    if row['distance'] < self.fingerprint_distance_cutoff]

        return minerals if minerals else None

    def _set_distance_matrix(self, structure: IStructure):
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
            fingerprint = get_structure_fingerprint(
                structure, use_distance_cutoffs=False)

        data['distance'] = data['fingerprint'].apply(
            lambda x: get_fingerprint_distance(x, fingerprint))

        self._mineral_db = data.sort_values(by='distance')
        self._structure = structure


def _get_row_data(row: Dict) -> Dict[Text, Any]:
    """Utility function to extract mineral data from pandas `DataFrame` row."""
    return {'type': row['mineral'], 'distance': row['distance'],
            'structure': row['structure']}
