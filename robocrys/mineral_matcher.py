"""
This module provides tools for matching structures to known mineral class.
"""

from itertools import islice
from pkg_resources import resource_filename

from pymatgen.analysis.aflow_prototypes import AflowPrototypeMatcher
from matminer.utils.io import load_dataframe_from_json

from robocrys.fingerprint import get_fingerprint, get_fingerprint_distance


class MineralMatcher(object):
    """Class to match a structure to a mineral name.

    Uses a precomputed database of minerals and their fingerprints, extracted
    from the AFLOW prototype database. For more information on this database
    see reference [aflow]_:

    .. [aflow] Mehl, M. J., Hicks, D., Toher, C., Levy, O., Hanson, R. M., Hart,
               G., & Curtarolo, S. (2017), The AFLOW library of crystallographic
               prototypes: part 1. Computational Materials Science, 136,
               S1-S828. doi: 10.1016/j.commatsci.2017.01.017
    """

    def __init__(self):
        db_file = resource_filename('robocrys', 'mineral_db.json.gz')
        self.mineral_db = load_dataframe_from_json(db_file)

    def get_best_mineral_name(self, structure):
        """Gets the "best" mineral name for a structure.

        Uses a combination of AFLOW prototype matching and fingerprinting to
        get the best guess. If a structure is not a perfect match but similar
        to a known mineral, "-like" will be added to the mineral name. If the
        structure is a good match to a mineral but contains a different number
        of element types than the mineral prototype, "-derived" will be added
        to the mineral name.

        The AFLOW structure prototypes are detailed in reference [aflow]_.

        Args:
            structure (Structure): A pymatgen `Structure` object to match.

        Return:
            (str or None): The mineral name as a string. If no match found,
            `None` will be returned
        """

        # The algorithm works as follows:
        #
        # 1. Check for AFLOW match. If single match return mineral name.
        # 2. If multiple matches, return the one with the smallest fingerprint
        #    distance.
        # 3. If no AFLOW match, get fingerprints within tolerance. If there are
        #    any matches, take the one with the smallest distance. Add "-like"
        #    to the mineral name.
        # 4. If no fingerprints within tolerance, check get fingerprints without
        #    constraining the number of species types. If any matches, take the
        #    best one and add "-derived" to the mineral name.

        self._set_distance_matrix(structure)  # pre-calculate distance matrix

        aflow_matches = self.get_aflow_matches(structure)

        if aflow_matches:
            return aflow_matches[0]['mineral']

        fingerprint_matches = self.get_fingerprint_matches(structure)
        if fingerprint_matches:
            return "{}-like".format(fingerprint_matches[0]['mineral'])

        fingerprint_derived = self.get_fingerprint_matches(
            structure, match_n_sp=False)
        if fingerprint_derived:
            return "{}-derived".format(fingerprint_derived[0]['mineral'])

        return None

    def get_aflow_matches(self, structure, initial_ltol=0.2, initial_stol=0.3,
                          initial_angle_tol=5.):
        """Gets minerals for a structure by matching to AFLOW prototypes.

        Overrides pymatgen `AflowPrototypeMatcher` to only return matches to
        prototypes with mineral names.

        Tolerance parameters are passed to a pymatgen `StructureMatcher` object.
        The tolerances are gradually decreased until only a single match is
        found (if possible).

        The AFLOW structure prototypes are detailed in reference [aflow]_.

        Args:
            structure (Structure): A pymatgen `Structure` object to match.
            initial_ltol (float, optional): The fractional length tolerance.
            initial_stol (float, optional): The site coordinate tolerance.
            initial_angle_tol (float, optional): The angle tolerance.

        Returns:
            (list or None): A list of dictionaries, sorted by how close the
            match is, with the keys 'mineral', 'distance', 'structure'.
            Distance is the euclidean distance between the structure and
            prototype fingerprints. If no match was found within the tolerances,
            `None` will be returned.
        """
        self._set_distance_matrix(structure)

        # redefine AflowPrototypeMatcher._match_prototype function to run over
        # our custom pandas DataFrame of AFLOW prototypes. This DataFrame only
        # contains entries from the AFLOW database with mineral names. We
        # have also pre-calculated the fingerprints and distances to make this
        # quicker.
        def _match_prototype(structure_matcher, s):
            tags = []
            for index, row in self.mineral_db_.iterrows():
                p = row['structure']
                m = structure_matcher.fit_anonymous(p, s)
                if m:
                    tags.append(_get_row_data(row))
            return tags

        matcher = AflowPrototypeMatcher(initial_ltol=initial_ltol,
                                        initial_stol=initial_stol,
                                        initial_angle_tol=initial_angle_tol)
        matcher._match_prototype = _match_prototype

        return matcher.get_prototypes(structure)

    def get_fingerprint_matches(self, structure, distance_cutoff=0.4,
                                max_n_matches=None, match_n_sp=True):
        """Gets minerals for a structure by matching to AFLOW fingerprints.

        Only AFLOW prototypes with mineral names are considered. The AFLOW
        structure prototypes are detailed in reference [aflow]_.

        Args:
            structure (Structure): A structure to match.
            distance_cutoff (float, optional): Cutoff to determine how similar a
                match must be to be returned. The distance is measured between
                the structural fingerprints in euclidean space.
            max_n_matches (int, optional): Maximum number of matches to return.
                Set to `None` to return all matches within the cutoff.
            match_n_sp (bool, optional): Whether the structure and mineral must
                have the same number of species. Defaults to True.

        Returns:
            (list or None): A list of dictionaries, sorted by how close the
            match is, with the keys 'mineral', 'distance', 'structure'.
            Distance is the euclidean distance between the structure and
            prototype fingerprints. If no match was found within the tolerances,
            `None` will be returned.
        """
        self._set_distance_matrix(structure)

        mineral_db = self.mineral_db_

        if match_n_sp:
            ntypesp = structure.ntypesp
            mineral_db = mineral_db[mineral_db['ntypesp'] == ntypesp]

        num_rows = mineral_db.shape[0]
        max_n_matches = max_n_matches if max_n_matches else num_rows
        max_n_matches = num_rows if max_n_matches > num_rows else max_n_matches

        minerals = [_get_row_data(row)
                    for i, row in islice(mineral_db.iterrows(), max_n_matches)
                    if row['distance'] < distance_cutoff]

        return minerals if minerals else None

    def _set_distance_matrix(self, structure):
        """Utility function to calculate distance between structure and minerals.

        First checks to see if the distances have already been calculated for
        the structure. If not, the distances are stored in a class variable
        for use by other class methods.

        Args:
            structure (structure): A pymatgen `Structure` object.
        """
        if (hasattr(self, 'structure_') and self.structure_ == structure and
                hasattr(self, 'mineral_db_')):
            return

        data = self.mineral_db.copy()
        fingerprint = get_fingerprint(structure)

        data['distance'] = data['fingerprint'].apply(
            lambda x: get_fingerprint_distance(x, fingerprint))

        self.mineral_db_ = data.sort_values(by='distance')
        self.structure_ = structure


def _get_row_data(row):
    """Utility function to extract mineral data from pandas `DataFrame` row."""
    return {'mineral': row['mineral'], 'distance': row['distance'],
            'structure': row['structure']}
