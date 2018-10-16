import numpy as np

from pymatgen.core.structure import IStructure
from matminer.featurizers.structure import SiteStatsFingerprint


def get_site_fingerprints(structure, as_dict=True,
                          fingerprint_preset='CrystalNNFingerprint_ops'):
    """Gets the fingerprint for all sites in a structure.

    Args:
        structure (Structure): A pymatgen Structure object.
        as_dict (bool, optional): Whether to return the fingerprints as a
            dictionary of `{'op': val}`. Defaults to `True`.
        fingerprint_preset (str, optional): The preset to use when calculating
            the fingerprint. See the `SiteStatsFingerprint` class in matminer
            for more details.

    Returns:
        (list of dict or `np.ndarray`): The fingerprint for all sites in the
        structure. If `as_dict == True`, the data will be returned as a list
        of dictionaries containing the order parameters as::

            [{'op': val}]

        for each site. If `as_dict == False`, the data will be returned as a
        numpy array, containing the fingerprint for each site as::

            [site_id][op_id]
    """
    ssf = SiteStatsFingerprint.from_preset(fingerprint_preset, stats=None)

    # transpose fingerprint from [op_type][site] to [site][op_type]
    site_fingerprints = np.array(ssf.featurize(structure)).T

    if as_dict:
        labels = ssf.feature_labels()
        site_fingerprints = [dict(zip(labels, x)) for x in site_fingerprints]

    return site_fingerprints


def get_structure_fingerprint(structure,
                              fingerprint_preset='CrystalNNFingerprint_ops',
                              stats=('mean', 'std_dev')):
    """Gets the fingerprint for a structure.

    Args:
        structure (Structure): A pymatgen `Structure` object.
        fingerprint_preset (str, optional): The preset to use when calculating
            the fingerprint. See the `SiteStatsFingerprint` class in matminer
            for more details.
        stats (str or iterable, optional): The stats to include in fingerprint.
            See the `SiteStatsFingerprint` class in matminer for more details.

    Returns:
        (numpy.ndarray): The structure fingerprint as a `numpy.ndarray`.
    """
    ssf = SiteStatsFingerprint.from_preset(fingerprint_preset, stats=stats)
    return np.array(ssf.featurize(structure))


def get_fingerprint_distance(structure_a, structure_b):
    """Gets the euclidean distance between the fingerprints of two structures.

    Args:
        structure_a (Structure or list-like): The first structure or
            fingerprint. Can be provided as a pymatgen `Structure` object or a
            fingerprint as a list, tuple or `numpy.ndarray`. If provided as a
            `Structure`, the fingerprint will be calculated first, so generally
            it is quicker to pre-calculate the fingerprint if comparing against
            multiple structures in turn.
        structure_b (Structure or list-like): The second structure or
            fingerprint. Can be provided as a pymatgen `Structure` object or a
            fingerprint as a list, tuple or `numpy.ndarray`. If provided as a
            `Structure`, the fingerprint will be calculated first, so generally
            it is quicker to pre-calculate the fingerprint if comparing against
            many structures in turn.

    Returns:
        (numpy.ndarray): The euclidean distance between fingerprints as a
        `numpy.ndarray`.
    """
    if issubclass(type(structure_a), IStructure):
        fingerprint_a = get_structure_fingerprint(structure_a)
    else:
        fingerprint_a = np.array(structure_a)

    if issubclass(type(structure_b), IStructure):
        fingerprint_b = get_structure_fingerprint(structure_b)
    else:
        fingerprint_b = np.array(structure_b)

    dist = np.linalg.norm(fingerprint_a - fingerprint_b)
    return dist
