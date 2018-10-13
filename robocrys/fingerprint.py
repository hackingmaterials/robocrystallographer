import numpy as np

from pymatgen.core.structure import IStructure
from matminer.featurizers.structure import SiteStatsFingerprint


def get_fingerprint(structure, fingerprint_preset='CrystalNNFingerprint_ops',
                    stats=('mean', 'std_dev')):
    """Gets the fingerprint for a structure.

    Args:
        structure (Structure): A pymatgen `Structure` object.
        fingerprint_preset (str): The preset to use when calculating the
            fingerprint. See the `SiteStatsFingerprint` class in matminer for
            more details.
        stats (str or iterable): The stats to include in fingerprint. See
            the `SiteStatsFingerprint` class in matminer for more details.

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
        fingerprint_a = get_fingerprint(structure_a)
    else:
        fingerprint_a = np.array(structure_a)

    if issubclass(type(structure_b), IStructure):
        fingerprint_b = get_fingerprint(structure_b)
    else:
        fingerprint_b = np.array(structure_b)

    dist = np.linalg.norm(fingerprint_a - fingerprint_b)
    return dist
