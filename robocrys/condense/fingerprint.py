from typing import Union, List, Dict, Optional, Tuple, Iterable

import numpy as np
from matminer.featurizers.site import CrystalNNFingerprint
from matminer.featurizers.structure import SiteStatsFingerprint

from pymatgen.core.structure import IStructure


def get_site_fingerprints(structure: IStructure,
                          as_dict: bool = True,
                          preset: str = 'CrystalNNFingerprint_ops'
                          ) -> Union[List[Dict[str, int]], np.ndarray]:
    """Gets the fingerprint for all sites in a structure.

    Args:
        structure: A structure.
        as_dict: Whether to return the fingerprints as a dictionary of
            ``{'op': val}``. Defaults to ``True``.
        preset: The preset to use when calculating the fingerprint. See
            :class:`matminer.featurizers.structure.SiteStatsFingerprint``
            for more details.

    Returns:
        The fingerprint for all sites in the structure. If ``as_dict == True``,
        the data will be returned as a :obj:`list` of :obj:`dict` containing the
        order parameters as::

            [{'op': val}]

        for each site. If ``as_dict == False``, the data will be returned as a
        :class:`numoy.ndarray` containing the fingerprint for each site as::

            [site_index][op_index]
    """
    ssf = SiteStatsFingerprint.from_preset(preset, stats=None)

    # transpose fingerprint from [op_type][site] to [site][op_type]
    site_fingerprints = np.array(ssf.featurize(structure)).T

    if as_dict:
        labels = ssf.feature_labels()
        site_fingerprints = [dict(zip(labels, x)) for x in site_fingerprints]

    return site_fingerprints


def get_structure_fingerprint(structure: IStructure,
                              preset: str = 'CrystalNNFingerprint_ops',
                              stats: Optional[Tuple[str]] = ('mean', 'std_dev'),
                              use_distance_cutoffs: bool = True
                              ) -> np.ndarray:
    """Gets the fingerprint for a structure.

    Args:
        structure: A structure.
        preset: The preset to use when calculating the fingerprint. See
            :class:`matminer.featurizers.structure.SiteStatsFingerprint``
            for more details.
        stats: The stats to include in fingerprint. See
            :class:`matminer.featurizers.structure.SiteStatsFingerprint``
            for more details.
        use_distance_cutoffs: Whether to use distance cutoffs when calculating
            the structure fingerprint.

    Returns:
        The structure fingerprint as a :class:`numpy.ndarray`.
    """
    # TODO: Add distance_cutoff option to matminer so we can user preset arg
    # currently don't use SiteStatsFingerprint.from_preset as we need to pass in
    # distance_cutoffs param
    if use_distance_cutoffs:
        ssf = SiteStatsFingerprint(
            CrystalNNFingerprint.from_preset("ops", cation_anion=False),
            stats=stats)
    else:
        ssf = SiteStatsFingerprint(
            CrystalNNFingerprint.from_preset("ops", cation_anion=False,
                                             distance_cutoffs=None),
            stats=stats)
    return np.array(ssf.featurize(structure))


def get_fingerprint_distance(structure_a: Union[IStructure, Iterable],
                             structure_b: Union[IStructure, Iterable]
                             ) -> float:
    """Gets the euclidean distance between the fingerprints of two structures.

    Args:
        structure_a: The first structure or fingerprint. Can be provided as a
            structure or a fingerprint. If provided as a structure, the
            fingerprint will be calculated first, so generally it is quicker
            to pre-calculate the fingerprint if comparing against multiple
            structures in turn.
        structure_b: The second structure or fingerprint. Can be provided as a
            structure or a fingerprint. If provided as a structure, the
            fingerprint will be calculated first, so generally it is quicker
            to pre-calculate the fingerprint if comparing against multiple
            structures in turn.

    Returns:
        The euclidean distance between fingerprints as a :class:`numpy.ndarray`.
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
