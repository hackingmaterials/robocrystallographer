"""
This module implements tools for calculating the statistics of structure
features across many structures.
"""

from collections import Counter, defaultdict
from statistics import mean
from typing import Dict, Iterable

from robocrys.stats.adapter import StatisticsAdapter


class StructureStatistics(object):
    """Class to calculate stats of structure features in a set of structures.

    Args:
        condensed_structures: A collection of condensed structures. Supports
            abstract databases of condensed structures provided they are
            iterable. For example a pymongo Cursor object.
        inc_oxi_state: Whether to include the oxidation state in the element
            string (if available)).
    """

    def __init__(self, condensed_structures: Iterable[Dict],
                 inc_oxi_state: bool = False):
        self.stat_adapters = [StatisticsAdapter(cs, inc_oxi_state=inc_oxi_state)
                              for cs in condensed_structures]

        self.element_counts = Counter()
        for sa in self.stat_adapters:
            self.element_counts += Counter(sa.species)

    def get_element_connectivity_probabilities(self):
        connectivity_counts = defaultdict(Counter)

        for sa in self.stat_adapters:
            for element, counts in sa.get_element_connectivity_count(
                    normalize=True).items():
                connectivity_counts[element] += Counter(counts)

        return {el: _divide_counter(counts, self.element_counts[el])
                for el, counts in connectivity_counts.items()}

    def get_element_dimensionality_probabilities(self):
        dimensionality_counts = defaultdict(Counter)

        for sa in self.stat_adapters:
            counts = Counter([sa.dimensionality])
            for element in sa.species:
                dimensionality_counts[element] += counts

        return {el: _divide_counter(counts, self.element_counts[el])
                for el, counts in dimensionality_counts.items()}

    def get_element_component_dimensionality_probabilities(self):
        component_counts = Counter()
        dimensionality_counts = defaultdict(Counter)

        for sa in self.stat_adapters:
            for component in sa.component_makeup:
                counts = Counter([sa.components[component]['dimensionality']])
                for element in sa.component_species[component]:
                    dimensionality_counts[element] += counts

                component_counts += Counter(sa.component_species[component])

        return {el: _divide_counter(counts, component_counts[el])
                for el, counts in dimensionality_counts.items()}

    def get_element_corner_sharing_octahedral_tilt_angles(self):
        angles = defaultdict(list)

        for sa in self.stat_adapters:
            el_angles = sa.get_element_corner_sharing_octahedra_tilt_angles()
            for element, angle in el_angles.items():
                angles[element].append(angle)

        return {el: mean(angles) for el, angles in angles.items()}


def _divide_counter(counter, denominator):
    return {k: v/denominator for k, v in counter.items()}
