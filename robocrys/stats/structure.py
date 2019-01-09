"""
This module implements tools for calculating the statistics of structure
features across many structures.
"""

from collections import Counter, defaultdict
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

    def get_element_connectivity_probabilities(self):
        element_counts = Counter()
        connectivity_counts = defaultdict(Counter)

        for sa in self.stat_adapters:
            element_counts += Counter(set(sa.elements.values()))
            for element, counts in sa.get_element_connectivity_count(
                    normalize=True).items():
                connectivity_counts[element] += Counter(counts)

        def divide_counter(counter, denominator):
            return {k: v/denominator for k, v in counter.items()}

        return {el: divide_counter(counts, element_counts[el])
                for el, counts in connectivity_counts.items()}
