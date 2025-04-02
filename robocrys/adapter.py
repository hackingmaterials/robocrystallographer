"""This module implements a class to resolve the symbolic references in condensed
structure data.
"""

from __future__ import annotations

from typing import Any


class BaseAdapter:
    """Base adapter class for facilitating access to condensed structure data.

    Attributes:
        elements: The site elements.

    Args:
        condensed_structure: The condensed structure data, formatted as produced
            by :meth:`robocrys.condense.StructureCondenser.condense_structure`.
    """

    def __init__(self, condensed_structure: dict[str, Any]):
        self._condensed_structure = condensed_structure
        self.elements = {
            site_index: site_data["element"]
            for site_index, site_data in self.sites.items()
        }

    def get_distance_details(
        self, from_site: int, to_sites: int | list[int]
    ) -> list[float]:
        """Gets the bond lengths between two sets of sites.

        Args:
            from_site: An inequivalent site index.
            to_sites: One or more inequivalent site indices.

        Returns:
            The distances between the sites.
        """
        if isinstance(to_sites, int):
            # If only one to_site is provided turn it into a list
            to_sites = [to_sites]

        return [
            distance
            for to_site in to_sites
            for distance in self.distances[from_site][to_site]
        ]

    def get_angle_details(
        self, from_site: int, to_sites: int | list[int], connectivity: str
    ) -> list[float]:
        """Gets the connectivity angles between two sets of sites.

        Args:
            from_site: An inequivalent site index.
            to_sites: One or more inequivalent site indices.
            connectivity: The site connectivity type. I.e. "corner", "edge", or
                "face".

        Returns:
            The distances between the sites.
        """
        if isinstance(to_sites, int):
            # If only one to_site is provided turn it into a list
            to_sites = [to_sites]

        return [
            angle
            for to_site in to_sites
            for angle in self.angles[from_site][to_site][connectivity]
        ]

    @property
    def mineral(self) -> dict[str, str | int | bool]:
        """The mineral data.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure["mineral"]

    @property
    def formula(self) -> str:
        """The structure formula.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure["formula"]

    @property
    def spg_symbol(self) -> str:
        """The space group symbol.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure["spg_symbol"]

    @property
    def crystal_system(self) -> str:
        """The crystal system.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure["crystal_system"]

    @property
    def dimensionality(self) -> int:
        """The overall dimensionality.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure["dimensionality"]

    @property
    def sites(self) -> dict[int, dict[str, Any]]:
        """The site data.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure["sites"]

    @property
    def distances(self) -> dict[int, dict[int, list[float]]]:
        """The distance data.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure["distances"]

    @property
    def angles(self) -> dict[int, dict[int, dict[str, list[float]]]]:
        """The angle data.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure["angles"]

    @property
    def components(self) -> dict[int, dict[str, Any]]:
        """The component data.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure["components"]

    @property
    def component_makeup(self) -> list[int]:
        """The component makeup of the structure.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure["component_makeup"]

    @property
    def is_vdw_heterostructure(self) -> bool:
        """Whether the structure is a vdw heterostructure.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return bool(self._condensed_structure["vdw_heterostructure_info"])
