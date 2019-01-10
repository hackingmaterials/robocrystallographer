"""
This module implements a class to resolve the symbolic references in condensed
structure data.
"""
from collections import Counter, defaultdict
from statistics import mean

from typing import Dict, Any

from pymatgen.core.periodic_table import get_el_sp
from robocrys.adapter import BaseAdapter
from robocrys.util import defaultdict_to_dict


class StatisticsAdapter(BaseAdapter):
    """Class to facilitate calculating statistics from condensed structures.

    Attributes:
        species: A list of species present in the structure.

    Args:
        condensed_structure: The condensed structure data, formatted as produced
            by :meth:`robocrys.condense.StructureCondenser.condense_structure`.
        inc_oxi_state: Whether to include the oxidation state in the element
            string (if available)).
    """

    def __init__(self, condensed_structure: Dict[str, Any],
                 inc_oxi_state: bool = False):
        super().__init__(condensed_structure)

        if not inc_oxi_state:
            self.elements = {site_index: get_el_sp(element).name
                             for site_index, element in self.elements.items()}
        self.species = list(set(self.elements.values()))
        self.component_species = {
            c_index: list(set([self.elements[s] for s in c_data['sites']]))
            for c_index, c_data in self.components.items()}

    def get_site_connectivities(self):
        """Get the connectivity types for each site.

        Returns:
            The connectivity types for each site, returned as a dict of::

                {connectivity: bool, etc}

            where bool is whether the site is connected to other polyhedra via
            that connectivity.

        """

        def connectivity_types(nnn_sites):
            connectivities = {'corner': False, 'edge': False, 'face': False}
            for connectivity, sites in nnn_sites.items():
                if any(self.sites[site]['poly_formula'] for site in sites):
                    connectivities[connectivity] = True
            return connectivities

        return {site_index: connectivity_types(site_data['nnn'])
                for site_index, site_data in self.sites.items()}

    def get_element_count(self) -> Dict[str, int]:
        """Get the number and types of elements in the structure.

        Returns:
            The number and types of elements in the structure. Formatted as::

                {element: count}
        """
        all_elements = [self.elements[site]
                        for component_index in self.component_makeup
                        for site in self.components[component_index]['sites']]

        return dict(Counter(all_elements))

    def get_element_connectivity_count(self, normalize: bool = False,
                                       ) -> Dict[str, Dict[str, int]]:
        """Count the number of connectivity types for each element.

        Returns:
            The number of connectivity types for each element in the structure.
            Formatted as::

                {element: {connectivity: count}}
        """
        connectivities = self.get_site_connectivities()
        poly_sites = [site for component_index in self.component_makeup
                      for site in self.components[component_index]['sites']
                      if self.sites[site]['poly_formula']]
        poly_elements = [self.elements[site] for site in poly_sites]
        n_elements = self.get_element_count()

        counts = defaultdict(lambda: defaultdict(int))
        for element, site_index in zip(poly_elements, poly_sites):
            for connectivity, present in connectivities[site_index].items():
                # uses the fact that bool will be cast to int
                count = present / n_elements[element] if normalize else present
                counts[element][connectivity] += count

        return defaultdict_to_dict(counts)

    def get_element_corner_sharing_octahedra_tilt_angles(
            self) -> Dict[str, Dict[str, int]]:
        """Gets the average tilt angles for corner sharing octahedral sites.

        Returns:
            The average tilt angles for each element in the structure, if
            that element is found in a corner sharing octahedral geometry.
            Formatted as::

                {element: tilt_angle}
        """

        oct_sites = [site for component_index in self.component_makeup
                     for site in self.components[component_index]['sites']
                     if self.sites[site]['geometry']['type'] == 'octahedral' and
                     'corner' in self.sites[site]['nnn']]

        el_angles = defaultdict(list)
        for site in oct_sites:
            nnn_sites = [
                nnn_site for nnn_site in self.sites[site]['nnn']['corner']
                if self.sites[nnn_site]['geometry']['type'] == 'octahedral']
            angles = [abs(180 - angle) for nnn_site in nnn_sites for
                      angle in self.angles[site][nnn_site]['corner']]
            if angles:
                el_angles[self.elements[site]].extend(angles)

        return {el: mean(angles) for el, angles in el_angles.items()}

