"""
This module implements a class to resolve the symbolic references in condensed
structure data.
"""

from collections import namedtuple, defaultdict
from typing import Dict, Any, List, Union

from pymatgen.core.periodic_table import get_el_sp

ComponentDetails = namedtuple('ComponentDetails',
                              ['formula', 'count', 'dimensionality',
                               'molecule_name', 'orientation', 'index'])

ComponentGroup = namedtuple('ComponentGroup',
                            ['formula', 'dimensionality', 'count', 'components',
                             'molecule_name'])

SiteGroup = namedtuple('SiteGroup', ['element', 'count', 'sites'])

NeighborSiteDetails = namedtuple('NeighborSiteDetails',
                                 ['element', 'count', 'sites', 'sym_label'])

NextNeighborSiteDetails = namedtuple('NextNeighborSiteDetails',
                                     ['element', 'count', 'geometry', 'sites',
                                      'sym_label', 'connectivity',
                                      'poly_formula'])


class DescriptionAdapter(object):
    """Class to facilitate pulling data from the condensed structure dictionary.

    Attributes:
        elements: The site elements.
        sym_labels: The symmetry labels as strings.
        use_iupac_ordering (bool, optional): Whether to order formulas
            by the iupac "electronegativity" series, defined in
            Table VI of "Nomenclature of Inorganic Chemistry (IUPAC
            Recommendations 2005)". This ordering effectively follows
            the groups and rows of the periodic table, except the
            Lanthanides, Actanides and hydrogen. If set to ``False``, the
            elements will be ordered according to the electronegativity values.

    Args:
        condensed_structure: The condensed structure data, formatted as produced
            by :meth:`robocrys.condense.StructureCondenser.condense_structure`.
    """

    def __init__(self, condensed_structure: Dict[str, Any],
                 use_iupac_ordering: bool = True):
        self._condensed_structure = condensed_structure
        self.use_iupac_ordering = use_iupac_ordering

        self.elements = {site_index: site_data['element']
                         for site_index, site_data in self.sites.items()}
        self.sym_labels = {site_index: self.get_sym_label(site_index)
                           for site_index in self.sites.keys()}

    def get_nearest_neighbor_details(self, site_index: int,
                                     group: bool = False
                                     ) -> List[NeighborSiteDetails]:
        """Gets a summary of all the nearest neighbors to a site.

        Args:
            site_index: An inequivalent site index.
            group: Whether to group all nearest neighbor sites
                with the same element together.

        Returns:
            A :obj:`list` of ``NeighborSiteDetails`` objects, each with the
            attributes:

            - ``element`` (``str``): The element of the nearest neighbor site.
            - ``count`` (``int``): The number of sites of this type.
            - ``sym_label`` (``str``): The symmetry label.
            - ``sites`` (``list[int]``): The site indices representing this
                nearest neighbor. Can be more than one site if
                ``group_by_element=True``.
        """
        nn_sites = self.sites[site_index]['nn']

        nn_dict = defaultdict(list)
        for nn_site in set(nn_sites):
            element = self.sites[nn_site]['element']
            labels = self.sym_labels[nn_site]
            identity = (element,) if group else (element, labels)

            nn_dict[identity].append(
                {'count': nn_sites.count(nn_site),
                 'labels': labels,
                 'site': nn_site})

        nn_details = []
        for identity, nn_group in nn_dict.items():
            sites = [nn_site['site'] for nn_site in nn_group]
            nn_details.append(NeighborSiteDetails(
                element=identity[0],
                sites=sites,
                count=sum([nn_site['count'] for nn_site in nn_group]),
                sym_label=self.get_sym_label(sites)))

        return sorted(nn_details, key=self._site_order)

    def get_next_nearest_neighbor_details(self, site_index: int,
                                          group: bool = False
                                          ) -> List[NextNeighborSiteDetails]:
        """Gets a summary of all the next nearest neighbors to a site.

        We only get the summaries for next nearest neighbor sites that have a
        geometry type listed in :attr:`robocrys.util.connected_geometries` and
        have a ``poly_formula``.

        Args:
            site_index: An inequivalent site index.
            group: Whether to group together all next nearest neighbor sites
                with the same element, connectivity and geometry but different
                symmetry labels.

        Returns:
            A :obj:`list` of ``NextNeighborSiteDetails`` objects, each with the
            attributes:

            - ``element`` (``str``): The element of the next nearest neighbor
                site.
            - ``connectivity`` (``str``): The connectivity type to this site.
            - ``geometry`` (``str``): The geometry type of the next nearest
                neighbor.
            - ``count`` (``int``): The number of sites of this type.
            - ``sym_label`` (``str``): The symmetry label.
            - ``sites`` (``list[int]``): The site indices representing this
                next nearest neighbor. Can be more than one site if
                ``group=True``.
            - ``poly_formula`` (``str``): The polyhedral formula.
        """
        nnn = self.sites[site_index]['nnn']

        # get a list of tuples of (nnn_site_index, connectivity)
        con_data = [(nnn_site_index, connectivity) for connectivity, sites in
                    nnn.items() for nnn_site_index in set(sites)]

        nnn_dict = defaultdict(list)
        for nnn_site, connectivity in con_data:
            poly_formula = self.sites[nnn_site]['poly_formula']
            if not poly_formula:
                # only interested in describing the connectivity to other
                # polyhedral sites of interest.
                continue

            element = self.sites[nnn_site]['element']
            labels = self.sym_labels[nnn_site]
            geometry = self.sites[nnn_site]['geometry']['type']

            if group:
                identity = (element, connectivity, geometry)
            else:
                identity = (element, connectivity, geometry, labels)

            nnn_dict[identity].append(
                {'count': nnn[connectivity].count(nnn_site),
                 'labels': labels,
                 'site': nnn_site,
                 'poly_formula': poly_formula})

        nnn_details = []
        for identity, nnn_group in nnn_dict.items():
            sites = [nnn_site['site'] for nnn_site in nnn_group]
            nnn_details.append(NextNeighborSiteDetails(
                element=identity[0],
                connectivity=identity[1],
                geometry=identity[2],
                sites=sites,
                poly_formula=nnn_group[0]['poly_formula'],
                count=sum([nn_site['count'] for nn_site in nnn_group]),
                sym_label=self.get_sym_label(sites)))

        return sorted(nnn_details, key=self._site_order)

    def get_distance_details(self, from_site: int,
                             to_sites: Union[int, List[int]]) -> List[float]:
        """Gets the bond lengths between two sets of sites.

        Args:
            from_site: An inequivalent site index.
            to_sites: One ore more inequivalent site indices.

        Returns:
            The distances between the sites.
        """
        if isinstance(to_sites, int):
            # If only one to_site is provided turn it into a list
            to_sites = [to_sites]

        return [distance for to_site in to_sites
                for distance in self.distances[from_site][to_site]]

    def get_angle_details(self, from_site: int, to_sites: Union[int, List[int]],
                          connectivity: str) -> List[float]:
        """Gets the connectivity angles between two sets of sites.

        Args:
            from_site: An inequivalent site index.
            to_sites: One ore more inequivalent site indices.
            connectivity: The site connectivity type. I.e. "corner", "edge", or
                "face".

        Returns:
            The distances between the sites.
        """
        if isinstance(to_sites, int):
            # If only one to_site is provided turn it into a list
            to_sites = [to_sites]

        return [angle for to_site in to_sites
                for angle in self.angles[from_site][to_site][connectivity]]

    def get_component_details(self) -> List[ComponentDetails]:
        """Gets a summary of all components.

        Returns:
            A :obj:`list` of ``ComponentDetails`` objects, each with the
            attributes:

            - ``count`` (``int``): The number of these components in the
                structure.
            - ``formula`` (``str``): The component formula.
            - ``dimensionality`` (``int``): The component dimensionality.
            - ``molecule_name`` (``str`` or ``None``): The molecule name if
                applicable, else ``None``.
            - ``orientation`` (``tuple[int]``): The component orientation.
            - ``index`` (``list[int]``): The component inequivalent index.
        """
        component_details = []

        for index in set(self.component_makeup):
            component_details.append(ComponentDetails(
                count=self.component_makeup.count(index),
                formula=self.components[index]['formula'],
                dimensionality=self.components[index]['dimensionality'],
                molecule_name=self.components[index]['molecule_name'],
                orientation=self.components[index]['orientation'],
                index=index))

        return sorted(component_details, key=_component_order)

    def get_component_groups(self) -> List[ComponentGroup]:
        """Gets a summary of all components groups.

        Returns:
            The components, grouped together by formula, dimensionality and
            molecule name. The data will be returned as a :obj:`list` of
            ``ComponentGroup`` objects, each with the attributes:

            - ``count`` (``int``): The total number of components in this group.
            - ``formula`` (``str``): The formula of the components..
            - ``dimensionality`` (``int``): The dimensionality of the
                components.
            - ``molecule_name`` (``str`` or ``None``): The molecule name if
                applicable, else ``None``.
            - ``components`` (``list[ComponentDetails]``): The components
                in the group.
        """
        component_details = self.get_component_details()

        grouped_components = defaultdict(list)
        for component in component_details:
            identity = (component.dimensionality, component.formula,
                        component.molecule_name)
            grouped_components[identity].append(component)

        component_group_details = []
        for identity, group in grouped_components.items():
            component_group_details.append(ComponentGroup(
                count=sum(component.count for component in group),
                dimensionality=identity[0],
                formula=identity[1],
                molecule_name=identity[2],
                components=sorted(group, key=_component_order)))

        return sorted(component_group_details, key=_component_order)

    def get_component_site_groups(self, component_index: int
                                  ) -> List[SiteGroup]:
        """Gets a summary of the sites in a component.

        Returns:
            The sites, grouped together by element. The data will be returned
            as a :obj:`list` of ``SiteGroup`` objects, each with the attributes:

            - ``count`` (``int``): The total number of sites in this group.
            - ``element`` (``str``): The site element.
            - ``sites`` (``list[int]``): A list of site indices in this group.
        """
        sites = self.components[component_index]['sites']

        grouped_sites = defaultdict(list)
        for site_index in sites:
            grouped_sites[self.elements[site_index]].append(site_index)

        site_groups = []
        for element, group in grouped_sites.items():
            site_groups.append(SiteGroup(
                count=sum(sites.count(site_index) for site_index in group),
                element=element,
                sites=group))

        return sorted(site_groups, key=self._site_order)

    def get_sym_label(self, site_indices: Union[int, List[int]]) -> str:
        """Convert site indices into a formatted symmetry label.

        Args:
            site_indices:  THe site indices.

        Returns:
            The formatted symmetry label. E.g., if the set of symmetry labels
            for the sites looks like ``(1, 2)``, the symmetry label will be
            ``(1,2)``.
        """
        if isinstance(site_indices, int):
            site_indices = [site_indices]

        all_labels = sorted([label for site_index in site_indices for label in
                             self.sites[site_index]['sym_labels']])
        return "({})".format(",".join(map(str, sorted(all_labels))))

    @property
    def mineral(self) -> Dict[str, Union[str, int, bool]]:
        """The mineral data.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure['mineral']

    @property
    def formula(self) -> str:
        """The structure formula.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure['formula']

    @property
    def spg_symbol(self) -> str:
        """The space group symbol.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure['spg_symbol']

    @property
    def crystal_system(self) -> str:
        """The crystal system.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure['crystal_system']

    @property
    def dimensionality(self) -> int:
        """The overall dimensionality.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure['dimensionality']

    @property
    def sites(self) -> Dict[int, Dict[str, Any]]:
        """The site data.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure['sites']

    @property
    def distances(self) -> Dict[int, Dict[int, List[float]]]:
        """The distance data.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure['distances']

    @property
    def angles(self) -> Dict[int, Dict[int, Dict[str, List[float]]]]:
        """The angle data.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure['angles']

    @property
    def components(self) -> Dict[int, Dict[str, Any]]:
        """The component data.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure['components']

    @property
    def component_makeup(self) -> List[int]:
        """The component makeup of the structure.

        See :meth:`robocrys.condense.StructureCondenser.condense_structure` for
        more details.
        """
        return self._condensed_structure['component_makeup']

    def _site_order(self, s: Union[SiteGroup,
                                   NeighborSiteDetails,
                                   NextNeighborSiteDetails]):
        """Utility function to help sort NeighborSiteDetails and SiteGroups."""
        specie = get_el_sp(s.element)
        x = specie.iupac_ordering if self.use_iupac_ordering else specie.X

        if isinstance(s, NeighborSiteDetails):
            return [x, s.count, s.sym_label, s.sites]
        elif isinstance(s, NextNeighborSiteDetails):
            return [s.connectivity, s.geometry, s.count, x, s.poly_formula,
                    s.sym_label, s.sites]
        else:
            return [x, s.count, s.sites]


def _component_order(c: Union[ComponentDetails, ComponentGroup]):
    """Utility function to help sort ComponentDetails and ComponentGroups."""
    mn = c.molecule_name if c.molecule_name else 'z'

    if isinstance(c, ComponentDetails):
        ori = c.orientation if c.orientation else (0, 0, 0)
        return [mn, c.dimensionality, c.formula, ori, c.count]
    else:
        return [mn, c.dimensionality, c.formula, c.count]
