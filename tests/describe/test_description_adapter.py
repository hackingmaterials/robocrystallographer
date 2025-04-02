from __future__ import annotations

from pytest import approx

from robocrys.describe.adapter import DescriptionAdapter
from robocrys.util.tests import RobocrysTest


class TestDescriptionAdapter(RobocrysTest):
    """Class to test the description adapter functionality."""

    def setUp(self):
        tin_dioxide = self.get_condensed_structure("SnO2")
        self.tin_dioxide_da = DescriptionAdapter(tin_dioxide)

        mapi = self.get_condensed_structure("mapi")
        self.mapi_da = DescriptionAdapter(mapi)

    def test_attributes(self):
        assert self.mapi_da.mineral["type"] == "Orthorhombic Perovskite"
        assert self.mapi_da.mineral["distance"] == -1
        assert self.mapi_da.formula == "CH3NH3PbI3"
        assert self.mapi_da.spg_symbol == "Pnma"
        assert self.mapi_da.crystal_system == "orthorhombic"
        assert self.mapi_da.dimensionality == 3
        assert self.mapi_da.sites
        assert self.mapi_da.distances
        assert self.mapi_da.angles
        assert self.mapi_da.components
        assert self.mapi_da.component_makeup
        assert self.mapi_da.elements[0] == "H+"
        assert self.mapi_da.sym_labels[0] == "(1)"

    def test_get_nearest_neighbor_details(self):
        """Check getting nearest neighbor summary for all neighbors."""
        all_details = self.tin_dioxide_da.get_nearest_neighbor_details(0)

        assert len(all_details) == 1
        details = all_details[0]
        assert details.element == "O2-"
        assert details.count == 6
        assert details.sites == [2]
        assert details.sym_label == "(1)"

        # test merging of sym labels
        all_details = self.mapi_da.get_nearest_neighbor_details(24, group=True)
        details = all_details[0]
        assert details.element == "I-"
        assert details.count == 6
        assert details.sites == [32, 36]
        assert details.sym_label == "(1,2)"

    def test_get_next_nearest_neighbor_details(self):
        all_details = self.tin_dioxide_da.get_next_nearest_neighbor_details(0)

        assert len(all_details) == 2
        details = all_details[0]
        assert details.element == "Sn4+"
        assert details.count == 8
        assert details.sites == [0]
        assert details.geometry == "octahedral"
        assert details.connectivity == "corner"
        assert details.poly_formula == "O6"
        assert details.sym_label == "(1)"

    def test_get_distance_details(self):
        # test get distance using int works
        distances = self.tin_dioxide_da.get_distance_details(0, 2)
        assert len(distances), 3
        assert distances[0] == approx(2.0922101061490546)

        # test get distance using list works
        distances = self.tin_dioxide_da.get_distance_details(0, [2])
        assert len(distances), 3
        assert distances[0] == approx(2.0922101061490546)

        # test getting multiple distances
        distances = self.mapi_da.get_distance_details(44, [0, 8])
        assert len(distances), 4
        assert distances[0] == approx(1.0386222568611572)

    def test_get_angle_details(self):
        # test get angles using int works
        distances = self.tin_dioxide_da.get_angle_details(0, 0, "corner")
        assert len(distances), 8
        assert distances[0] == approx(129.18849530149342)

        # test get angles using list works
        distances = self.tin_dioxide_da.get_angle_details(0, [0], "corner")
        assert len(distances), 8
        assert distances[0] == approx(129.18849530149342)

    def test_get_component_details(self):
        """Check getting component details."""
        all_details = self.mapi_da.get_component_details()
        assert len(all_details) == 2

        details = all_details[0]
        assert details.formula == "CH3NH3"
        assert details.count == 4
        assert details.dimensionality == 0
        assert details.molecule_name == "methylammonium"
        assert details.orientation is None
        assert details.index == 0

    def test_get_component_summary(self):
        """Check getting the component summaries."""
        all_groups = self.mapi_da.get_component_groups()
        assert len(all_groups) == 2

        group = all_groups[0]
        assert group.count == 4
        assert group.formula == "CH3NH3"
        assert group.dimensionality == 0
        assert group.molecule_name == "methylammonium"

        details = group.components[0]
        assert details.formula == "CH3NH3"
        assert details.count == 4
        assert details.dimensionality == 0
        assert details.molecule_name == "methylammonium"
        assert details.orientation is None
        assert details.index == 0

    def test_component_site_groups(self):
        """Check getting the SiteGroups in a component."""
        all_groups = self.mapi_da.get_component_site_groups(1)
        assert len(all_groups) == 2
        group = all_groups[1]
        assert group.element == "I-"
        assert group.count == 2
        assert group.sites == [32, 36]

    def test_get_sym_label(self):
        """Test getting symmetry labels."""
        assert self.mapi_da.get_sym_label(0) == "(1)"

        # test using list to get sym label
        assert self.mapi_da.get_sym_label([0]) == "(1)"

        # test combining multiple labels
        assert self.mapi_da.get_sym_label([0, 24]) == "(1,1)"
