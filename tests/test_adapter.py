from __future__ import annotations

from pytest import approx

from robocrys.adapter import BaseAdapter
from robocrys.util.tests import RobocrysTest


class TestDescriptionAdapter(RobocrysTest):
    """Class to test the base adapter functionality."""

    def setUp(self):
        tin_dioxide = self.get_condensed_structure("SnO2")
        self.tin_dioxide_ba = BaseAdapter(tin_dioxide)

        mapi = self.get_condensed_structure("mapi")
        self.mapi_ba = BaseAdapter(mapi)

    def test_attributes(self):
        assert self.mapi_ba.mineral["type"] == "Orthorhombic Perovskite"
        assert self.mapi_ba.mineral["distance"] == -1
        assert self.mapi_ba.formula == "CH3NH3PbI3"
        assert self.mapi_ba.spg_symbol == "Pnma"
        assert self.mapi_ba.crystal_system == "orthorhombic"
        assert self.mapi_ba.dimensionality == 3
        assert self.mapi_ba.sites
        assert self.mapi_ba.distances
        assert self.mapi_ba.angles
        assert self.mapi_ba.components
        assert self.mapi_ba.component_makeup
        assert self.mapi_ba.elements[0] == "H+"

    def test_get_distance_details(self):
        # test get distance using int works
        distances = self.tin_dioxide_ba.get_distance_details(0, 2)
        assert len(distances), 3
        assert distances[0] == approx(2.0922101061490546)

        # test get distance using list works
        distances = self.tin_dioxide_ba.get_distance_details(0, [2])
        assert len(distances), 3
        assert distances[0] == approx(2.0922101061490546)

        # test getting multiple distances
        distances = self.mapi_ba.get_distance_details(44, [0, 8])
        assert len(distances), 4
        assert distances[0] == approx(1.0386222568611572)

    def test_get_angle_details(self):
        # test get angles using int works
        distances = self.tin_dioxide_ba.get_angle_details(0, 0, "corner")
        assert len(distances), 8
        assert distances[0] == approx(129.18849530149342)

        # test get angles using list works
        distances = self.tin_dioxide_ba.get_angle_details(0, [0], "corner")
        assert len(distances), 8
        assert distances[0] == approx(129.18849530149342)
