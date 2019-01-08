from robocrys.adapter import BaseAdapter
from robocrys.tests import RobocrysTest


class TestDescriptionAdapter(RobocrysTest):
    """Class to test the base adapter functionality."""

    def setUp(self):
        tin_dioxide = self.get_condensed_structure("SnO2")
        self.tin_dioxide_ba = BaseAdapter(tin_dioxide)

        mapi = self.get_condensed_structure("mapi")
        self.mapi_ba = BaseAdapter(mapi)

    def test_attributes(self):
        self.assertEqual(self.mapi_ba.mineral['type'],
                         "Orthorhombic Perovskite")
        self.assertEqual(self.mapi_ba.mineral['distance'], -1)
        self.assertEqual(self.mapi_ba.formula, "CH3NH3PbI3")
        self.assertEqual(self.mapi_ba.spg_symbol, "Pnma")
        self.assertEqual(self.mapi_ba.crystal_system, "orthorhombic")
        self.assertEqual(self.mapi_ba.dimensionality, 3)
        self.assertTrue(self.mapi_ba.sites)
        self.assertTrue(self.mapi_ba.distances)
        self.assertTrue(self.mapi_ba.angles)
        self.assertTrue(self.mapi_ba.components)
        self.assertTrue(self.mapi_ba.component_makeup)
        self.assertEqual(self.mapi_ba.elements[0], 'H+')

    def test_get_distance_details(self):
        # test get distance using int works
        distances = self.tin_dioxide_ba.get_distance_details(0, 2)
        self.assertTrue(len(distances), 3)
        self.assertAlmostEqual(distances[0], 2.0922101061490546)

        # test get distance using list works
        distances = self.tin_dioxide_ba.get_distance_details(0, [2])
        self.assertTrue(len(distances), 3)
        self.assertAlmostEqual(distances[0], 2.0922101061490546)

        # test getting multiple distances
        distances = self.mapi_ba.get_distance_details(44, [0, 8])
        self.assertTrue(len(distances), 4)
        self.assertAlmostEqual(distances[0], 1.0386222568611572)

    def test_get_angle_details(self):
        # test get angles using int works
        distances = self.tin_dioxide_ba.get_angle_details(0, 0, 'corner')
        self.assertTrue(len(distances), 8)
        self.assertAlmostEqual(distances[0], 129.18849530149342)

        # test get angles using list works
        distances = self.tin_dioxide_ba.get_angle_details(0, [0], 'corner')
        self.assertTrue(len(distances), 8)
        self.assertAlmostEqual(distances[0], 129.18849530149342)
