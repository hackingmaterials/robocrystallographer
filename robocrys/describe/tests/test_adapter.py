from robocrys.describe.adapter import DescriptionAdapter
from robocrys.util import RobocrysTest


class TestDescriptionAdapter(RobocrysTest):
    """Class to test the description adapter functionality."""

    def setUp(self):
        tin_dioxide = self.get_condensed_structure("SnO2")
        self.tin_dioxide_da = DescriptionAdapter(tin_dioxide)

        mapi = self.get_condensed_structure("mapi")
        self.mapi_da = DescriptionAdapter(mapi)

    def test_attributes(self):
        self.assertEqual(self.mapi_da.mineral['type'],
                         "Orthorhombic Perovskite")
        self.assertEqual(self.mapi_da.mineral['distance'], -1)
        self.assertEqual(self.mapi_da.formula, "CH3NH3PbI3")
        self.assertEqual(self.mapi_da.spg_symbol, "Pnma")
        self.assertEqual(self.mapi_da.crystal_system, "orthorhombic")
        self.assertEqual(self.mapi_da.dimensionality, 3)
        self.assertTrue(self.mapi_da.sites)
        self.assertTrue(self.mapi_da.distances)
        self.assertTrue(self.mapi_da.angles)
        self.assertTrue(self.mapi_da.components)
        self.assertTrue(self.mapi_da.component_makeup)
        self.assertEqual(self.mapi_da.elements[0], 'H+')
        self.assertEqual(self.mapi_da.sym_labels[0], '(1)')

    def test_get_nearest_neighbor_details(self):
        """Check getting nearest neighbor summary for all neighbors."""
        all_details = self.tin_dioxide_da.get_nearest_neighbor_details(0)

        self.assertEqual(len(all_details), 1)
        details = all_details[0]
        self.assertEqual(details.element, "O2-")
        self.assertEqual(details.count, 6)
        self.assertEqual(details.sites, [2])
        self.assertEqual(details.sym_label, '(1)')

        # test merging of sym labels
        all_details = self.mapi_da.get_nearest_neighbor_details(
            24, group=True)
        details = all_details[0]
        self.assertEqual(details.element, "I-")
        self.assertEqual(details.count, 6)
        self.assertEqual(details.sites, [32, 36])
        self.assertEqual(details.sym_label, '(1,2)')

    def test_get_next_nearest_neighbor_details(self):
        all_details = self.tin_dioxide_da.get_next_nearest_neighbor_details(0)

        self.assertEqual(len(all_details), 2)
        details = all_details[0]
        self.assertEqual(details.element, "Sn4+")
        self.assertEqual(details.count, 8)
        self.assertEqual(details.sites, [0])
        self.assertEqual(details.geometry, 'octahedral')
        self.assertEqual(details.connectivity, 'corner')
        self.assertEqual(details.poly_formula, 'O6')
        self.assertEqual(details.sym_label, '(1)')

    def test_get_distance_details(self):
        # test get distance using int works
        distances = self.tin_dioxide_da.get_distance_details(0, 2)
        self.assertTrue(len(distances), 3)
        self.assertAlmostEqual(distances[0], 2.0922101061490546)

        # test get distance using list works
        distances = self.tin_dioxide_da.get_distance_details(0, [2])
        self.assertTrue(len(distances), 3)
        self.assertAlmostEqual(distances[0], 2.0922101061490546)

        # test getting multiple distances
        distances = self.mapi_da.get_distance_details(44, [0, 8])
        self.assertTrue(len(distances), 4)
        self.assertAlmostEqual(distances[0], 1.0386222568611572)

    def test_get_angle_details(self):
        # test get angles using int works
        distances = self.tin_dioxide_da.get_angle_details(0, 0, 'corner')
        self.assertTrue(len(distances), 8)
        self.assertAlmostEqual(distances[0], 129.18849530149342)

        # test get angles using list works
        distances = self.tin_dioxide_da.get_angle_details(0, [0], 'corner')
        self.assertTrue(len(distances), 8)
        self.assertAlmostEqual(distances[0], 129.18849530149342)

    def test_get_component_details(self):
        """Check getting component details"""
        all_details = self.mapi_da.get_component_details()
        self.assertEqual(len(all_details), 2)

        details = all_details[0]
        self.assertEqual(details.formula, "CH3NH3")
        self.assertEqual(details.count, 4)
        self.assertEqual(details.dimensionality, 0)
        self.assertEqual(details.molecule_name, "methylammonium")
        self.assertEqual(details.orientation, None)
        self.assertEqual(details.index, 0)

    def test_get_component_summary(self):
        """Check getting the component summaries."""
        all_groups = self.mapi_da.get_component_groups()
        self.assertEqual(len(all_groups), 2)

        group = all_groups[0]
        self.assertEqual(group.count, 4)
        self.assertEqual(group.formula, 'CH3NH3')
        self.assertEqual(group.dimensionality, 0)
        self.assertEqual(group.molecule_name, "methylammonium")

        details = group.components[0]
        self.assertEqual(details.formula, "CH3NH3")
        self.assertEqual(details.count, 4)
        self.assertEqual(details.dimensionality, 0)
        self.assertEqual(details.molecule_name, "methylammonium")
        self.assertEqual(details.orientation, None)
        self.assertEqual(details.index, 0)

    def test_component_site_groups(self):
        """Check getting the SiteGroups in a component."""
        all_groups = self.mapi_da.get_component_site_groups(1)
        self.assertEqual(len(all_groups), 2)
        group = all_groups[1]
        self.assertEqual(group.element, 'I-')
        self.assertEqual(group.count, 2)
        self.assertEqual(group.sites, [32, 36])

    def test_get_sym_label(self):
        """Test getting symmetry labels."""
        self.assertEqual(self.mapi_da.get_sym_label(0), '(1)')

        # test using list to get sym label
        self.assertEqual(self.mapi_da.get_sym_label([0]), '(1)')

        # test combining multiple labels
        self.assertEqual(self.mapi_da.get_sym_label([0, 24]), '(1,1)')
