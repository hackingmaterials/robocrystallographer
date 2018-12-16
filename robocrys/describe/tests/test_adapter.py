from robocrys.describe.adapter import DescriptionAdapter
from robocrys.util import RobocrysTest


class TestStructureCondenser(RobocrysTest):
    """Class to test mineral matching functionality."""

    def setUp(self):
        tin_dioxide = self.get_condensed_structure("SnO2")
        self.da = DescriptionAdapter(tin_dioxide)

    def test_get_nearest_neighbor_details(self):
        """Check getting nearest neighbor summary for all neighbors."""
        data = self.da.get_nearest_neighbor_details(0)

        self.assertTrue("O2-" in data)
        self.assertEqual(data["O2-"]["count"], 6)
        self.assertEqual(data["O2-"]['groups'][0]['count'], 6)
        self.assertEqual(data["O2-"]['groups'][0]['sym_label'], '(1)')
        self.assertAlmostEqual(data["O2-"]['groups'][0]['distances'][0],
                               2.0922101061490546)

    def test_get_component_summary(self):
        """Check getting the component summaries."""
        data = self.da.get_component_details()

        self.assertEqual(len(data), 1)
        self.assertEqual(data[0].count, 1)
        self.assertEqual(data[0].formula, 'SnO2')
        self.assertEqual(data[0].dimensionality, 3)
        self.assertEqual(data[0].molecule_name, None)
        self.assertEqual(data[0].orientation, None)



    # def test_get_next_nearest_neighbor_details(self):
    #     """Check getting next nearest neighbor details."""
    #     info = self.da.get_next_nearest_neighbor_details(5)
    #
    #     self.assertTrue("octahedral" in info["Sn"])
    #     self.assertTrue('corner-sharing' in info["Sn"]['octahedral'])
    #     self.assertEqual(info["Sn"]['octahedral']['corner-sharing']['n_sites'],
    #                      8)
    #     self.assertAlmostEqual(
    #         info["Sn"]["octahedral"]['corner-sharing']["angles"][0],
    #         130.16984393647132)
    #     self.assertEqual(
    #         len(info["Sn"]["octahedral"]['corner-sharing']["angles"]),
    #         8)
