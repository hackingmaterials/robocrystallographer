from pymatgen.analysis.local_env import CrystalNN

from robocrys.site import SiteAnalyzer
from robocrys.util import RobocrysTest


class TestSiteAnalyzer(RobocrysTest):
    """Class to test site analysis functionality."""

    def setUp(self):
        cnn = CrystalNN()
        self.tin_dioxide = cnn.get_bonded_structure(
            self.get_structure("tin_dioxide"))
        self.ba_n = cnn.get_bonded_structure(
            self.get_structure("BaN2"))

    def test_init(self):
        """Test to check SiteDescriber can be initialised"""
        sa = SiteAnalyzer(self.tin_dioxide)
        self.assertNotEqual(sa, None,
                            msg="tin dioxide site analyzer could not be init")

        # check different structure
        sa = SiteAnalyzer(self.ba_n)
        self.assertNotEqual(sa, None,
                            msg="BaN2 site analyzer could not be initialized")

    def test_get_site_geometry(self):
        """Test site geometry description."""
        sa = SiteAnalyzer(self.tin_dioxide)
        geom_data = sa.get_site_geometry(0)
        self.assertEqual(geom_data['type'], "trigonal planar")
        self.assertAlmostEqual(geom_data['likeness'], 0.65087198357)

        geom_data = sa.get_site_geometry(4)
        self.assertEqual(geom_data['type'], "octahedral")
        self.assertAlmostEqual(geom_data['likeness'], 0.92981547798)

        # check different structure
        sa = SiteAnalyzer(self.ba_n)
        geom_data = sa.get_site_geometry(0)
        self.assertEqual(geom_data['type'], "octahedral")
        self.assertAlmostEqual(geom_data['likeness'], 0.323784117039)

    def test_get_nearest_neighbor_info(self):
        """Check getting nearest neighbor information."""
        sa = SiteAnalyzer(self.tin_dioxide)
        info = sa._get_nearest_neighbor_info(0)

        self.assertEqual(len(info), 3)
        self.assertEqual(info[0]["element"], "Sn")
        self.assertEqual(info[0]["sym_id"], 4)
        self.assertAlmostEqual(info[0]["dist"], 2.0857160137)

        # check different structure
        sa = SiteAnalyzer(self.ba_n)
        info = sa._get_nearest_neighbor_info(0)
        self.assertEqual(len(info), 6)
        self.assertEqual(info[0]["element"], "N")
        self.assertEqual(info[0]["sym_id"], 0)
        self.assertAlmostEqual(info[0]["dist"], 1.2619877178483)

    def test_get_next_nearest_neighbor_info(self):
        """Check getting next nearest neighbor information."""
        sa = SiteAnalyzer(self.tin_dioxide)
        info = sa._get_next_nearest_neighbor_info(0)

        self.assertEqual(len(info), 15)
        self.assertEqual(info[0]["element"], "O")
        self.assertEqual(info[0]["connectivity"], "corner-sharing")
        self.assertEqual(info[0]["geometry"]["type"], "trigonal planar")

        info = sa._get_next_nearest_neighbor_info(5)
        self.assertEqual(info[0]["element"], 'Sn')
        self.assertEqual(info[0]["connectivity"], "corner-sharing")
        self.assertEqual(info[0]["geometry"]["type"], "octahedral")
        self.assertEqual(len(info[0]['angles']), 1)
        self.assertAlmostEqual(info[0]['angles'][0], 130.16984393647132)

        # check different structure
        sa = SiteAnalyzer(self.ba_n)
        info = sa._get_next_nearest_neighbor_info(0)
        self.assertEqual(len(info), 50)
        self.assertEqual(info[5]["element"], "N")
        self.assertEqual(info[5]["connectivity"], "edge-sharing")
        self.assertAlmostEqual(info[5]["geometry"]["type"], "octahedral")
        self.assertEqual(len(info[5]['angles']), 2)
        self.assertAlmostEqual(info[5]['angles'][0], 83.91397867959587)

    def test_get_next_nearest_neighbor_summary(self):
        """Check getting next nearest neighbor summary for all neighbors."""
        sa = SiteAnalyzer(self.tin_dioxide)
        info = sa.get_next_nearest_neighbor_summary(5)
        self.assertEqual(info["Sn"]["n_sites"], 12)
        self.assertTrue("corner-sharing" in info["Sn"])
        self.assertEqual(info["Sn"]['corner-sharing']['n_sites'], 8)
        self.assertEqual(info["Sn"]['corner-sharing']['geometries'][0],
                         'octahedral')
        self.assertAlmostEqual(info["Sn"]["corner-sharing"]["angles"][0],
                               130.16984393647132)
        self.assertEqual(len(info["Sn"]["corner-sharing"]["angles"]), 8)

    def test_get_nearest_neighbor_summary(self):
        """Check getting nearest neighbor summary for all neighbors."""
        sa = SiteAnalyzer(self.tin_dioxide)
        info = sa.get_nearest_neighbor_summary(5)
        self.assertTrue("O" in info)
        self.assertEqual(info["O"]["n_sites"], 6)
        self.assertEqual(info["O"]['sym_groups'][0]['n_sites'], 1)
        self.assertEqual(info["O"]['sym_groups'][0]['sym_id'], 2)
        self.assertAlmostEqual(info["O"]['sym_groups'][0]['dists'][0],
                               2.051317869078315)
