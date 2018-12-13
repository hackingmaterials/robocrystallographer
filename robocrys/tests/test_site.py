from pymatgen.analysis.local_env import CrystalNN

from robocrys.site import SiteAnalyzer, geometries_match, nn_summaries_match, \
    nnn_summaries_match
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
        self.assertEqual(info[0]["inequiv_id"], 4)
        self.assertAlmostEqual(info[0]["dist"], 2.0857160137)

        # check different structure
        sa = SiteAnalyzer(self.ba_n)
        info = sa._get_nearest_neighbor_info(0)
        self.assertEqual(len(info), 6)
        self.assertEqual(info[0]["element"], "N")
        self.assertEqual(info[0]["inequiv_id"], 0)
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
        info = sa.get_next_nearest_neighbor_data(5)
        self.assertTrue("octahedral" in info["Sn"])
        self.assertTrue('corner-sharing' in info["Sn"]['octahedral'])
        self.assertEqual(info["Sn"]['octahedral']['corner-sharing']['n_sites'],
                         8)
        self.assertAlmostEqual(
            info["Sn"]["octahedral"]['corner-sharing']["angles"][0],
            130.16984393647132)
        self.assertEqual(
            len(info["Sn"]["octahedral"]['corner-sharing']["angles"]),
            8)

    def test_get_nearest_neighbor_summary(self):
        """Check getting nearest neighbor summary for all neighbors."""
        sa = SiteAnalyzer(self.tin_dioxide)
        info = sa.get_nearest_neighbor_data(5)
        self.assertTrue("O" in info)
        self.assertEqual(info["O"]["n_sites"], 6)
        self.assertEqual(info["O"]['inequiv_groups'][0]['n_sites'], 1)
        self.assertEqual(info["O"]['inequiv_groups'][0]['inequiv_id'], 2)
        self.assertAlmostEqual(info["O"]['inequiv_groups'][0]['dists'][0],
                               2.051317869078315)

    def test_equivalent_sites(self):
        """Check equivalent sites instance variable set correctly."""
        sa = SiteAnalyzer(self.ba_n, use_symmetry=False)
        self.assertEqual(sa.equivalent_sites, [0, 0, 0, 0, 4, 4])

        # test using symmetry to determine equivalent sites
        sa = SiteAnalyzer(self.ba_n, use_symmetry=True)
        self.assertEqual(sa.equivalent_sites.tolist(), [0, 0, 0, 0, 4, 4])

        # test symprec option works
        sa = SiteAnalyzer(self.ba_n, use_symmetry=True, symprec=0.0001)
        self.assertEqual(sa.equivalent_sites.tolist(), [0, 1, 1, 0, 4, 4])

    def test_get_inequivalent_site_ids(self):
        sa = SiteAnalyzer(self.ba_n, use_symmetry=False)
        inequiv_ids = sa.get_inequivalent_site_ids(list(range(6)))
        self.assertEqual(inequiv_ids, [0, 4])

        # test using symmetry to determine inequivalent sites.
        sa = SiteAnalyzer(self.ba_n, use_symmetry=True)
        inequiv_ids = sa.get_inequivalent_site_ids(list(range(6)))
        self.assertEqual(inequiv_ids, [0, 4])

        # test symprec option
        sa = SiteAnalyzer(self.ba_n, use_symmetry=True, symprec=0.0001)
        inequiv_ids = sa.get_inequivalent_site_ids(list(range(6)))
        self.assertEqual(inequiv_ids, [0, 1, 4])

    def test_geometries_match(self):
        """Test geometry matching function."""
        sa = SiteAnalyzer(self.ba_n, use_symmetry=True)
        geom_a = sa.get_site_geometry(0)
        geom_b = sa.get_site_geometry(1)
        geom_c = sa.get_site_geometry(4)

        self.assertTrue(geometries_match(geom_a, geom_b))
        self.assertFalse(geometries_match(geom_a, geom_c))

        # test likeness tol works
        self.assertFalse(geometries_match(geom_a, geom_b, likeness_tol=1e-10))

    def test_nn_summaries_match(self):
        """Test nearest neighbour summary matching function."""
        sa = SiteAnalyzer(self.ba_n, use_symmetry=True)
        nn_a = sa.get_nearest_neighbor_data(0, split_into_groups=False)
        nn_b = sa.get_nearest_neighbor_data(1, split_into_groups=False)
        nn_c = sa.get_nearest_neighbor_data(4, split_into_groups=False)

        self.assertTrue(nn_summaries_match(nn_a, nn_b))
        self.assertFalse(nn_summaries_match(nn_a, nn_c))

        # test bond dist tol works
        self.assertFalse(nn_summaries_match(nn_a, nn_b, bond_dist_tol=1e-10))

        # test not matching bond dists
        self.assertTrue(nn_summaries_match(nn_a, nn_b, bond_dist_tol=1e-10,
                                           match_bond_dists=False))

    def test_nnn_summaries_match(self):
        """Test nearest neighbour summary matching function."""
        sa = SiteAnalyzer(self.ba_n, use_symmetry=True)
        nnn_a = sa.get_next_nearest_neighbor_data(0)
        nnn_b = sa.get_next_nearest_neighbor_data(1)
        nnn_c = sa.get_next_nearest_neighbor_data(4)

        self.assertTrue(nnn_summaries_match(nnn_a, nnn_b))
        self.assertFalse(nnn_summaries_match(nnn_a, nnn_c))

        # test bond angle tol works
        self.assertFalse(nnn_summaries_match(
            nnn_a, nnn_b, bond_angle_tol=1e-10))

        # test not matching bond angles
        self.assertTrue(nnn_summaries_match(
            nnn_a, nnn_b, bond_angle_tol=1e-10, match_bond_angles=False))
