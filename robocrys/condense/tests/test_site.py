from pymatgen.analysis.local_env import CrystalNN

from robocrys.condense.site import (SiteAnalyzer, geometries_match,
                                    nn_summaries_match, nnn_summaries_match)
from robocrys.util import RobocrysTest


class TestSiteAnalyzer(RobocrysTest):
    """Class to test site analysis functionality."""

    def setUp(self):
        cnn = CrystalNN()
        self.tin_dioxide = cnn.get_bonded_structure(
            self.get_structure("SnO2"))
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

    def test_equivalent_sites(self):
        """Check equivalent sites instance variable set correctly."""
        sa = SiteAnalyzer(self.ba_n, use_symmetry_equivalent_sites=False)
        self.assertEqual(sa.equivalent_sites, [0, 0, 0, 0, 4, 4])

        # test using symmetry to determine equivalent sites
        sa = SiteAnalyzer(self.ba_n, use_symmetry_equivalent_sites=True)
        self.assertIsInstance(sa.equivalent_sites, list)
        self.assertEqual(sa.equivalent_sites, [0, 0, 0, 0, 4, 4])

        # test symprec option works
        sa = SiteAnalyzer(self.ba_n, use_symmetry_equivalent_sites=True,
                          symprec=0.0001)
        self.assertEqual(sa.equivalent_sites, [0, 1, 1, 0, 4, 4])

    def test_symmetry_labels(self):
        """Check equivalent sites instance variable set correctly."""
        sa = SiteAnalyzer(self.ba_n, use_symmetry_equivalent_sites=False)
        self.assertEqual(sa.symmetry_labels, [1, 1, 1, 1, 1, 1])

    def test_get_site_geometry(self):
        """Test site geometry description."""
        sa = SiteAnalyzer(self.tin_dioxide)
        geom_data = sa.get_site_geometry(0)
        self.assertEqual(geom_data['type'], "octahedral")
        self.assertAlmostEqual(geom_data['likeness'], 0.9349776258427136 )

        geom_data = sa.get_site_geometry(4)
        self.assertEqual(geom_data['type'], "trigonal planar")
        self.assertAlmostEqual(geom_data['likeness'], 0.6050243049545359)

        # check different structure

        sa = SiteAnalyzer(self.ba_n)
        geom_data = sa.get_site_geometry(0)
        self.assertEqual(geom_data['type'], "6-coordinate")
        self.assertAlmostEqual(geom_data['likeness'], 1)

    def test_get_nearest_neighbors(self):
        """Check getting nearest neighbors."""
        sa = SiteAnalyzer(self.tin_dioxide)
        info = sa.get_nearest_neighbors(4)

        self.assertEqual(len(info), 3)
        self.assertEqual(info[0]["element"], "Sn4+")
        self.assertEqual(info[0]["inequiv_index"], 0)
        self.assertAlmostEqual(info[0]["dist"], 2.0922101061490546)

        info = sa.get_nearest_neighbors(0, inc_inequivalent_site_index=False)
        self.assertTrue('inequiv_index' not in info[0])

        # check different structure without oxi state
        sa = SiteAnalyzer(self.ba_n)
        info = sa.get_nearest_neighbors(0)
        self.assertEqual(len(info), 6)
        self.assertEqual(info[0]["element"], "N")
        self.assertEqual(info[0]["inequiv_index"], 0)
        self.assertAlmostEqual(info[0]["dist"], 1.2619877178483)

    def test_get_next_nearest_neighbors(self):
        """Check getting next nearest neighbors."""
        sa = SiteAnalyzer(self.tin_dioxide)
        info = sa.get_next_nearest_neighbors(5)

        self.assertEqual(len(info), 14)
        self.assertEqual(info[0]["element"], "O2-")
        self.assertEqual(info[0]["connectivity"], "edge")
        self.assertEqual(info[0]["geometry"]["type"], "trigonal planar")
        self.assertEqual(info[0]["inequiv_index"], 2)

        info = sa.get_next_nearest_neighbors(
            0, inc_inequivalent_site_index=False)
        self.assertTrue('inequiv_index' not in info[0])

        info = sa.get_next_nearest_neighbors(0)
        self.assertEqual(info[0]["element"], 'Sn4+')
        self.assertEqual(info[0]["connectivity"], "edge")
        self.assertEqual(info[0]["geometry"]["type"], "octahedral")
        self.assertEqual(len(info[0]['angles']), 2)
        self.assertAlmostEqual(info[0]['angles'][0], 101.62287790513848 )

        # check different structure without oxi state
        sa = SiteAnalyzer(self.ba_n)
        info = sa.get_next_nearest_neighbors(0)
        self.assertEqual(len(info), 36)
        self.assertEqual(info[5]["element"], "N")
        self.assertEqual(info[5]["connectivity"], "edge")
        self.assertEqual(info[5]["geometry"]["type"], "6-coordinate")
        self.assertEqual(len(info[5]['angles']), 2)
        self.assertAlmostEqual(info[5]['angles'][0], 83.91397867959587)

    def test_get_site_summary(self):
        """Test getting the site summary"""
        sa = SiteAnalyzer(self.tin_dioxide)
        data = sa.get_site_summary(0)
        self.assertEqual(data['element'], 'Sn4+')
        self.assertEqual(data['geometry']['type'], 'octahedral')
        self.assertAlmostEqual(data['geometry']['likeness'],
                               0.9349776258427136)
        self.assertEqual(len(data['nn']), 6)
        self.assertEqual(len(data['nnn']['corner']), 8)
        self.assertEqual(len(data['nnn']['edge']), 2)
        self.assertEqual(data['poly_formula'], 'O6')
        self.assertEqual(data['sym_labels'], (1, ))

    def test_get_bond_distance_summary(self):
        """Test getting the bond distance summary"""
        sa = SiteAnalyzer(self.tin_dioxide)
        data = sa.get_bond_distance_summary(0)

        self.assertEqual(len(data[2]), 6)
        self.assertAlmostEqual(data[2][0], 2.0922101061490546)

    def test_connectivity_angle_summary(self):
        """Test getting the connectivity angle summary"""
        sa = SiteAnalyzer(self.tin_dioxide)
        data = sa.get_connectivity_angle_summary(0)

        self.assertEqual(len(data[0]['corner']), 8)
        self.assertEqual(len(data[0]['edge']), 4)
        self.assertAlmostEqual(data[0]['edge'][0], 101.62287790513848)

    def test_get_all_site_summaries(self):
        """Test getting all the site summaries."""
        sa = SiteAnalyzer(self.tin_dioxide)
        data = sa.get_all_site_summaries()
        self.assertEqual(len(data.keys()), 2)
        self.assertEqual(data[0]['element'], 'Sn4+')
        self.assertEqual(data[0]['geometry']['type'], 'octahedral')
        self.assertAlmostEqual(data[0]['geometry']['likeness'],
                               0.9349776258427136)
        self.assertEqual(len(data[0]['nn']), 6)
        self.assertEqual(len(data[0]['nnn']['corner']), 8)
        self.assertEqual(len(data[0]['nnn']['edge']), 2)
        self.assertEqual(data[0]['poly_formula'], 'O6')
        self.assertEqual(data[0]['sym_labels'], (1, ))

    def test_get_all_bond_distance_summaries(self):
        sa = SiteAnalyzer(self.tin_dioxide)
        data = sa.get_all_bond_distance_summaries()
        self.assertEqual(len(data[0][2]), 6)
        self.assertEqual(len(data[2][0]), 3)
        self.assertAlmostEqual(data[0][2][0], 2.0922101061490546)

    def test_get_all_connectivity_angle_summaries(self):
        sa = SiteAnalyzer(self.tin_dioxide)
        data = sa.get_all_connectivity_angle_summaries()
        self.assertEqual(len(data[0][0]['corner']), 8)
        self.assertEqual(len(data[0][0]['edge']), 4)
        self.assertAlmostEqual(data[0][0]['edge'][0], 101.62287790513848)

    def test_get_inequivalent_site_indices(self):
        sa = SiteAnalyzer(self.ba_n, use_symmetry_equivalent_sites=False)
        inequiv_indices = sa.get_inequivalent_site_indices(list(range(6)))
        self.assertEqual(inequiv_indices, [0, 4])

        # test using symmetry to determine inequivalent sites.
        sa = SiteAnalyzer(self.ba_n, use_symmetry_equivalent_sites=True)
        inequiv_indices = sa.get_inequivalent_site_indices(list(range(6)))
        self.assertEqual(inequiv_indices, [0, 4])

        # test symprec option
        sa = SiteAnalyzer(self.ba_n, use_symmetry_equivalent_sites=True,
                          symprec=0.000001)
        inequiv_indices = sa.get_inequivalent_site_indices(list(range(6)))
        self.assertEqual(inequiv_indices, [0, 1, 4])

    def test_geometries_match(self):
        """Test geometry matching function."""
        sa = SiteAnalyzer(self.tin_dioxide, use_symmetry_equivalent_sites=False)
        geom_a = sa.get_site_geometry(0)
        geom_b = sa.get_site_geometry(1)
        geom_c = sa.get_site_geometry(4)

        self.assertTrue(geometries_match(geom_a, geom_b))
        self.assertFalse(geometries_match(geom_a, geom_c))

    def test_nn_summaries_match(self):
        """Test nearest neighbour summary matching function."""
        sa = SiteAnalyzer(self.ba_n, use_symmetry_equivalent_sites=True)
        nn_a = sa.get_nearest_neighbors(0)
        nn_b = sa.get_nearest_neighbors(1)
        nn_c = sa.get_nearest_neighbors(4)

        self.assertTrue(nn_summaries_match(nn_a, nn_b))
        self.assertFalse(nn_summaries_match(nn_a, nn_c))

        # test bond dist tol works
        self.assertFalse(nn_summaries_match(nn_a, nn_b, bond_dist_tol=1e-10))

        # test not matching bond dists
        self.assertTrue(nn_summaries_match(nn_a, nn_b, bond_dist_tol=1e-10,
                                           match_bond_dists=False))

    def test_nnn_summaries_match(self):
        """Test nearest neighbour summary matching function."""
        sa = SiteAnalyzer(self.ba_n, use_symmetry_equivalent_sites=True)
        nnn_a = sa.get_next_nearest_neighbors(0)
        nnn_b = sa.get_next_nearest_neighbors(3)
        nnn_c = sa.get_next_nearest_neighbors(4)

        self.assertTrue(nnn_summaries_match(nnn_a, nnn_b))
        self.assertFalse(nnn_summaries_match(nnn_a, nnn_c))

        # test not matching bond angles
        self.assertTrue(nnn_summaries_match(
            nnn_a, nnn_b, bond_angle_tol=1e-10, match_bond_angles=False))
