import unittest

from robocrys import SiteDescriber

from pymatgen.core.structure import Structure


class TestMineralMatcher(unittest.TestCase):
    """Class to test mineral matching functionality."""

    def setUp(self):

        self.tin_dioxide = Structure(
            [3.24, 0, 0, 0, 4.83, 0, 0, 0, 4.84],
            ['O', 'O', 'O', 'O', 'Sn', 'Sn'],
            [[0.5, 0.19, 0.80], [0.5, 0.80, 0.19], [0, 0.30, 0.30],
             [0, 0.69, 0.69], [0.5, 0.50, 0.50], [0, 0, 0]]
        )

        self.ba_n = Structure(
            [4.16, -0.02, 0.94, 1.76, 3.77, 0.94, 0.01, 0.01, 7.34],
            ['N', 'N', 'N', 'N', 'Ba', 'Ba'],
            [[0.15, 0.45, 0.45], [0.55, 0.85, 0.05], [0.45, 0.15, 0.95],
             [0.85, 0.55, 0.55], [0.79, 0.21, 0.25], [0.21, 0.79, 0.75]]
        )

    def test_init(self):
        """Test to check SiteDescriber can be initialised"""
        describer = SiteDescriber(self.tin_dioxide)
        self.assertNotEqual(describer, None,
                            msg="tin dioxide site describer could not be init")

        # check different structure
        describer = SiteDescriber(self.ba_n)
        self.assertNotEqual(describer, None,
                            msg="BaN2 site describer could not be initialized")

    def test_get_site_geometry(self):
        """Test site geometry description."""
        describer = SiteDescriber(self.tin_dioxide)
        geom = describer.get_site_geometry(0)
        self.assertEqual(geom, "trigonal planar")

        geom = describer.get_site_geometry(4)
        self.assertEqual(geom, "octahedral")

        # test distorted tolerance
        geom = describer.get_site_geometry(0, distorted_tol=0.8)
        self.assertEqual(geom, "distorted trigonal planar")

        # check different structure
        describer = SiteDescriber(self.ba_n)
        geom = describer.get_site_geometry(0)
        self.assertEqual(geom, "distorted octahedral")

    def test_get_nearest_neighbour_info(self):
        """Check getting nearest neighbour information."""
        describer = SiteDescriber(self.tin_dioxide)
        info = describer.get_nearest_neighbour_info(0)

        self.assertEqual(len(info), 3)
        self.assertEqual(info[0]["element"], "Sn")
        self.assertEqual(info[0]["sym_id"], 4)
        self.assertAlmostEqual(info[0]["dist"], 2.0857160137)

        # check different structure
        describer = SiteDescriber(self.ba_n)
        info = describer.get_nearest_neighbour_info(0)
        self.assertEqual(len(info), 6)
        self.assertEqual(info[0]["element"], "N")
        self.assertEqual(info[0]["sym_id"], 0)
        self.assertAlmostEqual(info[0]["dist"], 1.2619877178483)

    def test_get_bond_length_description(self):
        """Check getting bond length descriptions"""
        # content liable to change so just check function runs without error
        # and returns *something*.
        describer = SiteDescriber(self.tin_dioxide)
        desc = describer.get_bond_length_description(0, "Sn")
        self.assertNotEqual(desc, None)

        desc = describer.get_bond_length_description(4, "O")
        self.assertNotEqual(desc, None)

        # check different structure
        describer = SiteDescriber(self.ba_n)
        desc = describer.get_bond_length_description(0, "Ba")
        self.assertNotEqual(desc, None)

    def test_get_site_description(self):
        """Check getting site description."""
        # content liable to change so just check function runs without error
        # and returns *something*.
        describer = SiteDescriber(self.tin_dioxide)
        desc = describer.get_site_description(0)
        self.assertNotEqual(desc, None)

        desc = describer.get_site_description(4)
        self.assertNotEqual(desc, None)

        # check output changes when turning of bond length information
        desc_bond_lengths = describer.get_site_description(
            4, describe_bond_lengths=False)
        self.assertNotEqual(desc, desc_bond_lengths)

        # check different structure
        describer = SiteDescriber(self.ba_n)
        desc = describer.get_site_description(0)
        self.assertNotEqual(desc, None)


if __name__ == '__main__':
    unittest.main()
