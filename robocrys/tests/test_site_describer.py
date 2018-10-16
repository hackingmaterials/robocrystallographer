import unittest

from robocrys import SiteDescriber

from pymatgen.core.structure import Structure


class TestMineralMatcher(unittest.TestCase):
    """Class to test mineral matching functionality."""

    def setUp(self):

        tin_dioxide = Structure(
            [3.24, 0, 0, 0, 4.83, 0, 0, 0, 4.84],
            ['O', 'O', 'O', 'O', 'Sn', 'Sn'],
            [[0.5, 0.19, 0.80], [0.5, 0.80, 0.19], [0, 0.30, 0.30],
             [0, 0.69, 0.69], [0.5, 0.50, 0.50], [0, 0, 0]]
        )

        self.describer = SiteDescriber(tin_dioxide)

    def test_get_site_geometry(self):
        """Test site geometry description."""
        geom = self.describer.get_site_geometry(0)
        self.assertEqual(geom, "trigonal planar")

        geom = self.describer.get_site_geometry(4)
        self.assertEqual(geom, "octahedral")

        # test distorted tolerance
        geom = self.describer.get_site_geometry(0, distorted_tol=0.8)
        self.assertEqual(geom, "distorted trigonal planar")

    def test_get_nearest_neighbour_info(self):
        """Check getting nearest neighbour information."""
        info = self.describer.get_nearest_neighbour_info(0)

        self.assertEqual(len(info), 3)
        self.assertEqual(info[0]["element"], "Sn")
        self.assertEqual(info[0]["sym_id"], 0)
        self.assertAlmostEqual(info[0]["dist"], 2.0857160137)

    def test_get_site_description(self):
        """Check getting site description."""
        info = self.describer.get_site_description(0)
        print(info)


if __name__ == '__main__':
    unittest.main()
