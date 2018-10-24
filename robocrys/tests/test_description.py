"""Currently this is just a backup and should be ignored"""

import unittest

from robocrys import SiteDescriber

from pymatgen.core.structure import Structure


class TestDescriptionMethods(unittest.TestCase):
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
