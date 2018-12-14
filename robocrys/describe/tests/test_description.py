"""Currently this is just a backup and should be ignored"""

from robocrys.condense.site import SiteDescriber
from robocrys.util import RobocrysTest


class TestDescriptionMethods(RobocrysTest):
    """Class to test mineral matching functionality."""

    def setUp(self):
        self.tin_dioxide = self.get_structure("tin_dioxide")

        self.ba_n = self.get_structure("BaN2")

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
        desc = describer._get_site_description(0)
        self.assertNotEqual(desc, None)

        desc = describer._get_site_description(4)
        self.assertNotEqual(desc, None)

        # check output changes when turning of bond length information
        desc_bond_lengths = describer._get_site_description(
            4, describe_bond_lengths=False)
        self.assertNotEqual(desc, desc_bond_lengths)

        # check different structure
        describer = SiteDescriber(self.ba_n)
        desc = describer._get_site_description(0)
        self.assertNotEqual(desc, None)
