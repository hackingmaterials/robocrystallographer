from robocrys import StructureDescriber
from robocrys.util import RobocrysTest


class TestDescriptionMethods(RobocrysTest):
    """Class to test mineral matching functionality."""

    def setUp(self):
        self.tin_dioxide = self.get_condensed_structure("SnO2")
        self.mapi = self.get_condensed_structure("mapi")

    def test_describe(self):
        """Broad tests to check the right information is in the description."""
        # test general
        d = StructureDescriber(describe_oxidation_states=True,
                               describe_symmetry_labels=True,
                               return_parts=False,
                               bond_length_decimal_places=2,
                               latexify=False)
        description = d.describe(self.tin_dioxide)
        self.assertTrue("Rutile" in description)
        self.assertTrue("SnO2" in description)
        self.assertTrue("tetragonal" in description)
        self.assertTrue("P4_2/mnm" in description)
        self.assertTrue("three-dimensional" in description)
        self.assertTrue("Sn(1)4+" in description)
        self.assertTrue("equivalent" in description)
        self.assertTrue("corner" in description)
        self.assertTrue("edge" in description)
        self.assertTrue("Sn(1)–O(1)" in description)
        self.assertTrue("2.09" in description)

        # test different settings
        d = StructureDescriber(describe_oxidation_states=False,
                               describe_symmetry_labels=True,
                               return_parts=False,
                               bond_length_decimal_places=4,
                               latexify=False)
        description = d.describe(self.tin_dioxide)
        self.assertTrue("Sn(1)" in description)
        self.assertTrue("Sn(1)–O(1)" in description)
        self.assertTrue("2.0922" in description)

        # test different settings
        d = StructureDescriber(describe_oxidation_states=True,
                               describe_symmetry_labels=False,
                               return_parts=False,
                               bond_length_decimal_places=2,
                               latexify=True)
        description = d.describe(self.tin_dioxide)
        self.assertTrue(r"Sn^{4+}" in description)
        self.assertTrue("Sn–O" in description)

        # test return parts
        d = StructureDescriber(describe_oxidation_states=True,
                               describe_symmetry_labels=True,
                               return_parts=True,
                               bond_length_decimal_places=2,
                               latexify=False)
        description = d.describe(self.tin_dioxide)
        self.assertTrue("Rutile" in description['mineral'])
        self.assertTrue("SnO2" in description['mineral'])
        self.assertTrue("tetragonal" in description['mineral'])
        self.assertTrue("P4_2/mnm" in description['mineral'])

        self.assertTrue("three-dimensional" in description['component_makeup'])
        self.assertTrue("Sn(1)4+" in description['components'])
        self.assertTrue("equivalent" in description['components'])
        self.assertTrue("corner" in description['components'])
        self.assertTrue("edge" in description['components'])
        self.assertTrue("Sn(1)–O(1)" in description['components'])
        self.assertTrue("2.09" in description['components'])

    def test_grammar_and_punctuation(self):
        """Check common grammatical errors are not present"""
        d = StructureDescriber()
        description = d.describe(self.tin_dioxide)
        self.assertTrue(".." not in description)
        self.assertTrue("  " not in description)
        self.assertTrue(". ." not in description)

        description = d.describe(self.mapi)
        self.assertTrue(".." not in description)
        self.assertTrue("  " not in description)
        self.assertTrue(". ." not in description)
