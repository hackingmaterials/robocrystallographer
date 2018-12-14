from robocrys.condense.structure import StructureCondenser
from robocrys.util import RobocrysTest


class TestStructureCondenser(RobocrysTest):
    """Class to test mineral matching functionality."""

    def setUp(self):
        self.tin_dioxide = self.get_structure("tin_dioxide")
        self.tin_dioxide.add_oxidation_state_by_guess()

    def test_init(self):
        sc = StructureCondenser()
        self.assertNotEqual(sc, None)

    def test_condense_structure(self):
        """Test structure condensing."""
        sc = StructureCondenser(symprec=0.1)
        data = sc.condense_structure(self.tin_dioxide)

        # Test standard
        self.assertEqual(data['formula'], 'SnO2')
        self.assertEqual(data['spg'], 'P4_2/mnm')
        self.assertEqual(data['mineral']['type'], 'Rutile')
        self.assertEqual(data['dimensionality'], 3)
        self.assertEqual(data['components'][0]['dimensionality'], 3)
        self.assertEqual(data['components'][0]['orientation'], None)
        self.assertEqual(data['components'][0]['count'], 1)
        self.assertEqual(data['components'][0]['formula'], 'SnO2')
        self.assertTrue("sites" in data['components'][0])
        self.assertEqual(len(data['components'][0]['sites']), 2)

        # TODO: add a test for merging of similar sites
