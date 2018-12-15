from robocrys.condense.mineral import MineralMatcher
from robocrys.util import RobocrysTest


class TestMineralMatcher(RobocrysTest):
    """Class to test mineral matching functionality."""

    def setUp(self):
        self.matcher = MineralMatcher()
        self.tin_dioxide = self.get_structure("SnO2")
        self.double_perov = self.get_structure("double_perovskite")

    def test_get_aflow_matches(self):
        """Test AFLOW prototype matching."""
        matches = self.matcher.get_aflow_matches(self.tin_dioxide)

        self.assertEqual(len(matches), 1,
                         msg="number of matches is not equal to 1")
        self.assertEqual(matches[0]['type'], 'Rutile',
                         msg="SnO2 mineral name incorrect")
        self.assertAlmostEqual(matches[0]['distance'], 0.15047694852244528,
                               msg="SnO2 fingerprint distance does not match")
        self.assertTrue('structure' in matches[0] and matches[0]['structure'],
                        msg="SnO2 structure not present in match dictionary")

    def test_get_fingerprint_matches(self):
        """Test fingerprint based matching."""
        matches = self.matcher.get_fingerprint_matches(self.tin_dioxide)
        self.assertEqual(len(matches), 4,
                         msg="number of matches is not equal to 1")
        self.assertEqual(matches[0]['type'], 'Hydrophilite',
                         msg="SnO2 mineral name incorrect")
        self.assertAlmostEqual(matches[0]['distance'], 0.1429748846147379,
                               msg="SnO2 fingerprint distance does not match")
        self.assertTrue('structure' in matches[0] and matches[0]['structure'],
                        msg="SnO2 structure not present in match dictionary")

        # test fingerprint only matches same number of atoms
        matches = self.matcher.get_fingerprint_matches(self.double_perov)
        self.assertEqual(matches, None)

        # test fingerprint can match different number of atoms
        matches = self.matcher.get_fingerprint_matches(self.double_perov,
                                                       match_n_sp=False)
        self.assertEqual(len(matches), 1,
                         msg="double perovskite number of matches not correct")
        self.assertEqual(matches[0]['type'], '(Cubic) Perovskite',
                         msg="Double perovskite mineral name incorrect")
        self.assertAlmostEqual(
            matches[0]['distance'], 0.11697185,
            msg="double perovskite fingerprint distance does not match")
        self.assertTrue(
            'structure' in matches[0] and matches[0]['structure'],
            msg="perovskite structure not present in match dictionary")

    def test_get_best_mineral_name(self):
        """Test mineral name matching."""
        mineral_data = self.matcher.get_best_mineral_name(self.tin_dioxide)
        self.assertEqual(mineral_data['type'], 'Rutile')
        self.assertEqual(mineral_data['distance'], -1.)
        self.assertEqual(mineral_data['n_species_type_match'], True)

        mineral_data = self.matcher.get_best_mineral_name(self.double_perov)
        self.assertEqual(mineral_data['type'], '(Cubic) Perovskite')
        self.assertAlmostEqual(mineral_data['distance'], 0.116971854532)
        self.assertEqual(mineral_data['n_species_type_match'], False)
