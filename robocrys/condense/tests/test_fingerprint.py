from robocrys.condense.fingerprint import (get_site_fingerprints,
                                           get_structure_fingerprint,
                                           get_fingerprint_distance)
from robocrys.util import RobocrysTest


class TestFingerprint(RobocrysTest):
    """Class to test fingerprint functions."""

    def setUp(self):
        self.fe = self.get_structure("iron")

    def test_get_site_fingerprints(self):
        """Test to check site fingerprinting."""
        finger = get_site_fingerprints(self.fe)[0]
        self.assertAlmostEqual(finger['body-centered cubic CN_8'], 0.576950507)

        # check as_dict option
        finger = get_site_fingerprints(self.fe, as_dict=False)[0]
        self.assertAlmostEqual(finger[30], 0.576950507)

    def test_get_structure_fingerprint(self):
        """Test to check structure fingerprinting."""
        fingerprint = get_structure_fingerprint(self.fe)
        self.assertAlmostEqual(fingerprint[4], 1.98432036e-03)

        # test stats option
        fingerprint = get_structure_fingerprint(self.fe, stats=('mean',))
        self.assertAlmostEqual(fingerprint[31], 2.51322893e-01)

        # test preset options – reenable once fixed
        # fingerprint = get_structure_fingerprint(
        #     self.fe, preset='CrystalNNFingerprint_cn')
        # self.assertAlmostEqual(fingerprint[2], 1.98432036e-03)

    def test_get_fingerprint_distance(self):
        """Tests to check getting fingerprint distance."""
        finger_1 = [0, 0, 0, 1]
        finger_2 = [1, 0, 0, 0]
        dist = get_fingerprint_distance(finger_1, finger_2)
        self.assertAlmostEqual(dist, 1.4142135623730951)

        # test automatic conversion from structure to fingerprint
        dist = get_fingerprint_distance(self.fe, self.fe)
        self.assertAlmostEqual(dist, 0.)

        # test one structure one fingerprint
        finger_1 = get_structure_fingerprint(self.fe)
        dist = get_fingerprint_distance(self.fe, finger_1)
        self.assertAlmostEqual(dist, 0.)
