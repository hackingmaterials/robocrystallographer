import unittest

from pymatgen.core.structure import Structure

from robocrys.fingerprint import get_fingerprint, get_fingerprint_distance


class TestFingerprint(unittest.TestCase):
    """Class to test fingerprint functions."""

    def setUp(self):
        self.fe = Structure(
            [2.33, 0.0, -0.82, -1.16, 2.01, -0.82, 0.0, 0.0, 2.47], ['Fe'],
            [[0.0, 0.0, 0.0]]
        )

    def test_get_fingerprint(self):
        """Test to check fingerprinting."""
        fingerprint = get_fingerprint(self.fe)
        self.assertAlmostEqual(fingerprint[4], 1.98432036e-03)

        # test stats option
        fingerprint = get_fingerprint(self.fe, stats=('mean',))
        self.assertAlmostEqual(fingerprint[31], 2.51322893e-01)

        # test preset options
        fingerprint = get_fingerprint(
            self.fe, fingerprint_preset='CrystalNNFingerprint_cn')
        self.assertAlmostEqual(fingerprint[2], 1.98432036e-03)

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
        finger_1 = get_fingerprint(self.fe)
        dist = get_fingerprint_distance(self.fe, finger_1)
        self.assertAlmostEqual(dist, 0.)


if __name__ == '__main__':
    unittest.main()