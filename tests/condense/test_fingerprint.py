from __future__ import annotations

from pytest import approx

from robocrys.condense.fingerprint import (
    get_fingerprint_distance,
    get_site_fingerprints,
    get_structure_fingerprint,
)
from robocrys.util.tests import RobocrysTest


class TestFingerprint(RobocrysTest):
    """Class to test fingerprint functions."""

    def setUp(self):
        self.fe = self.get_structure("iron")

    def test_get_site_fingerprints(self):
        """Test to check site fingerprinting."""
        finger = get_site_fingerprints(self.fe)[0]
        assert finger["body-centered cubic CN_8"] == approx(0.576950507)

        # check as_dict option
        finger = get_site_fingerprints(self.fe, as_dict=False)[0]
        assert finger[30] == approx(0.576950507)

    def test_get_structure_fingerprint(self):
        """Test to check structure fingerprinting."""
        fingerprint = get_structure_fingerprint(self.fe)
        assert fingerprint[4] == approx(1.98432036e-03)

        # test stats option
        fingerprint = get_structure_fingerprint(self.fe, stats=("mean",))
        assert fingerprint[31] == approx(2.51322893e-01)

        # test preset options - re-enable once fixed
        # fingerprint = get_structure_fingerprint(
        #     self.fe, preset="CrystalNNFingerprint_cn"
        # )
        # assert fingerprint[2] == approx(1.98432036e-03)

    def test_get_fingerprint_distance(self):
        """Tests to check getting fingerprint distance."""
        finger_1 = [0, 0, 0, 1]
        finger_2 = [1, 0, 0, 0]
        dist = get_fingerprint_distance(finger_1, finger_2)
        assert dist == approx(1.4142135623730951)

        # test automatic conversion from structure to fingerprint
        dist = get_fingerprint_distance(self.fe, self.fe)
        assert dist == approx(0.0)

        # test one structure one fingerprint
        finger_1 = get_structure_fingerprint(self.fe)
        dist = get_fingerprint_distance(self.fe, finger_1)
        assert dist == approx(0.0)
