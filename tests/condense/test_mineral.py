from __future__ import annotations

from pytest import approx

from robocrys.condense.mineral import MineralMatcher
from robocrys.util.tests import RobocrysTest


class TestMineralMatcher(RobocrysTest):
    """Class to test mineral matching functionality."""

    def setUp(self):
        self.matcher = MineralMatcher()
        self.tin_dioxide = self.get_structure("SnO2")
        self.double_perov = self.get_structure("double_perovskite")

    def test_get_aflow_matches(self):
        """Test AFLOW prototype matching."""
        matches = self.matcher.get_aflow_matches(self.tin_dioxide)

        assert len(matches) == 1, "number of matches is not equal to 1"
        assert matches[0]["type"] == "Rutile", "SnO2 mineral name incorrect"
        assert matches[0]["distance"] == approx(
            0.15047694852244528
        ), "SnO2 fingerprint distance does not match"
        assert matches[0].get(
            "structure"
        ), "SnO2 structure not present in match dictionary"

    def test_get_fingerprint_matches(self):
        """Test fingerprint based matching."""
        matches = self.matcher.get_fingerprint_matches(self.tin_dioxide)
        assert len(matches) == 4, "number of matches is not equal to 1"
        assert matches[0]["type"] == "Hydrophilite", "SnO2 mineral name incorrect"
        assert matches[0]["distance"] == approx(
            0.1429748846147379
        ), "SnO2 fingerprint distance does not match"
        assert matches[0].get(
            "structure"
        ), "SnO2 structure not present in match dictionary"

        # test fingerprint only matches same number of atoms
        matches = self.matcher.get_fingerprint_matches(self.double_perov)
        assert matches is None

        # test fingerprint can match different number of atoms
        matches = self.matcher.get_fingerprint_matches(
            self.double_perov, match_n_sp=False
        )
        assert len(matches) == 1, "double perovskite number of matches not correct"
        assert (
            matches[0]["type"] == "(Cubic) Perovskite"
        ), "Double perovskite mineral name incorrect"
        assert matches[0]["distance"] == approx(
            0.11697185
        ), "double perovskite fingerprint distance does not match"
        assert matches[0].get(
            "structure"
        ), "perovskite structure not present in match dictionary"

    def test_get_best_mineral_name(self):
        """Test mineral name matching."""
        mineral_data = self.matcher.get_best_mineral_name(self.tin_dioxide)
        assert mineral_data["type"] == "Rutile"
        assert mineral_data["distance"] == -1.0
        assert mineral_data["n_species_type_match"] is True

        mineral_data = self.matcher.get_best_mineral_name(self.double_perov)
        assert mineral_data["type"] == "(Cubic) Perovskite"
        assert mineral_data["distance"] == approx(0.116971854532)
        assert mineral_data["n_species_type_match"] is False
