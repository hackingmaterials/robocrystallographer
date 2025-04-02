from __future__ import annotations

from robocrys import StructureDescriber
from robocrys.util.tests import RobocrysTest


class TestDescriptionMethods(RobocrysTest):
    """Class to test mineral matching functionality."""

    def setUp(self):
        self.tin_dioxide = self.get_condensed_structure("SnO2")
        self.mapi = self.get_condensed_structure("mapi")

    def test_describe(self):
        """Broad tests to check the right information is in the description."""
        # test general
        d = StructureDescriber(
            describe_oxidation_states=True,
            describe_symmetry_labels=True,
            return_parts=False,
            bond_length_decimal_places=2,
            fmt="raw",
        )
        description = d.describe(self.tin_dioxide)
        assert "Rutile" in description
        assert "SnO2" in description
        assert "tetragonal" in description
        assert "P4_2/mnm" in description
        assert "Sn(1)4+" in description
        assert "equivalent" in description
        assert "corner" in description
        assert "edge" in description
        assert "Sn(1)-O(1)" in description
        assert "2.09" in description

        # test different settings
        d = StructureDescriber(
            describe_oxidation_states=False,
            describe_symmetry_labels=True,
            return_parts=False,
            bond_length_decimal_places=4,
            fmt="raw",
        )
        description = d.describe(self.tin_dioxide)
        assert "Sn(1)" in description
        assert "Sn(1)-O(1)" in description
        assert "2.0922" in description

        # test different settings
        d = StructureDescriber(
            describe_oxidation_states=True,
            describe_symmetry_labels=False,
            return_parts=False,
            bond_length_decimal_places=2,
            fmt="latex",
        )
        description = d.describe(self.tin_dioxide)
        assert "Sn^{4+}" in description
        assert "Sn-O" in description

        # test return parts
        d = StructureDescriber(
            describe_oxidation_states=True,
            describe_symmetry_labels=True,
            return_parts=True,
            bond_length_decimal_places=2,
            fmt="raw",
        )
        description = d.describe(self.tin_dioxide)
        assert "Rutile" in description["mineral"]
        assert "SnO2" in description["mineral"]
        assert "tetragonal" in description["mineral"]
        assert "P4_2/mnm" in description["mineral"]

        assert description["component_makeup"] == ""
        assert "Sn(1)4+" in description["components"]
        assert "equivalent" in description["components"]
        assert "corner" in description["components"]
        assert "edge" in description["components"]
        assert "Sn(1)-O(1)" in description["components"]
        assert "2.09" in description["components"]

    def test_grammar_and_punctuation(self):
        """Check common grammatical errors are not present."""
        d = StructureDescriber()
        description = d.describe(self.tin_dioxide)
        assert ".." not in description
        assert "  " not in description
        assert ". ." not in description

        description = d.describe(self.mapi)
        assert ".." not in description
        assert "  " not in description
        assert ". ." not in description
