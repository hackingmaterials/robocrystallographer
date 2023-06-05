"""Nominal tests for command line functions. Note this script does not actually
test the main method with argparsing, just the robocrystallographer function.
"""

from robocrys.cli import robocrystallographer
from robocrys.tests import RobocrysTest


class TestCommandLineInterface(RobocrysTest):
    """Class to test CLI functionality."""

    def setUp(self):
        self.tin_dioxide = self.get_structure("SnO2")

    def test_robocrystallographer(self):
        # check not passing any arguments
        description = robocrystallographer(self.tin_dioxide)
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

        # check passing condense and describe kwargs
        description = robocrystallographer(
            self.tin_dioxide,
            condenser_kwargs={"use_symmetry_equivalent_sites": True},
            describer_kwargs={
                "describe_symmetry_labels": False,
                "bond_length_decimal_places": 4,
            },
        )
        assert "Sn4+" in description
        assert "Sn-O" in description
        assert "2.0922" in description
