"""Nominal tests for command line functions. Note this script does not actually
test the main method with argparsing, just the robocrystallographer function.
"""
from __future__ import annotations

from robocrys.cli import robocrystallographer, main
from robocrys.util.tests import RobocrysTest
import os
import pytest

from pymatgen.core import SETTINGS as PMG_SETTINGS
from robocrys.cli import MPRestError

_mp_api_key = os.environ.get("MP_API_KEY") or PMG_SETTINGS.get("PMG_MAPI_KEY")

class TestCommandLineInterface(RobocrysTest):
    """Class to test CLI functionality."""

    def setUp(self):
        self.tin_dioxide = self.get_structure("SnO2")

    def test_robocrystallographer(self):
        # check not passing any arguments
        description = robocrystallographer(self.tin_dioxide)
        assert all(
            attr in description for attr in [
                "Rutile","SnO2","tetragonal","P4_2/mnm","Sn(1)4+","equivalent","corner","edge","Sn(1)-O(1)","2.09"
            ]
        )

        # check passing condense and describe kwargs
        description = robocrystallographer(
            self.tin_dioxide,
            condenser_kwargs={"use_symmetry_equivalent_sites": True},
            describer_kwargs={
                "describe_symmetry_labels": False,
                "bond_length_decimal_places": 4,
            },
        )
        assert all(
            attr in description for attr in ["Sn4+", "Sn-O", "2.0922"]
        )

    @pytest.fixture(autouse=True)
    def capsys(self, capsys):
        self.capsys = capsys
        
    @pytest.mark.skipif((not _mp_api_key) or (MPRestError is None), reason="No MP API key set, or `mp_api` is not installed.")
    def test_robocrystallographer_mp_api(self):
        
        import sys

        with pytest.MonkeyPatch.context() as monke:
            monke.setattr(sys,"argv",[sys.argv[0], "mp-856"])
            main()
    
        stdout, _ = self.capsys.readouterr()
        assert all(
            attr in stdout for attr in [
                "Rutile","SnO","tetragonal","P4","/mnm","Sn(1)","equivalent","corner","edge","Sn(1)-O(1)","2.06"
            ]
        )