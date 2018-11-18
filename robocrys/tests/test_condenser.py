import unittest

from robocrys import MineralMatcher

from pymatgen.core.structure import Structure

from pprint import pprint
from robocrys.condenser import StructureCondenser


class TestStructureCondenser(unittest.TestCase):
    """Class to test mineral matching functionality."""

    def setUp(self):
        self.tin_dioxide = Structure(
            [3.24, 0, 0, 0, 4.83, 0, 0, 0, 4.83],
            ['O', 'O', 'O', 'O', 'Sn', 'Sn'],
            [[0.5, 0.19, 0.80], [0.5, 0.80, 0.19], [0, 0.30, 0.30],
             [0, 0.69, 0.69], [0.5, 0.50, 0.50], [0, 0, 0]]
        )
        self.tin_dioxide.add_oxidation_state_by_guess()

        self.double_perov = Structure(
            [-5.75, -5.75, -0.0, -5.75, -0.0, -5.75, 0.0, -5.75, -5.75],
            ['Cs', 'Cs', 'Ag', 'Bi', 'Br', 'Br', 'Br', 'Br', 'Br', 'Br'],
            [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5],
             [-0.0, 0.0, -0.0], [0.75, 0.25, 0.25], [0.75, 0.25, 0.75],
             [0.75, 0.75, 0.25], [0.25, 0.75, 0.75], [0.25, 0.75, 0.25],
             [0.25, 0.25, 0.75]]
        )

    def test_init(self):
        sc = StructureCondenser()
        self.assertNotEqual(sc, None)

    def test_condense_structure(self):
        """Test structure condensing."""
        sc = StructureCondenser(symprec=0.1)

        data = sc.condense_structure(self.tin_dioxide)
