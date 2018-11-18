import unittest

from pymatgen.analysis.dimensionality import get_structure_component_info
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from robocrys import MineralMatcher

from pymatgen.core.structure import Structure

from pprint import pprint

from robocrys.component import get_sym_inequiv_components
from robocrys.condenser import StructureCondenser


class TestComponent(unittest.TestCase):
    """Class to test component related functions."""

    def setUp(self):
        self.tin_dioxide = Structure(
            [3.24, 0, 0, 0, 4.83, 0, 0, 0, 4.83],
            ['O', 'O', 'O', 'O', 'Sn', 'Sn'],
            [[0.5, 0.19, 0.80], [0.5, 0.80, 0.19], [0, 0.30, 0.30],
             [0, 0.69, 0.69], [0.5, 0.50, 0.50], [0, 0, 0]]
        )

        self.double_perov = Structure(
            [-5.75, -5.75, -0.0, -5.75, -0.0, -5.75, 0.0, -5.75, -5.75],
            ['Cs', 'Cs', 'Ag', 'Bi', 'Br', 'Br', 'Br', 'Br', 'Br', 'Br'],
            [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5],
             [-0.0, 0.0, -0.0], [0.75, 0.25, 0.25], [0.75, 0.25, 0.75],
             [0.75, 0.75, 0.25], [0.25, 0.75, 0.75], [0.25, 0.75, 0.25],
             [0.25, 0.25, 0.75]]
        )

    def test_get_sym_inequiv_components(self):
        """Test getting symmetrically inequivalent structure components."""
        sga = SpacegroupAnalyzer(self.tin_dioxide, symprec=0.1)
        bonded_structure = CrystalNN().get_bonded_structure(self.tin_dioxide)
        components = get_structure_component_info(
            bonded_structure, inc_orientation=True, inc_site_ids=True)

        inequiv_comp = get_sym_inequiv_components(
            sga.get_symmetry_dataset()['equivalent_atoms'], components)

        pprint(inequiv_comp[0]['site_ids'])
