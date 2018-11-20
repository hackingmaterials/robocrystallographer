import os
import unittest

from pymatgen.analysis.dimensionality import get_structure_components
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from robocrys.component import (get_sym_inequiv_components,
                                filter_molecular_components,
                                get_reconstructed_structure)

from monty.serialization import loadfn


test_dir = os.path.join(os.path.dirname(__file__))


class TestComponent(unittest.TestCase):
    """Class to test component related functions."""

    def setUp(self):
        mapi = loadfn(os.path.join(test_dir, 'mapi.json.gz'))
        self.mapi = CrystalNN().get_bonded_structure(mapi)
        self.mapi_components = get_structure_components(
            self.mapi, inc_molecule_graph=True, inc_site_ids=True)

    def test_get_sym_inequiv_components(self):
        """Test getting symmetrically inequivalent structure components."""
        sga = SpacegroupAnalyzer(self.mapi.structure, symprec=0.01)

        inequiv_comp = get_sym_inequiv_components(self.mapi_components, sga)

        self.assertEqual(len(inequiv_comp), 2)
        self.assertEqual(inequiv_comp[0]['count'], 4)
        self.assertEqual(inequiv_comp[0]['inequivalent_ids'],
                         (0, 8, 12, 44, 20, 28))

        self.assertEqual(inequiv_comp[1]['count'], 1)
        self.assertEqual(inequiv_comp[1]['inequivalent_ids'],
                         (32, 24, 36))

    def test_filter_molecular_components(self):
        """Test filtering of molecular components."""
        mol_comps, comps = filter_molecular_components(self.mapi_components)
        mol_dimen = list(set([c['dimensionality'] for c in mol_comps]))
        other_dimen = list(set([c['dimensionality'] for c in comps]))

        self.assertEqual(len(mol_comps), 4)
        self.assertEqual(len(mol_dimen), 1)
        self.assertEqual(mol_dimen[0], 0)

        self.assertEqual(len(comps), 1)
        self.assertTrue(0 not in other_dimen)

    def test_get_reconstructed_structure(self):
        structure = get_reconstructed_structure(
            self.mapi_components, simplify_molecules=False)

        # check the reconstructred structure matches the original using
        # pymatgen's structure matcher.
        sm = StructureMatcher(scale=False, primitive_cell=False, ltol=1e-4,
                              stol=1e-4, angle_tol=1e-4)
        self.assertTrue(sm.fit(structure, self.mapi.structure))

        # Test the structure matches when we simplify molecules.
        # To do this we ignore all C, H, and N atoms.
        structure = get_reconstructed_structure(
            self.mapi_components, simplify_molecules=True)

        sm = StructureMatcher(scale=False, primitive_cell=False, ltol=1e-4,
                              stol=1e-4, angle_tol=1e-4,
                              ignored_species=['C', 'H', 'N'])
        self.assertTrue(sm.fit(structure, self.mapi.structure))
