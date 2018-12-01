import os
import unittest

from pymatgen.analysis.dimensionality import get_structure_components
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.testing import PymatgenTest

from robocrys.component import (get_sym_inequiv_components,
                                filter_molecular_components,
                                get_reconstructed_structure,
                                get_formula_from_components,
                                get_formula_inequiv_components)

from monty.serialization import loadfn

from robocrys.util import RobocrysTest

test_dir = os.path.join(os.path.dirname(__file__))


class TestComponent(RobocrysTest):
    """Class to test component related functions."""

    def setUp(self):
        cnn = CrystalNN()

        self.mapi = cnn.get_bonded_structure(self.get_structure("mapi"))
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

    def test_get_comp_inequiv_components(self):
        """Test getting compositionally inequivalent structure components."""

        # print(self.mapi_components)
        inequiv_comp = get_formula_inequiv_components(
            self.mapi_components, use_iupac_formula=True)

        # print(inequiv_comp)
        self.assertEqual(len(inequiv_comp), 2)
        self.assertEqual(inequiv_comp[0]['count'], 4)
        self.assertEqual(inequiv_comp[0]['formula'], "CNH6")

        self.assertEqual(inequiv_comp[1]['count'], 4)
        self.assertEqual(inequiv_comp[1]['formula'], 'PbI3')

        # test non-iupac formula
        inequiv_comp = get_formula_inequiv_components(
            self.mapi_components, use_iupac_formula=False)
        self.assertEqual(inequiv_comp[0]['count'], 4)
        self.assertEqual(inequiv_comp[0]['formula'], "H6CN")

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

    def test_get_formula_from_components(self):
        formula = get_formula_from_components(self.mapi_components)
        self.assertTrue(formula, "CNH6PbI3")

        # check non-iupac ordering works
        formula = get_formula_from_components(self.mapi_components,
                                              use_iupac_formula=False)
        self.assertEqual(formula, "H6CNPbI3")

        # test multiple groups of different numbers of compositions
        s = CrystalNN().get_bonded_structure(self.get_structure("CuH8CN5Cl3"))
        comps = get_structure_components(s)
        formula = get_formula_from_components(comps)
        self.assertEqual(formula, "(CuCN4HCl)2(NH2)2(H2)3(HCl)4")

        # test putting molecules first
        s = CrystalNN().get_bonded_structure(self.get_structure("ZrCuH8C2NCl6"))
        comps = get_structure_components(s)
        formula = get_formula_from_components(comps, molecules_first=False)
        self.assertEqual(formula, "ZrCuCl6C2NH8")

        formula = get_formula_from_components(comps, molecules_first=True)
        self.assertEqual(formula, "C2NH8ZrCuCl6")

