import os

from pymatgen.analysis.dimensionality import get_structure_components
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from robocrys.condense.component import (
    get_sym_inequiv_components, filter_molecular_components,
    get_reconstructed_structure, get_formula_from_components,
    get_formula_inequiv_components, components_are_vdw_heterostructure,
    get_vdw_heterostructure_information, get_component_formula_and_factor,
    get_component_formula, components_are_isomorphic,
    get_structure_inequiv_components)
from robocrys.util import RobocrysTest

test_dir = os.path.join(os.path.dirname(__file__))


class TestComponent(RobocrysTest):
    """Class to test component related functions."""

    def setUp(self):
        cnn = CrystalNN()

        self.mapi = cnn.get_bonded_structure(self.get_structure("mapi"))
        self.mapi_components = get_structure_components(
            self.mapi, inc_molecule_graph=True, inc_site_ids=True,
            inc_orientation=True)

        self.vdw_hetero = cnn.get_bonded_structure(self.get_structure("MoWS4"))
        self.vdw_hetero_components = get_structure_components(
            self.vdw_hetero, inc_molecule_graph=True, inc_site_ids=True,
            inc_orientation=True)

    def test_get_component_formula_and_factor(self):
        """Test getting the component formula and factor."""
        formula, factor = get_component_formula_and_factor(
            self.mapi_components[0], use_common_formulas=True,
            use_iupac_formula=True)
        self.assertEqual(formula, "CH3NH3")
        self.assertEqual(factor, 1)

        formula, factor = get_component_formula_and_factor(
            self.mapi_components[4], use_common_formulas=True,
            use_iupac_formula=True)
        self.assertEqual(formula, "PbI3")
        self.assertEqual(factor, 4)

        # test without common formulas
        formula, factor = get_component_formula_and_factor(
            self.mapi_components[0], use_common_formulas=False,
            use_iupac_formula=True)
        self.assertEqual(formula, "CNH6")
        self.assertEqual(factor, 1)

        # test without common formulas and without iupac formula
        formula, factor = get_component_formula_and_factor(
            self.mapi_components[0], use_common_formulas=False,
            use_iupac_formula=False)
        self.assertEqual(formula, "H6CN")
        self.assertEqual(factor, 1)

    def test_get_component_formula(self):
        """Test getting the component formula."""
        formula = get_component_formula(
            self.mapi_components[0], use_common_formulas=True,
            use_iupac_formula=True)
        self.assertEqual(formula, "CH3NH3")

        formula = get_component_formula(
            self.mapi_components[4], use_common_formulas=True,
            use_iupac_formula=True)
        self.assertEqual(formula, "PbI3")

        # test without common formulas
        formula = get_component_formula(
            self.mapi_components[0], use_common_formulas=False,
            use_iupac_formula=True)
        self.assertEqual(formula, "CNH6")

        # test without common formulas and without iupac formula
        formula = get_component_formula(
            self.mapi_components[0], use_common_formulas=False,
            use_iupac_formula=False)
        self.assertEqual(formula, "H6CN")

    def test_get_sym_inequiv_components(self):
        """Test getting symmetrically inequivalent structure components."""
        sga = SpacegroupAnalyzer(self.mapi.structure, symprec=0.01)

        inequiv_comp = get_sym_inequiv_components(self.mapi_components, sga)

        self.assertEqual(len(inequiv_comp), 2)
        self.assertEqual(inequiv_comp[0]['count'], 4)
        self.assertEqual(inequiv_comp[1]['count'], 1)

    def test_get_comp_inequiv_components(self):
        """Test getting compositionally inequivalent structure components."""

        inequiv_comp = get_formula_inequiv_components(self.mapi_components)

        self.assertEqual(len(inequiv_comp), 2)
        self.assertEqual(inequiv_comp[0]['count'], 4)
        self.assertEqual(inequiv_comp[0]['formula'], "CH3NH3")
        self.assertEqual(inequiv_comp[1]['count'], 4)
        self.assertEqual(inequiv_comp[1]['formula'], 'PbI3')

        # Test not using common_formulas but with iupac formula
        inequiv_comp = get_formula_inequiv_components(
            self.mapi_components, use_iupac_formula=True,
            use_common_formulas=False)
        self.assertEqual(inequiv_comp[0]['count'], 4)
        self.assertEqual(inequiv_comp[0]['formula'], "CNH6")

        # test non-iupac formula
        inequiv_comp = get_formula_inequiv_components(
            self.mapi_components, use_iupac_formula=False,
            use_common_formulas=False)
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
        formula = get_formula_from_components(
            self.mapi_components, use_common_formulas=True,
            use_iupac_formula=True)
        self.assertTrue(formula, "CH3NH3PbI3")

        # check not using common formulas
        formula = get_formula_from_components(
            self.mapi_components, use_common_formulas=False,
            use_iupac_formula=True)
        self.assertTrue(formula, "CNH6PbI3")

        # check non-iupac ordering works
        formula = get_formula_from_components(
            self.mapi_components, use_iupac_formula=False,
            use_common_formulas=False)
        self.assertEqual(formula, "H6CNPbI3")

        # test multiple groups of different numbers of compositions
        s = CrystalNN().get_bonded_structure(self.get_structure("CuH8CN5Cl3"))
        comps = get_structure_components(s)
        formula = get_formula_from_components(comps, use_iupac_formula=True,
                                              use_common_formulas=True)
        self.assertEqual(formula, "(CuCN4HCl)2(NH2)2(H2)3(HCl)4")

        # test putting molecules first
        s = CrystalNN().get_bonded_structure(self.get_structure("ZrCuH8C2NCl6"))
        comps = get_structure_components(s)
        formula = get_formula_from_components(
            comps, molecules_first=False, use_iupac_formula=True,
            use_common_formulas=True)
        self.assertEqual(formula, "ZrCuCl6(CH3)2NH2")

        formula = get_formula_from_components(
            comps, molecules_first=True, use_iupac_formula=True,
            use_common_formulas=True)
        self.assertEqual(formula, "(CH3)2NH2ZrCuCl6")

    def test_components_are_vdw_heterostructure(self):
        result = components_are_vdw_heterostructure(self.vdw_hetero_components)
        self.assertTrue(result)

        result = components_are_vdw_heterostructure(self.mapi_components)
        self.assertFalse(result)

    def test_get_vdw_heterostructure_information(self):
        data = get_vdw_heterostructure_information(self.vdw_hetero_components,
                                                   inc_ordered_components=True,
                                                   inc_intercalants=True)
        self.assertEqual(len(data['ordered_components']), 4)
        self.assertAlmostEqual(
            data['ordered_components'][0][
                'structure_graph'].structure.frac_coords[0][0], 0.33330876)
        self.assertAlmostEqual(
            data['ordered_components'][3][
                'structure_graph'].structure.frac_coords[0][0], 0.6666924)
        self.assertEqual(data['repeating_unit'], ['MoS2', 'WS2'])
        self.assertEqual(data['num_repetitions'], 2)
        self.assertEqual(data['intercalants'], [])

        # test error catching
        self.assertRaises(ValueError, get_vdw_heterostructure_information,
                          self.mapi_components)
        self.assertRaises(KeyError, get_vdw_heterostructure_information,
                          get_structure_components(self.vdw_hetero))

    def test_components_are_isomorphic(self):
        # check two CH3NH3 components are isomorphic
        result = components_are_isomorphic(self.mapi_components[0],
                                           self.mapi_components[1])
        self.assertTrue(result)

        # check CH3NH3 and PbI3 are not isomorphic
        result = components_are_isomorphic(self.mapi_components[0],
                                           self.mapi_components[4])
        self.assertFalse(result)

    def test_get_structure_inequiv_components(self):
        inequiv_comp = get_structure_inequiv_components(
            self.mapi_components, use_structure_graph=False)

        self.assertEqual(len(inequiv_comp), 2)
        self.assertEqual(inequiv_comp[0]['count'], 4)
        # TODO: Uncomment out inequivalent ids when functionality implemented
        # self.assertEqual(inequiv_comp[0]['inequivalent_ids'],
        #                  (0, 8, 12, 44, 20, 28))

        self.assertEqual(inequiv_comp[1]['count'], 1)
        # self.assertEqual(inequiv_comp[1]['inequivalent_ids'],
        #                  (32, 24, 36))

        # test using graph/fingerprint matching
        inequiv_comp = get_structure_inequiv_components(
            self.mapi_components, use_structure_graph=True)

        self.assertEqual(len(inequiv_comp), 2)
        self.assertEqual(inequiv_comp[0]['count'], 4)
        # self.assertEqual(inequiv_comp[0]['inequivalent_ids'],
        #                  (0, 8, 12, 44, 20, 28))

        self.assertEqual(inequiv_comp[1]['count'], 1)
        # self.assertEqual(inequiv_comp[1]['inequivalent_ids'],
        #                  (32, 24, 36))
