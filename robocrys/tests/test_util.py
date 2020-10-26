import os
import unittest

from pymatgen.core.periodic_table import get_el_sp
from robocrys.util import get_el, get_formatted_el, \
    load_condensed_structure_json, superscript_number, unicodeify_spacegroup, \
    htmlify_spacegroup

test_dir = os.path.dirname(os.path.abspath(__file__))


class TestDescriptionMethods(unittest.TestCase):
    """Class to test utility functions."""

    def test_common_formulas(self):
        from robocrys.util import common_formulas
        self.assertEqual(common_formulas["H6CN"], "CH3NH3")

    def test_connected_geometries(self):
        from robocrys.util import connected_geometries
        self.assertTrue("octahedral" in connected_geometries)

    def test_geometry_to_polyhedra(self):
        from robocrys.util import geometry_to_polyhedra
        self.assertEqual(geometry_to_polyhedra["octahedral"], "octahedra")

    def test_polyhedra_plurals(self):
        from robocrys.util import polyhedra_plurals
        self.assertEqual(polyhedra_plurals["pentagonal pyramid"],
                         "pentagonal pyramids")

    def test_dimensionality_to_shape(self):
        from robocrys.util import dimensionality_to_shape
        self.assertEqual(dimensionality_to_shape[3], "framework")

    def test_superscript_number(self):
        self.assertEqual(superscript_number("+23"), "⁺²³")

    def test_get_al(self):
        """Test getting element names"""
        species = get_el_sp("Sn2+")
        self.assertEqual(get_el(species), "Sn")

        element = get_el_sp("Sn")
        self.assertEqual(get_el(element), "Sn")

        self.assertEqual(get_el(50), "Sn")

    def test_get_formatted_el(self):
        """Test getting formatted element strings."""
        species = get_el_sp("Sn2+")
        form_el = get_formatted_el(species, "")
        self.assertEqual(form_el, "Sn2+")

        element = get_el_sp("Sn")
        form_el = get_formatted_el(element, "")
        self.assertEqual(form_el, "Sn")

        form_el = get_formatted_el(get_el(50), "")
        self.assertEqual(form_el, "Sn")

        form_el = get_formatted_el("Sn2+", "(1,2)", use_oxi_state=True,
                                   use_sym_label=True, fmt="raw")
        self.assertEqual(form_el, "Sn(1,2)2+")

        form_el = get_formatted_el("Sn2+", "(1,2)", use_oxi_state=False,
                                   use_sym_label=True, fmt="raw")
        self.assertEqual(form_el, "Sn(1,2)")

        form_el = get_formatted_el("Sn2+", "(1,2)", use_oxi_state=True,
                                   use_sym_label=False, fmt="raw")
        self.assertEqual(form_el, "Sn2+")

        form_el = get_formatted_el("Sn2+", "(1,2)", use_oxi_state=False,
                                   use_sym_label=False, fmt="raw")
        self.assertEqual(form_el, "Sn")

        form_el = get_formatted_el("Sn2+", "(1,2)", use_oxi_state=True,
                                   use_sym_label=True, fmt="latex")
        self.assertEqual(form_el, "Sn(1,2)^{2+}")

        form_el = get_formatted_el("Sn2+", "(1,2)", use_oxi_state=True,
                                   use_sym_label=True, fmt="html")
        self.assertEqual(form_el, "Sn(1,2)<sup>2+</sup>")

        form_el = get_formatted_el("Sn2+", "(1,2)", use_oxi_state=True,
                                   use_sym_label=True, fmt="unicode")
        self.assertEqual(form_el, "Sn(1,2)²⁺")

    def test_unicodeify_spacegroup(self):
        spg_symbol = unicodeify_spacegroup("P-42_1m")
        self.assertEqual(spg_symbol, "P̅42₁m")

    def test_htmlify_spacegroup(self):
        spg_symbol = htmlify_spacegroup("P-42_1m")
        self.assertEqual(spg_symbol, "P̅42<sub>1</sub>m")

    def test_load_condense_structure_json(self):
        condensed_structure = load_condensed_structure_json(
            os.path.join(test_dir, 'condensed_structures', 'SnO2.json.gz'))

        # check that the site keys are correctly converted from str to int
        site_keys = list(condensed_structure['sites'].keys())
        self.assertIsInstance(site_keys[0], int)

        component_keys = list(condensed_structure['components'].keys())
        self.assertIsInstance(component_keys[0], int)
