from __future__ import annotations

import unittest

from pymatgen.core.periodic_table import get_el_sp

from robocrys.util import (
    get_el,
    get_formatted_el,
    htmlify_spacegroup,
    load_condensed_structure_json,
    superscript_number,
    unicodeify_spacegroup,
)
from robocrys.util.tests import TEST_FILES_DIR

class TestDescriptionMethods(unittest.TestCase):
    """Class to test utility functions."""

    def test_common_formulas(self):
        from robocrys.util import common_formulas

        assert common_formulas["H6CN"] == "CH3NH3"

    def test_connected_geometries(self):
        from robocrys.util import connected_geometries

        assert "octahedral" in connected_geometries

    def test_geometry_to_polyhedra(self):
        from robocrys.util import geometry_to_polyhedra

        assert geometry_to_polyhedra["octahedral"] == "octahedra"

    def test_polyhedra_plurals(self):
        from robocrys.util import polyhedra_plurals

        assert polyhedra_plurals["pentagonal pyramid"] == "pentagonal pyramids"

    def test_dimensionality_to_shape(self):
        from robocrys.util import dimensionality_to_shape

        assert dimensionality_to_shape[3] == "framework"

    def test_superscript_number(self):
        assert superscript_number("+23") == "⁺²³"

    def test_get_al(self):
        """Test getting element names."""
        species = get_el_sp("Sn2+")
        assert get_el(species) == "Sn"

        element = get_el_sp("Sn")
        assert get_el(element) == "Sn"

        assert get_el(50) == "Sn"

    def test_get_formatted_el(self):
        """Test getting formatted element strings."""
        species = get_el_sp("Sn2+")
        form_el = get_formatted_el(species, "")
        assert form_el == "Sn2+"

        element = get_el_sp("Sn")
        form_el = get_formatted_el(element, "")
        assert form_el == "Sn"

        form_el = get_formatted_el(get_el(50), "")
        assert form_el == "Sn"

        for (species, symm, use_oxi, use_sym, fmt, expected_label) in [
            ("Sn2+", "(1,2)", True, True, "raw", "Sn(1,2)2+"),
            ("Sn2+", "(1,2)", False, True, "raw", "Sn(1,2)"),
            ("Sn2+", "(1,2)", True, False, "raw", "Sn2+"),
            ("Sn2+", "(1,2)", False, False, "raw", "Sn"),
            ("Sn2+", "(1,2)", True, True, "latex", "Sn(1,2)^{2+}"),
            ("Sn2+", "(1,2)", True, True, "html", "Sn(1,2)<sup>2+</sup>"),
            ("Sn2+", "(1,2)", True, True, "unicode", "Sn(1,2)²⁺"),
            ("Sn1.33-","(3)",True, True, "unicode","Sn(3)¹\u1427³³⁻")
        ]:

            assert get_formatted_el(
                species, symm, use_oxi_state=use_oxi, use_sym_label=use_sym, fmt=fmt
            ) == expected_label



    def test_unicodeify_spacegroup(self):
        spg_symbol = unicodeify_spacegroup("P-42_1m")
        assert spg_symbol == "P̅42₁m"

    def test_htmlify_spacegroup(self):
        spg_symbol = htmlify_spacegroup("P-42_1m")
        assert spg_symbol == "P̅42<sub>1</sub>m"

    def test_load_condense_structure_json(self):

        condensed_structure = load_condensed_structure_json(
            TEST_FILES_DIR / "condensed_structures" / "SnO2.json.gz"
        )

        # check that the site keys are correctly converted from str to int
        site_keys = list(condensed_structure["sites"].keys())
        assert isinstance(site_keys[0], int)

        component_keys = list(condensed_structure["components"].keys())
        assert isinstance(component_keys[0], int)
