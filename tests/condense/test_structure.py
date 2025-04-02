from __future__ import annotations

from pytest import approx

from robocrys.condense.condenser import StructureCondenser
from robocrys.util.tests import RobocrysTest


class TestStructureCondenser(RobocrysTest):
    """Class to test mineral matching functionality."""

    def setUp(self):
        self.tin_dioxide = self.get_structure("SnO2")

    def test_init(self):
        sc = StructureCondenser()
        assert sc is not None

    def test_condense_structure_default(self):
        """Test structure condensing."""
        sc = StructureCondenser()
        data = sc.condense_structure(self.tin_dioxide)

        assert data["mineral"]["type"] == "Rutile"
        assert data["mineral"]["simplified"] is False
        assert data["mineral"]["n_species_type_match"] is True
        assert data["mineral"]["distance"] == -1

        assert data["spg_symbol"] == "P4_2/mnm"
        assert data["crystal_system"] == "tetragonal"
        assert data["dimensionality"] == 3

        # check the right number of sites and that the site data is correct
        # for one site
        assert len(data["sites"].keys()) == 2
        assert data["sites"][0]["element"] == "Sn4+"
        assert data["sites"][0]["geometry"]["type"] == "octahedral"
        assert data["sites"][0]["geometry"]["likeness"] == approx(0.9349817375244279)
        assert len(data["sites"][0]["nn"]) == 6
        assert len(data["sites"][0]["nnn"]["corner"]) == 8
        assert len(data["sites"][0]["nnn"]["edge"]) == 2
        assert data["sites"][0]["poly_formula"] == "O6"
        assert data["sites"][0]["sym_labels"] == (1,)

        # check distances
        assert len(data["distances"][0][2]) == 6
        assert len(data["distances"][2][0]) == 3
        assert data["distances"][0][2][0] == approx(2.092210570377848)

        # check angles
        assert len(data["angles"][0][0]["corner"]) == 8
        assert len(data["angles"][0][0]["edge"]) == 4
        assert data["angles"][0][0]["edge"][0] == approx(101.62284671698572)

        # check nnn distances
        assert len(data["nnn_distances"][0][0]["corner"]) == 8
        assert len(data["nnn_distances"][0][0]["edge"]) == 2
        assert data["nnn_distances"][0][0]["edge"][0] == approx(3.24322132)

        # check components
        assert data["components"][0]["dimensionality"] == 3
        assert data["components"][0]["orientation"] is None
        assert data["components"][0]["formula"] == "SnO2"
        assert data["components"][0]["molecule_name"] is None
        assert data["components"][0]["sites"] == [0, 0, 2, 2, 2, 2]
        assert data["component_makeup"] == [0]

        # check vdw heterostructure information doesn't exist
        assert data["vdw_heterostructure_info"] is None

    def test_condense_structure_sym(self):
        """Test nothing changes when we use symmetry to reduce components."""
        sc = StructureCondenser(use_symmetry_equivalent_sites=True)
        data = sc.condense_structure(self.tin_dioxide)

        # check the right number of sites and that the site data is correct
        # for one site
        assert len(data["sites"].keys()) == 2
        assert data["sites"][0]["element"] == "Sn4+"
        assert data["sites"][0]["geometry"]["type"] == "octahedral"
        assert data["sites"][0]["geometry"]["likeness"] == approx(0.9349817375244279)
        assert len(data["sites"][0]["nn"]) == 6
        assert len(data["sites"][0]["nnn"]["corner"]) == 8
        assert len(data["sites"][0]["nnn"]["edge"]) == 2
        assert data["sites"][0]["poly_formula"] == "O6"
        assert data["sites"][0]["sym_labels"] == (1,)

        # check distances
        assert len(data["distances"][0][2]) == 6
        assert len(data["distances"][2][0]) == 3
        assert data["distances"][0][2][0] == approx(2.092210570377848)

        # check angles
        assert len(data["angles"][0][0]["corner"]) == 8
        assert len(data["angles"][0][0]["edge"]) == 4
        assert data["angles"][0][0]["edge"][0] == approx(101.62284671698572)

        # check nnn distances
        assert len(data["nnn_distances"][0][0]["corner"]) == 8
        assert len(data["nnn_distances"][0][0]["edge"]) == 2
        assert data["nnn_distances"][0][0]["edge"][0] == approx(3.24322132)

        # check components
        assert data["components"][0]["dimensionality"] == 3
        assert data["components"][0]["orientation"] is None
        assert data["components"][0]["formula"] == "SnO2"
        assert data["components"][0]["molecule_name"] is None
        assert data["components"][0]["sites"] == [0, 0, 2, 2, 2, 2]
        assert data["component_makeup"] == [0]
