from robocrys.condense.structure import StructureCondenser
from robocrys.util import RobocrysTest


class TestStructureCondenser(RobocrysTest):
    """Class to test mineral matching functionality."""

    def setUp(self):
        self.tin_dioxide = self.get_structure("SnO2")

    def test_init(self):
        sc = StructureCondenser()
        self.assertNotEqual(sc, None)

    def test_condense_structure_default(self):
        """Test structure condensing."""
        sc = StructureCondenser()
        data = sc.condense_structure(self.tin_dioxide)

        self.assertEqual(data['mineral']['type'], 'Rutile')
        self.assertEqual(data['mineral']['simplified'], False)
        self.assertEqual(data['mineral']['n_species_type_match'], True)
        self.assertEqual(data['mineral']['distance'], -1)

        self.assertEqual(data['spg_symbol'], 'P4_2/mnm')
        self.assertEqual(data['crystal_system'], 'tetragonal')
        self.assertEqual(data['dimensionality'], 3)

        # check the right number of sites and that the site data is correct
        # for one site
        self.assertEqual(len(data['sites'].keys()), 2)
        self.assertEqual(data['sites'][0]['element'], 'Sn4+')
        self.assertEqual(data['sites'][0]['geometry']['type'], 'octahedral')
        self.assertAlmostEqual(data['sites'][0]['geometry']['likeness'],
                               0.9349817375244279)
        self.assertEqual(len(data['sites'][0]['nn']), 6)
        self.assertEqual(len(data['sites'][0]['nnn']['corner']), 8)
        self.assertEqual(len(data['sites'][0]['nnn']['edge']), 2)
        self.assertEqual(data['sites'][0]['poly_formula'], 'O6')
        self.assertEqual(data['sites'][0]['sym_labels'], (1, ))

        # check distances
        self.assertEqual(len(data['distances'][0][2]), 6)
        self.assertEqual(len(data['distances'][2][0]), 3)
        self.assertAlmostEqual(data['distances'][0][2][0], 2.0922101061490546)

        # check angles
        self.assertEqual(len(data['angles'][0][0]['corner']), 8)
        self.assertEqual(len(data['angles'][0][0]['edge']), 4)
        self.assertAlmostEqual(data['angles'][0][0]['edge'][0],
                               101.62287790513848)

        # check components
        self.assertEqual(data['components'][0]['dimensionality'], 3)
        self.assertEqual(data['components'][0]['orientation'], None)
        self.assertEqual(data['components'][0]['formula'], 'SnO2')
        self.assertEqual(data['components'][0]['molecule_name'], None)
        self.assertEqual(data['components'][0]['sites'], [0, 2])
        self.assertEqual(data['component_makeup'], [0])

        # check vdw heterostructure information doesn't exist
        self.assertEqual(data['vdw_heterostructure_info'], None)

    def test_condense_structure_sym(self):
        """Test nothing changes when we use symmetry to reduce components."""
        sc = StructureCondenser(use_symmetry=True)
        data = sc.condense_structure(self.tin_dioxide)

        # check the right number of sites and that the site data is correct
        # for one site
        self.assertEqual(len(data['sites'].keys()), 2)
        self.assertEqual(data['sites'][0]['element'], 'Sn4+')
        self.assertEqual(data['sites'][0]['geometry']['type'], 'octahedral')
        self.assertAlmostEqual(data['sites'][0]['geometry']['likeness'],
                               0.9349817375244279)
        self.assertEqual(len(data['sites'][0]['nn']), 6)
        self.assertEqual(len(data['sites'][0]['nnn']['corner']), 8)
        self.assertEqual(len(data['sites'][0]['nnn']['edge']), 2)
        self.assertEqual(data['sites'][0]['poly_formula'], 'O6')
        self.assertEqual(data['sites'][0]['sym_labels'], (1, ))

        # check distances
        self.assertEqual(len(data['distances'][0][2]), 6)
        self.assertEqual(len(data['distances'][2][0]), 3)
        self.assertAlmostEqual(data['distances'][0][2][0], 2.0922101061490546)

        # check angles
        self.assertEqual(len(data['angles'][0][0]['corner']), 8)
        self.assertEqual(len(data['angles'][0][0]['edge']), 4)
        self.assertAlmostEqual(data['angles'][0][0]['edge'][0],
                               101.62287790513848)

        # check components
        self.assertEqual(data['components'][0]['dimensionality'], 3)
        self.assertEqual(data['components'][0]['orientation'], None)
        self.assertEqual(data['components'][0]['formula'], 'SnO2')
        self.assertEqual(data['components'][0]['molecule_name'], None)
        self.assertEqual(data['components'][0]['sites'], [0, 2])
        self.assertEqual(data['component_makeup'], [0])
