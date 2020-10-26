Condensed structure format
==========================

An example of the condensed structure JSON format is given below::

    {
        'mineral': {
            'type': 'Molybdenite',
            'distance': -1,
            'n_species_type_match': True,
            'simplified': False
        },
        'formula': 'MoS2',
        'spg_symbol': 'P6_3/mmc',
        'crystal_system': 'hexagonal',
        'dimensionality': 2,
        'sites': {
            0: {'element': 'Mo4+',
                'geometry': {
                    'distorted': True,
                    'type': 'pentagonal pyramidal'
                },
                'nn': [2, 2, 2, 2, 2, 2],
                'nnn': {'edge': [0, 0, 0, 0, 0, 0]},
                'poly_formula': 'S6',
                'sym_labels': (1,)
            },
            2: {'element': 'S2-',
                'geometry': {
                    'distorted': False,
                    'type': '3-coordinate'
                },
                'nn': [0, 0, 0],
                'nnn': {'corner': [2, 2, 2, ...],
                        'face': [2]},
                'poly_formula': None,
                'sym_labels': (1,)
            }
        },
        'distances': {
            0: {2: [2.42, 2.42, 2.42, ...]},
            2: {0: [2.42, 2.42, 2.42]}
        },
        'angles': {
            0: {0: {'edge': [82.60, 82.60, 82.60, ...]}},
            2: {2: {'corner': [135.20, 82.60, 135.20, ...],
                    'face': [80.70, 80.70, 80.70]}}},
        'nnn_distances': {
            0: {0: {'edge': [3.19, 3.19, 3.19, ...]}},
            2: {2: {'face': [3.13],
                    'corner': [3.19, 4.47, 3.19, ...]}}},
        'components': {
            0: {'formula':
                'MoS2',
                'sites': [0, 2]
                'dimensionality': 2,
                'molecule_name': None,
                'orientation': (0, 0, 1)
            }
        },
        'component_makeup': [0, 0],
        'vdw_heterostructure_info': {
            'intercalant_formulas': [],
            'num_repetitions': 2,
            'repeating_unit': ['MoS2']
        }
    }
