"""
Miscellaneous utility functions and common data.

Attributes:
    common_formulas: A set of common formulas. The keys to the data are strings
        from :obj:`pymatgen.core.composition.Composition.reduced_formula`.
    connected_geometries: A list of geometries that are considered
        "connectable" polyhedra. E.g. Their face-sharing, edge-sharing, etc
        properties are of interest.
    geometry_to_polyhedra: A mapping from geometry type (e.g. octahedral) to the
        plural polyhedra name (e.g. octahedra).
    dimensionality_to_shape: A mapping from dimensionality to the component
        shape.
"""

from __future__ import unicode_literals

import os
import unittest
from collections import defaultdict
from typing import Union, Dict, List

from monty.json import MontyDecoder
from monty.serialization import loadfn
from pkg_resources import resource_filename

from pymatgen import Element, Specie
from pymatgen.core.periodic_table import get_el_sp

common_formulas: Dict[str, str] = loadfn(
    resource_filename('robocrys.condense', 'formula_db.json.gz'))

connected_geometries: List[str] = [
    'tetrahedral', 'octahedral', 'trigonal pyramidal',
    'square pyramidal', 'trigonal bipyramidal',
    'pentagonal pyramidal', 'hexagonal pyramidal',
    'pentagonal bipyramidal', 'hexagonal bipyramidal',
    'cuboctahedral']

geometry_to_polyhedra: Dict[str, str] = {
    'octahedral': 'octahedra',
    'tetrahedral': 'tetrahedra',
    'trigonal pyramidal': 'trigonal pyramids',
    'square pyramidal': 'square pyramids',
    'trigonal bipyramidal': 'trigonal bipyramids',
    'pentagonal pyramidal': 'pentagonal pyramids',
    'hexagonal pyramidal': 'hexagonal pyramids',
    'pentagonal bipyramidal': 'pentagonal bipyramids',
    'hexagonal bipyramidal': 'hexagonal bipyramids',
    'cuboctahedral': 'cuboctahedra'}

dimensionality_to_shape: Dict[int, str] = {
    3: 'framework', 2: 'sheet', 1: 'ribbon', 0: 'cluster'}


def get_el(obj: Union[Element, Specie, str, int]) -> str:
    """Utility method to get an element str from a symbol, Element, or Specie.

    Args:
        obj: An arbitrary object. Spported objects are Element/Specie objects,
            integers (representing atomic numbers), or strings (element
            symbols or species strings).

    Returns:
        The element as a string.
    """
    if isinstance(obj, str):
        obj = get_el_sp(obj)

    if isinstance(obj, Element):
        return obj.name
    elif isinstance(obj, Specie):
        return obj.element.name
    elif isinstance(obj, int):
        return Element.from_Z(obj).name
    else:
        raise ValueError("Unsupported element type: {}.".format(type(obj)))


def defaultdict_to_dict(dictionary: defaultdict) -> Dict:
    """Recursively convert nested :obj:`defaultdict` to :obj:`dict`.

    Args:
        dictionary: A defaultdict.

    Returns:
        The defaultdict as a :obj:`dict`.
    """
    if isinstance(dictionary, defaultdict):
        dictionary = {k: defaultdict_to_dict(v) for k, v in dictionary.items()}
    return dictionary


def json_keys_to_int(x):
    """Convert dictionary keys from str to int if possible."""
    if isinstance(x, dict):
        return {int(k) if k.isdigit() else k: v for k, v in x.items()}


class RobocrysTest(unittest.TestCase):
    """Base test class providing access to common test data. """

    module_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(module_dir, "tests", "structures")
    condensed_structures_dir = os.path.join(
        module_dir, "tests", "condensed_structures")

    test_structures = {}
    for fn in os.listdir(structures_dir):
        test_structures[fn.split(".")[0]] = loadfn(os.path.join(
            structures_dir, fn), cls=MontyDecoder)

    test_condensed_structures = {}
    for fn in os.listdir(condensed_structures_dir):
        test_condensed_structures[fn.split(".")[0]] = loadfn(os.path.join(
            condensed_structures_dir, fn), cls=MontyDecoder,
            object_hook=json_keys_to_int)

    @classmethod
    def get_structure(cls, name):
        return cls.test_structures[name].copy()

    @classmethod
    def get_condensed_structure(cls, name):
        return cls.test_condensed_structures[name].copy()
