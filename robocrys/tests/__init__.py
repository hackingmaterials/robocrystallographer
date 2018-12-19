import os
import unittest

from monty.json import MontyDecoder
from monty.serialization import loadfn
from robocrys.util import load_condensed_structure_json


class RobocrysTest(unittest.TestCase):
    """Base test class providing access to common test data. """

    _module_dir = os.path.dirname(os.path.abspath(__file__))
    _structures_dir = os.path.join(_module_dir, "structures")
    _condensed_structures_dir = os.path.join(
        _module_dir, "condensed_structures")

    _test_structures = {}
    for _fn in os.listdir(_structures_dir):
        if ".json.gz" in _fn:
            _test_structures[_fn.split(".")[0]] = loadfn(os.path.join(
                _structures_dir, _fn), cls=MontyDecoder)

    _test_condensed_structures = {}
    for _fn in os.listdir(_condensed_structures_dir):
        if ".json.gz" in _fn:
            _test_condensed_structures[_fn.split(".")[0]] = \
                load_condensed_structure_json(os.path.join(
                    _condensed_structures_dir, _fn))

    @classmethod
    def get_structure(cls, name):
        return cls._test_structures[name].copy()

    @classmethod
    def get_condensed_structure(cls, name):
        return cls._test_condensed_structures[name].copy()
