from __future__ import unicode_literals
import unittest
import os

from monty.json import MontyDecoder
from monty.serialization import loadfn

"""
Common test support for robocrystallographer test scripts.
"""


class RobocrysTest(unittest.TestCase):
    """Robocrystallographer base test class providing common functions.

    The main use of this class is currently the get_sturcture method,
    which provides access to commonly used structures.
    """

    module_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(module_dir, "tests", "structures")

    test_structures = {}
    for fn in os.listdir(structures_dir):
        test_structures[fn.split(".")[0]] = loadfn(os.path.join(
            structures_dir, fn), cls=MontyDecoder)

    @classmethod
    def get_structure(cls, name):
        return cls.test_structures[name].copy()
