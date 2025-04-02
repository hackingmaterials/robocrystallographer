from __future__ import annotations

import os
from pathlib import Path
import unittest

from monty.json import MontyDecoder
from monty.serialization import loadfn

from robocrys.util import load_condensed_structure_json

TEST_DIR = Path(__file__).absolute().parent / ".." / ".." / "tests"
TEST_FILES_DIR = TEST_DIR / "test_files"


class RobocrysTest(unittest.TestCase):
    """Base test class providing access to common test data."""

    _structures_dir = TEST_FILES_DIR / "structures"
    _condensed_structures_dir = TEST_FILES_DIR / "condensed_structures"

    _test_structures = {}
    for _fn in os.listdir(_structures_dir):
        if ".json.gz" in _fn:
            _test_structures[_fn.split(".")[0]] = loadfn(
                _structures_dir / _fn, cls=MontyDecoder
            )

    _test_condensed_structures = {}
    for _fn in os.listdir(_condensed_structures_dir):
        if ".json.gz" in _fn:
            _test_condensed_structures[_fn.split(".")[0]] = (
                load_condensed_structure_json(_condensed_structures_dir / _fn)
            )

    @classmethod
    def get_structure(cls, name):
        return cls._test_structures[name].copy()

    @classmethod
    def get_condensed_structure(cls, name):
        return cls._test_condensed_structures[name].copy()
