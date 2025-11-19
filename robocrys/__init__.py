"""General imports.

isort:skip_file
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

from robocrys.condense.condenser import StructureCondenser
from robocrys.describe.describer import StructureDescriber
from robocrys.util import common_formulas

__all__ = ["__version__", "StructureDescriber", "StructureCondenser", "common_formulas"]

try:
    __version__ = version("robocrys")
except PackageNotFoundError:
    # package is not installed
    pass
