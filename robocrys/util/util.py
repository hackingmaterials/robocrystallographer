"""Miscellaneous utility functions and common data.

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

from __future__ import annotations

import re
from collections import defaultdict
from typing import Any, TYPE_CHECKING

from monty.json import MontyDecoder
from monty.serialization import loadfn
from importlib.resources import files as import_resource_file
from pymatgen.core.periodic_table import Element, Species, get_el_sp
from pymatgen.util.string import latexify_spacegroup

if TYPE_CHECKING:
    from pathlib import Path


def _get_common_formulas() -> dict[str, str]:
    """Retrieve common formula information from stored data."""
    all_formulas = loadfn(
        str(import_resource_file("robocrys.condense") / "formula_db.json.gz"), cls=None
    )
    all_formulas["alias"].update({k: k for k in all_formulas["no_alias"]})
    return all_formulas["alias"]


common_formulas: dict[str, str] = _get_common_formulas()

connected_geometries: list[str] = [
    "tetrahedral",
    "octahedral",
    "trigonal pyramidal",
    "square pyramidal",
    "trigonal bipyramidal",
    "pentagonal pyramidal",
    "hexagonal pyramidal",
    "pentagonal bipyramidal",
    "hexagonal bipyramidal",
    "cuboctahedral",
]

geometry_to_polyhedra: dict[str, str] = {
    "octahedral": "octahedra",
    "tetrahedral": "tetrahedra",
    "trigonal pyramidal": "trigonal pyramid",
    "square pyramidal": "square pyramid",
    "trigonal bipyramidal": "trigonal bipyramid",
    "pentagonal pyramidal": "pentagonal pyramid",
    "hexagonal pyramidal": "hexagonal pyramid",
    "pentagonal bipyramidal": "pentagonal bipyramid",
    "hexagonal bipyramidal": "hexagonal bipyramid",
    "cuboctahedral": "cuboctahedra",
}

polyhedra_plurals: dict[str, str] = {
    "octahedra": "octahedra",
    "tetrahedra": "tetrahedra",
    "trigonal pyramid": "trigonal pyramids",
    "square pyramid": "square pyramids",
    "trigonal bipyramid": "trigonal bipyramids",
    "pentagonal pyramid": "pentagonal pyramids",
    "hexagonal pyramid": "hexagonal pyramids",
    "pentagonal bipyramid": "pentagonal bipyramids",
    "hexagonal bipyramid": "hexagonal bipyramids",
    "cuboctahedra": "cuboctahedra",
}

dimensionality_to_shape: dict[int, str] = {
    3: "framework",
    2: "sheet",
    1: "ribbon",
    0: "cluster",
}


def get_el(obj: Element | Species | str | int) -> str:
    """Utility method to get an element str from a symbol, Element, or Specie.

    Args:
        obj: An arbitrary object. Supported objects are Element/Species objects,
            integers (representing atomic numbers), or strings (element
            symbols or species strings).

    Returns:
        The element as a string.
    """
    if isinstance(obj, str):
        obj = get_el_sp(obj)

    if isinstance(obj, Element):
        return obj.name
    if isinstance(obj, Species):
        return obj.element.name
    if isinstance(obj, int):
        return Element.from_Z(obj).name
    raise ValueError(f"Unsupported element type: {type(obj)}.")


def get_formatted_el(
    element: str,
    sym_label: str,
    use_oxi_state: bool = True,
    use_sym_label: bool = True,
    fmt: str = "raw",
):
    """Formats an element string.

    Performs a variety of functions, including:

    - Changing "Sn+0" to "Sn".
    - Inserting the symmetry label between the element and oxidation state, if
        required.
    - Removing the oxidation state if required.
    - Latexifying the element and oxidation state.
    - Unicodeifying the element and oxidation state.
    - Converting the element and oxidation state to html.

    Args:
        element: The element string (possibly including the oxidation state.
            E.g. "Sn" or "Sn2+".
        sym_label: The symmetry label. E.g. "(1)"
        use_oxi_state: Whether to include the oxidation state, if present.
        use_sym_label: Whether to use the symmetry label.
        fmt: How to format the element strings. Options are:

            - "raw" (default): Don't apply special formatting (e.g. "SnO2").
            - "unicode": Format super/subscripts using unicode characters
              (e.g. SnO₂).
            - "latex": Use LaTeX markup for formatting (e.g. "SnO$_2$").
            - "html": Use html markup for formatting (e.g. "SnO<sub>2</sub>").

    Returns:
        The formatted element string.
    """
    specie = get_el_sp(element)

    # NB: the `str(specie)` is only to make mypy happy:
    # `get_el_sp` can return a `DummySpecies` which has
    # no `name` but has `__str__`
    formatted_element = specie.name or str(specie)

    if use_sym_label:
        formatted_element += sym_label

    if (
        use_oxi_state
        and isinstance(specie, Species)
        and (_oxi_state := specie.oxi_state)
        and _oxi_state != 0
    ):
        sign = "+" if _oxi_state > 0 else "-"
        if _oxi_state % 1 == 0:
            oxi_state = f"{int(abs(_oxi_state)):d}{sign}"
        else:
            oxi_state = f"{abs(_oxi_state):.2f}{sign}"
        if fmt == "latex":
            oxi_state = f"^{{{oxi_state}}}"

        elif fmt == "unicode":
            oxi_state = superscript_number(oxi_state)

        elif fmt == "html":
            oxi_state = f"<sup>{oxi_state}</sup>"

        formatted_element += oxi_state

    return formatted_element


def superscript_number(string):
    """Converts a string containing numbers to superscript.

    Will only convert the numbers 0-9, and the + and - characters.

    Args:
        string: A string containing the numbers 0-9 or +/- characters.

    Returns:
        The superscript string.
    """
    subscript_unicode_map = {
        0: "⁰",
        1: "¹",
        2: "²",
        3: "³",
        4: "⁴",
        5: "⁵",
        6: "⁶",
        7: "⁷",
        8: "⁸",
        9: "⁹",
        "-": "⁻",
        "+": "⁺",
        ".": "\u1427",
    }

    for original_subscript, subscript_unicode in subscript_unicode_map.items():
        string = string.replace(str(original_subscript), subscript_unicode)

    return string


def unicodeify_spacegroup(spacegroup_symbol: str) -> str:
    """Formats a spacegroup using unicode symbols.

    E.g. Fd-3m -> Fd̅3m

    Args:
        spacegroup_symbol: A spacegroup symbol.

    Returns:
        The unicode formatted spacegroup symbol.
    """
    subscript_unicode_map = {
        0: "₀",
        1: "₁",
        2: "₂",
        3: "₃",
        4: "₄",
        5: "₅",
        6: "₆",
        7: "₇",
        8: "₈",
        9: "₉",
    }

    symbol = latexify_spacegroup(spacegroup_symbol)

    for number, unicode_number in subscript_unicode_map.items():
        symbol = symbol.replace("$_{" + str(number) + "}$", unicode_number)

    overline = "\u0305"  # u"\u0304" (macron) is also an option

    symbol = symbol.replace("$\\overline{", overline)
    symbol = symbol.replace("$", "")
    symbol = symbol.replace("{", "")
    symbol = symbol.replace("}", "")

    return symbol


def htmlify_spacegroup(spacegroup_symbol: str) -> str:
    """Formats a spacegroup using unicode symbols.

    E.g. P-42_1m -> P̅42<sub>1</sub>m

    Args:
        spacegroup_symbol: A spacegroup symbol.

    Returns:
        The html formatted spacegroup symbol.
    """
    overline = "\u0305"  # u"\u0304" (macron) is also an option
    symbol = re.sub(r"_(\d+)", r"<sub>\1</sub>", spacegroup_symbol)
    symbol = re.sub(r"-(\d)", rf"{overline}\1", symbol)
    return symbol


def defaultdict_to_dict(dictionary: defaultdict) -> dict:
    """Recursively convert nested :obj:`defaultdict` to :obj:`dict`.

    Args:
        dictionary: A defaultdict.

    Returns:
        The defaultdict as a :obj:`dict`.
    """
    if isinstance(dictionary, defaultdict):
        return {k: defaultdict_to_dict(v) for k, v in dictionary.items()}
    return dictionary


def load_condensed_structure_json(filename: str | Path) -> dict[str, Any]:
    """Load condensed structure data from a file.

    Args:
        filename: The filename.

    Returns:
        The condensed structure data.
    """

    # JSON does not support using integers a dictionary keys, therefore
    # manually convert dictionary keys from str to int if possible.
    def json_keys_to_int(x: Any) -> Any:
        if isinstance(x, dict):
            return {
                int(k) if k.isdigit() else k: json_keys_to_int(v) for k, v in x.items()
            }
        return x

    # For some reason, specifying `object_hook = json_keys_to_int` in `loadfn`
    # doesn't seem to work. This does reliably:
    return json_keys_to_int(loadfn(filename, cls=MontyDecoder))
