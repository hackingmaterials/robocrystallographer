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
import re
from collections import defaultdict
from typing import Union, Dict, List, Any

from monty.json import MontyDecoder
from monty.serialization import loadfn
from pkg_resources import resource_filename

from pymatgen import Element
from pymatgen.core.periodic_table import get_el_sp, Species
from pymatgen.util.string import latexify_spacegroup

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
    'trigonal pyramidal': 'trigonal pyramid',
    'square pyramidal': 'square pyramid',
    'trigonal bipyramidal': 'trigonal bipyramid',
    'pentagonal pyramidal': 'pentagonal pyramid',
    'hexagonal pyramidal': 'hexagonal pyramid',
    'pentagonal bipyramidal': 'pentagonal bipyramid',
    'hexagonal bipyramidal': 'hexagonal bipyramid',
    'cuboctahedral': 'cuboctahedra'}

polyhedra_plurals: Dict[str, str] = {
    'octahedra': 'octahedra',
    'tetrahedra': 'tetrahedra',
    'trigonal pyramid': 'trigonal pyramids',
    'square pyramid': 'square pyramids',
    'trigonal bipyramid': 'trigonal bipyramids',
    'pentagonal pyramid': 'pentagonal pyramids',
    'hexagonal pyramid': 'hexagonal pyramids',
    'pentagonal bipyramid': 'pentagonal bipyramids',
    'hexagonal bipyramid': 'hexagonal bipyramids',
    'cuboctahedra': 'cuboctahedra'}

dimensionality_to_shape: Dict[int, str] = {
    3: 'framework', 2: 'sheet', 1: 'ribbon', 0: 'cluster'}


def get_el(obj: Union[Element, Species, str, int]) -> str:
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
    elif isinstance(obj, Species):
        return obj.element.name
    elif isinstance(obj, int):
        return Element.from_Z(obj).name
    else:
        raise ValueError("Unsupported element type: {}.".format(type(obj)))


def get_formatted_el(element: str,
                     sym_label: str,
                     use_oxi_state: bool = True,
                     use_sym_label: bool = True,
                     fmt: str = "raw"):
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

    if isinstance(specie, Species):
        oxi_state = specie.oxi_state
        sign = '+' if oxi_state > 0 else '-'
        if oxi_state == 0:
            oxi_state = None
        elif oxi_state % 1 == 0:
            oxi_state = '{:d}{}'.format(int(abs(oxi_state)), sign)
        else:
            oxi_state = '{:+.2f}{}'.format(abs(oxi_state), sign)
    else:
        oxi_state = None

    formatted_element = specie.name

    if use_sym_label:
        formatted_element += sym_label

    if use_oxi_state and oxi_state:
        if fmt == "latex":
            oxi_state = "^{{{}}}".format(oxi_state)

        elif fmt == "unicode":
            oxi_state = superscript_number(oxi_state)

        elif fmt == "html":
            oxi_state = "<sup>{}</sup>".format(oxi_state)

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

    if '.' in string:
        # no unicode period exists
        return string

    subscript_unicode_map = {0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵',
                             6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹', "-": "⁻", "+": "⁺"}

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
    subscript_unicode_map = {0: "₀", 1: "₁", 2: "₂", 3: "₃", 4: "₄", 5: "₅",
                             6: "₆", 7: "₇", 8: "₈", 9: "₉"}

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
    symbol = re.sub(r"-(\d)", r"{}\1".format(overline), symbol)
    return symbol


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


def load_condensed_structure_json(filename: str) -> Dict[str, Any]:
    """Load condensed structure data from a file.

    Args:
        filename: The filename.

    Returns:
        The condensed structure data.
    """

    # Json does not support using integeras a dictionary keys, therefore
    # manually convert dictionary keys from str to int if possible.
    def json_keys_to_int(x):
        if isinstance(x, dict):
            return {int(k) if k.isdigit() else k: v for k, v in x.items()}

    return loadfn(filename, cls=MontyDecoder, object_hook=json_keys_to_int)
