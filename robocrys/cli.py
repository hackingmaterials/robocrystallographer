"""This module contains a script for using robocrys from the command line."""

from __future__ import annotations

import argparse
import sys
import warnings

from pymatgen.core.structure import Structure

try:
    from mp_api.client import MPRestError  # type:ignore[import-untyped]
except ImportError:
    MPRestError = None

from robocrys import StructureCondenser, StructureDescriber, __version__

__author__ = "Alex Ganose"
__maintainer__ = "Alex Ganose"
__email__ = "aganose@lbl.gov"
__date__ = "December 17, 2018"


def robocrystallographer(
    structure: Structure,
    condenser_kwargs: dict | None = None,
    describer_kwargs: dict | None = None,
) -> str:
    """Gets the robocrystallographer description of a structure.

    Args:
        structure: A structure.
        condenser_kwargs: Keyword arguments that will be passed to
            :obj:`robocrys.condense.StructureCondenser`.
        describer_kwargs: Keyword arguments that will be passed to
            :obj:`robocrys.describe.StructureDescriber`.

    Returns:
        The description.
    """
    condenser_kwargs = condenser_kwargs if condenser_kwargs else {}
    describer_kwargs = describer_kwargs if describer_kwargs else {}

    sc = StructureCondenser(**condenser_kwargs)
    describer = StructureDescriber(**describer_kwargs)

    if not any(hasattr(s, "oxi_state") for s in structure.composition.elements):
        try:
            structure.add_oxidation_state_by_guess(max_sites=-80)
        except ValueError:
            warnings.warn("Could not add oxidation states!")

    condensed_structure = sc.condense_structure(structure)
    description = describer.describe(condensed_structure)
    print(description)
    return description  # type: ignore[return-value]


def _get_parser():
    parser = argparse.ArgumentParser(
        description="robocrystallographer is a tool to generate crystal "
        "structure descriptions",
        epilog="Author: {}, Version: {}, Last updated: {}".format(
            __author__, __version__, __date__
        ),
    )

    parser.add_argument("filename", help="structure file or mpid")
    parser.add_argument(
        "-c",
        "--conventional",
        dest="use_conventional_cell",
        action="store_true",
        help="use the convention cell",
    )
    parser.add_argument(
        "-s",
        "--symmetry",
        action="store_true",
        dest="use_symmetry_equivalent_sites",
        help="use symmetry to determine inequivalent sites",
    )
    parser.add_argument(
        "--symprec", default=0.01, type=float, help="symmetry tolerance"
    )
    parser.add_argument(
        "--no-simplify",
        action="store_false",
        dest="simplify_molecules",
        help="don't simplify molecules when mineral matching",
    )
    parser.add_argument(
        "--no-iupac",
        action="store_false",
        dest="use_iupac_formula",
        help="don't use IUPAC formula ordering",
    )
    parser.add_argument(
        "--no-common-formulas",
        dest="use_common_formulas",
        action="store_false",
        help="don't use common formulas",
    )
    parser.add_argument(
        "--no-mineral",
        dest="describe_mineral",
        action="store_false",
        help="don't describe the mineral information",
    )
    parser.add_argument(
        "--no-makeup",
        dest="describe_component_makeup",
        action="store_false",
        help="don't describe the component makeup",
    )
    parser.add_argument(
        "--no-components",
        dest="describe_components",
        action="store_false",
        help="don't describe the components",
    )
    parser.add_argument(
        "--no-symmetry-labels",
        dest="describe_symmetry_labels",
        action="store_false",
        help="don't describe symmetry labels",
    )
    parser.add_argument(
        "--no-oxi",
        dest="describe_oxidation_states",
        action="store_false",
        help="don't describe oxidation states",
    )
    parser.add_argument(
        "--no-bond",
        dest="describe_bond_lengths",
        action="store_false",
        help="don't describe bond lengths",
    )
    parser.add_argument(
        "--precision",
        metavar="P",
        dest="bond_length_decimal_places",
        default=2,
        type=int,
        help="decimal places for bond lengths",
    )
    parser.add_argument(
        "--distorted-tol",
        metavar="T",
        dest="distorted_tol",
        default=0.6,
        type=float,
        help="order parameter below which sites are distorted",
    )
    parser.add_argument(
        "--anion-polyhedra",
        dest="only_describe_cation_polyhedra_connectivity",
        action="store_true",
        help="describe anion polyhedra connectivity",
    )
    parser.add_argument(
        "--verbose-bonds",
        dest="only_describe_bonds_once",
        action="store_false",
        help="describe bond lengths for each site",
    )
    parser.add_argument(
        "--format",
        dest="fmt",
        default="unicode",
        help="how to format the description (unicode [default], html, latex, raw)",
    )
    parser.add_argument(
        "--api-key",
        help="set the materials project API key. See: "
        "https://materialsproject.org/docs/api",
    )
    return parser


def main():
    args = _get_parser().parse_args()
    args_dict = vars(args)

    condenser_keys = [
        "use_conventional_cell",
        "use_symmetry_equivalent_sites",
        "symprec",
        "use_iupac_formula",
        "use_common_formulas",
    ]
    describer_keys = [
        "describe_mineral",
        "describe_component_makeup",
        "describe_components",
        "describe_symmetry_labels",
        "describe_oxidation_states",
        "describe_bond_lengths",
        "bond_length_decimal_places",
        "distorted_tol",
        "only_describe_cation_polyhedra_connectivity",
        "only_describe_bonds_once",
        "fmt",
    ]

    condenser_kwargs = {key: args_dict[key] for key in condenser_keys}
    describer_kwargs = {key: args_dict[key] for key in describer_keys}

    warnings.filterwarnings("ignore", category=UserWarning, module="pymatgen")
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    try:
        structure = Structure.from_file(args.filename)

    except FileNotFoundError:
        if args.filename[:3] in ["mp-", "mvc"]:

            if MPRestError is None:
                raise ImportError(
                    "Please install the Materials Project API using"
                    "`pip install mp_api`"
                )

            from mp_api.client import MPRester

            mpr = MPRester(args.api_key)

            try:
                structure = mpr.get_structure_by_material_id(args.filename)
            except IndexError:
                print("filename or mp-id not found.")
                sys.exit()
            except MPRestError as e:
                if "API_KEY is not supplied" in str(e):
                    print(
                        "Materials project API key not set. Use the the "
                        "--api-key option.\nSee robocrys -h for more details"
                    )
                    sys.exit()
                else:
                    raise e
        else:
            print(f"structure file '{args.filename}' not found.")
            sys.exit()

    if not structure.is_ordered:
        print(
            "disordered structures are not currently supported by "
            "robocrystallographer... exiting"
        )
        sys.exit()

    robocrystallographer(
        structure, condenser_kwargs=condenser_kwargs, describer_kwargs=describer_kwargs
    )
