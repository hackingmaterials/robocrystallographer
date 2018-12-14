import argparse
import sys
import warnings
import logging

from pymatgen.core.structure import Structure

from robocrys import StructureCondenser #, Describer
from robocrys.util import RoboPrinter

__author__ = "Alex Ganose"
__version__ = "0.0.1"
__maintainer__ = "Alex Ganose"
__email__ = "aganose@lbl.gov"
__date__ = "October 12, 2018"


def robocrystallographer(structure):
    sc = StructureCondenser(force_conventional_cell=True)
    #describer = Describer()

    try:
        print("adding oxidation states")
        structure.add_oxidation_state_by_guess(max_sites=-80)
    except ValueError:
        print("could not add oxidation states")
        pass

    condensed_structure = sc.condense_structure(structure)
    from pprint import pprint
    rp = RoboPrinter({float: "%.2f"}, width=250)
    rp.pprint(condensed_structure)
    # description = describer.describe(condensed_structure)

    # print(description)
    # return description


def _get_parser():
    parser = argparse.ArgumentParser(description="""
    robocrystallographer is a tool to generate structure descriptions""",
                                     epilog="""
    Author: {}
    Version: {}
    Last updated: {}""".format(__author__, __version__, __date__))

    parser.add_argument('filename', type=str,
                        help="structure file or mpid")

    # TODO: Options e.g. NN algo, Use fingerprint matching, other mineral
    # properties, force conv_cell
    return parser


def main():
    args = _get_parser().parse_args()

    # logging.basicConfig(filename='robocrys.log', level=logging.INFO,
    #                     filemode='w', format='%(message)s')
    # console = logging.StreamHandler()
    # logging.info(" ".join(sys.argv[:]))
    # logging.getLogger('').addHandler(console)

    warnings.filterwarnings("ignore", category=UserWarning,
                            module="pymatgen")
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    try:
        structure = Structure.from_file(args.filename)

    except FileNotFoundError:
        from pymatgen.ext.matproj import MPRester

        mpr = MPRester()

        try:
            structure = mpr.get_entry_by_material_id(
                args.filename, inc_structure='final').structure
        except IndexError:
            logging.error("filename or mp-id not found.")
            sys.exit()

    robocrystallographer(structure)


if __name__ == "__main__":
    main()
