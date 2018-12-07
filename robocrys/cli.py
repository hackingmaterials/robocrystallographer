import argparse
import logging
import sys
import warnings

import inflect
from pymatgen.core.structure import Structure

from robocrys import StructureCondenser
from robocrys.description import Describer

__author__ = "Alex Ganose"
__version__ = "0.0.1"
__maintainer__ = "Alex Ganose"
__email__ = "aganose@lbl.gov"
__date__ = "October 12, 2018"

en = inflect.engine()


def robocrystallographer(structure):
    sc = StructureCondenser()
    describer = Describer()

    condensed_structure = sc.condense_structure(structure)
    description = describer.describe(condensed_structure)

    logging.info(description)
    return description


def _get_parser():
    parser = argparse.ArgumentParser(description="""
    robocrystallographer is a tool to generate structure descriptions""",
                                     epilog="""
    Author: {}
    Version: {}
    Last updated: {}""".format(__author__, __version__, __date__))

    parser.add_argument('filename', type=str,
                        help="structure file or mpid")
    return parser


def main():
    args = _get_parser().parse_args()
    logging.basicConfig(filename='robocrys.log', level=logging.INFO,
                        filemode='w', format='%(message)s')
    console = logging.StreamHandler()
    logging.info(" ".join(sys.argv[:]))
    logging.getLogger('').addHandler(console)

    warnings.filterwarnings("ignore", category=UserWarning,
                            module="pymatgen")
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    try:
        structure = Structure.from_file(args.filename)

    except FileNotFoundError:
        from pymatgen.ext.matproj import MPRester

        mpr = MPRester()
        structure = mpr.get_entry_by_material_id(
            args.filename, inc_structure='final').structure

    robocrystallographer(structure)


if __name__ == "__main__":
    main()
