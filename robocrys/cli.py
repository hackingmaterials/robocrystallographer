import sys
import logging
import warnings
import argparse
import inflect

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from robocrys import MineralMatcher
from robocrys import SiteDescriber
from robocrys.fragment import get_structure_fragments

__author__ = "Alex Ganose"
__version__ = "0.0.1"
__maintainer__ = "Alex Ganose"
__email__ = "aganose@lbl.gov"
__date__ = "October 12, 2018"

en = inflect.engine()


def robocrystallographer(structure):

    mineral_matcher = MineralMatcher()
    mineral = mineral_matcher.get_best_mineral_name(structure)

    logging.info("{} is {} structured. ".format(
        structure.composition.reduced_formula, mineral))

    fragments = get_structure_fragments(structure)
    dimensionality = max([f['dimensionality'] for f in fragments])

    desc = "The structure is {} dimensional ".format(
        en.number_to_words(dimensionality))

    if dimensionality < 3:
        comps = [f['structure'].composition.reduced_formula for f in fragments]

        if dimensionality == 2:
            shape = "sheet"
        elif dimensionality == 1:
            shape = "ribbon"
        else:
            shape = "cluster"

        desc += "and consists of {} {} {} oriented in the {} direction. ".format(
            en.number_to_words(len(comps)), comps[0],
            en.plural(shape, len(comps)),
            tuple(map(int, fragments[0]['orientation'])))

        desc += "In each {}, ".format(shape)
    else:
        desc += ". "

    logging.info(desc)

    sga = SpacegroupAnalyzer(structure)
    structure = sga.get_symmetrized_structure()

    site_describer = SiteDescriber(structure)
    for i, list_sites in enumerate(structure.equivalent_indices):
        # very rough way of not overloading with information about bond lengths
        bond_lengths = i == len(structure.equivalent_indices) - 1
        logging.info(site_describer.get_site_description(
            list_sites[0], describe_bond_lengths=bond_lengths))


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
