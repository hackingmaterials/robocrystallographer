"""Pick 48 Materials Project structures for checking robocrys descriptions.

The first are the mp-ids from robocrys issues (mp-24850 and mp-504366 are no
longer on MP, so 8 remain). The other 40 are stable, experimentally
known materials spread over crystal system, number of sites, chemistry, and
dimensionality. Needs MP_API_KEY in the environment.

    python pick_check_structures.py check_structures.json.gz
"""

import random
import sys
import warnings
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

from monty.serialization import dumpfn
from mp_api.client import MPRester
from pymatgen.analysis.dimensionality import get_dimensionality_larsen
from pymatgen.analysis.local_env import CrystalNN

warnings.filterwarnings("ignore")

# mp-ids from robocrys issues #4, #6, #7, #8, #11, #14, #16, #17, #22 and the README
ISSUE_IDS = [
    "mp-13",
    "mp-1960",
    "mp-82",
    "mp-975272",
    "mp-7631",
    "mp-24850",
    "mp-18820",
    "mp-504366",
    "mp-744259",
    "mp-856",
]
N_OTHER = 40
SITE_BINS = [(1, 4), (5, 12), (13, 30), (31, 60)]
HALOGENS = {"F", "Cl", "Br", "I"}
CHALCOGENS = {"S", "Se", "Te"}
PNICTOGENS = {"N", "P", "As", "Sb"}
NONMETALS = {
    "H",
    "B",
    "C",
    "N",
    "O",
    "Si",
    "P",
    "S",
    "Se",
    "Te",
    "F",
    "Cl",
    "Br",
    "I",
    "As",
    "Sb",
}


def chemistry(elements):
    if "O" in elements:
        return "oxide"
    if elements & HALOGENS:
        return "halide"
    if elements & CHALCOGENS:
        return "chalcogenide"
    if elements & PNICTOGENS:
        return "pnictide"
    if not elements & NONMETALS:
        return "intermetallic"
    return "other"


def site_bin(nsites):
    return next(i for i, (lo, hi) in enumerate(SITE_BINS) if lo <= nsites <= hi)


def dimensionality(structure):
    return get_dimensionality_larsen(CrystalNN().get_bonded_structure(structure))


def main(out_path):
    with MPRester() as mpr:
        docs = mpr.materials.summary.search(
            energy_above_hull=(0, 0),
            theoretical=False,
            num_sites=(1, 60),
            fields=["material_id", "formula_pretty", "symmetry", "nsites", "elements"],
            num_chunks=40,
            chunk_size=1000,
        )
        print(f"{len(docs)} stable, experimentally known materials with 1-60 sites")

        # one candidate per (crystal system, site bin, chemistry), picked at random
        groups = defaultdict(list)
        for d in docs:
            if str(d.material_id) in ISSUE_IDS or d.nsites > 60:
                continue
            elements = {el.symbol for el in d.elements}
            system = getattr(
                d.symmetry.crystal_system, "value", d.symmetry.crystal_system
            )
            key = (str(system), site_bin(d.nsites), chemistry(elements))
            groups[key].append(str(d.material_id))
        rng = random.Random(0)
        candidates = {
            key: rng.choice(sorted(ids)) for key, ids in sorted(groups.items())
        }
        print(f"{len(candidates)} candidates, one per group")

        ids = ISSUE_IDS + list(candidates.values())
        found = mpr.materials.summary.search(
            material_ids=ids,
            fields=["material_id", "formula_pretty", "structure"],
            num_chunks=1,
            chunk_size=1000,
        )
        structures = {str(d.material_id): d.structure for d in found}
        formulas = {str(d.material_id): d.formula_pretty for d in found}

        issue_ids = [i for i in ISSUE_IDS if i in structures]
        gone = sorted(set(ISSUE_IDS) - set(issue_ids))
        print(f"issue mp-ids not on MP any more: {gone}")

    cand_ids = [i for i in candidates.values() if i in structures]
    with Pool() as pool:
        dims = dict(
            zip(cand_ids, pool.map(dimensionality, [structures[i] for i in cand_ids]))
        )

    # take every 0D, 1D, and 2D candidate first (they are rare), then fill up
    # round-robin over crystal systems, so each system gets a similar share
    key_of = {mp_id: key for key, mp_id in candidates.items()}
    picked = [i for i in cand_ids if dims[i] < 3][:12]
    by_system = defaultdict(list)
    for i in cand_ids:
        if i not in picked:
            by_system[key_of[i][0]].append(i)
    for ids_ in by_system.values():
        rng.shuffle(ids_)
    while len(picked) < N_OTHER and any(by_system.values()):
        for system in sorted(by_system):
            if by_system[system] and len(picked) < N_OTHER:
                picked.append(by_system[system].pop())

    final = issue_ids + [i for i in picked if i not in issue_ids]
    for i in final:
        info = (
            "issue"
            if i in issue_ids
            else "{} {} {} {}D".format(
                key_of[i][0], SITE_BINS[key_of[i][1]], key_of[i][2], dims[i]
            )
        )
        print(f"{i:>12}  {formulas[i]:<20} {len(structures[i]):>3} sites  {info}")
    dumpfn({i: structures[i] for i in final}, out_path)
    size = Path(out_path).stat().st_size // 1000
    print(f"wrote {len(final)} structures to {out_path} ({size} kB)")


if __name__ == "__main__":
    main(sys.argv[1])
