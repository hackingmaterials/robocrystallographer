"""Describe a set of structures, or compare two sets of descriptions.

    python describe_structures.py STRUCTURES OUT
    python describe_structures.py --compare OLD NEW

STRUCTURES is a json.gz file of {name: Structure}. OUT, OLD and NEW are json
files of {name: description}. The compare step prints, as markdown, the
sentences that differ in each changed description. Run both describe steps
with the same PYTHONHASHSEED.
"""

import difflib
import json
import sys


def describe(structures_path, out_path):
    from monty.serialization import loadfn

    from robocrys import StructureCondenser, StructureDescriber
    from robocrys.condense.molecule import MoleculeNamer

    # no online PubChem lookups, so the descriptions do not depend on the network
    MoleculeNamer.get_name_from_pubchem = lambda self, smiles: None

    condenser, describer = StructureCondenser(), StructureDescriber()
    descriptions = {}
    for name, structure in loadfn(structures_path).items():
        try:
            condensed = condenser.condense_structure(structure)
            descriptions[name] = describer.describe(condensed)
        except Exception as exc:
            descriptions[name] = f"Error: {type(exc).__name__}: {exc}"
    with open(out_path, "w") as f:
        json.dump(descriptions, f, indent=1)


def compare(old_path, new_path):
    with open(old_path) as f:
        old = json.load(f)
    with open(new_path) as f:
        new = json.load(f)
    changed = [name for name in old if old[name] != new.get(name)]
    print(f"## {len(changed)} of {len(old)} descriptions change\n")
    for name in changed:
        # show only the sentences that differ
        before = old[name].split(". ")
        after = (new.get(name) or "").split(". ")
        lines = [
            line
            for line in difflib.ndiff(before, after)
            if line.startswith(("- ", "+ "))
        ]
        print(f"### {name}\n\n```diff")
        print("\n".join(lines))
        print("```\n")


if __name__ == "__main__":
    if sys.argv[1] == "--compare":
        compare(sys.argv[2], sys.argv[3])
    else:
        describe(sys.argv[1], sys.argv[2])
