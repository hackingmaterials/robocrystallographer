"""
This script processes the Pubchem entries stored in a mongodb database for use
in robocrystallographer.

The aim to to produce a dictionary mapping SMILES string to compound name.
"""

from maggma.advanced_stores import MongograntStore
from monty.serialization import dumpfn
from tqdm import tqdm

name_keys = ['name_cas', 'name_iupac', 'name_traditional']

mp_pubchem = MongograntStore("rw:knowhere.lbl.gov/mp_pubchem", "mp_pubchem",
                             key="pubchem_id")
mp_pubchem.connect()
coll = mp_pubchem.collection

total_mols = coll.count()
pbar = tqdm(total=total_mols, desc="process")

results = coll.find()

smiles_to_names = {}
n_skipped = 0
for mol in results:
    skip = False
    smiles = mol['smiles_can']
    names = {key: mol[key] for key in name_keys if key in mol}

    if int(mol['xyz'].split('\n')[0]) > 50:
        # skip if too many atoms
        skip = True

    if not names:
        skip = True

    if not skip:
        smiles_to_names[smiles] = names
    else:
        n_skipped += 1

    pbar.update()

pbar.close()

print("total skipped: {}".format(n_skipped))
print("total saved: {}".format(len(smiles_to_names.keys())))

dumpfn(smiles_to_names, "smiles_to_names.json")
