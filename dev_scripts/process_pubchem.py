"""
This script processes the Pubchem entries stored in a mongodb database for use
in robocrystallographer.

The aim to to produce a dictionary mapping SMILES string to compound name.

The size of the dataset is tuned by the max_atoms parameter.

- max_atoms=60 gives a 2.1 GB gzipped dataset.
"""

from maggma.advanced_stores import MongograntStore
from monty.serialization import dumpfn
import pebble
from tqdm import tqdm
import bson

def process_batch(batch):
    batch = bson.decode_all(batch)

    mols = []
    for mol in batch:
        mol_dict = {}
        try:
            smiles = mol['smiles_can']
            names = {key: mol[key] for key in name_keys if key in mol}

            if int(mol['xyz'].split('\n')[0]) > max_atoms:
                # skip if too many atoms
                raise KeyError("atoms")

            if not names:
                raise KeyError("names")

            mol_dict[smiles] = names
            mols.append(mol_dict)
        except KeyError:
            mols.append(False)

    return mols


def task_done(future):
    results.append(future.result())
    pbar.update()


max_atoms = 12
batch_size = 5000
name_keys = ['name_iupac', 'name_traditional']

mp_pubchem = MongograntStore("rw:knowhere.lbl.gov/mp_pubchem", "mp_pubchem",
                             key="pubchem_id")
mp_pubchem.connect()
coll = mp_pubchem.collection

total_mols = coll.count()
batches = coll.find_raw_batches(batch_size=batch_size)
pbar = tqdm(total=total_mols/batch_size, desc="process")

results = []

with pebble.ProcessPool() as pool:
    for batch in batches:
        f = pool.schedule(process_batch, args=(batch, ))
        f.add_done_callback(task_done)

pbar.close()

n_skipped = 0
smiles_to_names = {}
for processed_batch in results:
    for mol in processed_batch:
        if not mol:
            n_skipped += 1
        else:
            smiles_to_names.update(mol)

print("total skipped: {}".format(n_skipped))
print("total saved: {}".format(len(smiles_to_names.keys())))

dumpfn(smiles_to_names, "smiles_to_names.json")
