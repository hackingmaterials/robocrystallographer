"""
This module implements a script for extracting molecule names, SMILES and xyz
files from the Pubchem database.

It uses multiprocessing to process many files simultaneously.

WARNING: While in progress the disk usage maybe fairly large (100's GB).

Expects a db.json file in the current directory detailing the connection details
for a mongodb database where the data will be sorted. The json file should
contain the keys:

- host
- port
- username
- password

The data will be inserted into the "pubchem" collection in the "pubchem"
database in the mongo database.
"""

import ftplib
import pebble
import os

import pybel
from pebble import ProcessExpired
from pymongo import MongoClient

from tqdm import tqdm
from monty.serialization import loadfn
from concurrent.futures import TimeoutError

keys = ['PUBCHEM_OPENEYE_CAN_SMILES', 'PUBCHEM_OPENEYE_ISO_SMILES',
        'PUBCHEM_IUPAC_CAS_NAME', 'PUBCHEM_IUPAC_NAME',
        'PUBCHEM_IUPAC_TRADITIONAL_NAME', 'PUBCHEM_IUPAC_SYSTEMATIC_NAME',
        'PUBCHEM_IUPAC_OPENEYE_NAME']

key_map = {'PUBCHEM_OPENEYE_CAN_SMILES': 'smiles_can',
           'PUBCHEM_OPENEYE_ISO_SMILES': 'smiles_iso',
           'PUBCHEM_IUPAC_CAS_NAME': 'name_cas',
           'PUBCHEM_IUPAC_NAME': 'name_iupac',
           'PUBCHEM_IUPAC_TRADITIONAL_NAME': 'name_traditional',
           'PUBCHEM_IUPAC_SYSTEMATIC_NAME': 'name_systematic',
           'PUBCHEM_IUPAC_OPENEYE_NAME': 'name_openeye'}

tmp_dir = "tmp"

db_settings = loadfn("db.json")
db = MongoClient(**db_settings, connect=True)


def process_sdf_file(filename):
    skipped = 0
    pubchem_molecules = []
    for i, mol in enumerate(pybel.readfile('sdf', filename)):
        try:
            pubchem_id = mol.data['PUBCHEM_COMPOUND_CID']
            xyz = mol.write(format="xyz")

            data = {'pubchem_id': pubchem_id,
                    'xyz': xyz}
            for key in keys:
                if key in mol.data:
                    data[key_map[key]] = mol.data[key]

            pubchem_molecules[pubchem_id] = data

        except KeyError:
            skipped += 1

    os.remove(filename)
    return pubchem_molecules, skipped


list_ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov',)
list_ftp.login()

list_ftp.cwd("pubchem/Compound/CURRENT-Full/SDF")
files = [f for f in list_ftp.nlst() if '.sdf.gz' in f]
total_files = len(files)

tqdm_files = tqdm(files, desc="download")

if not os.path.exists(tmp_dir):
    os.mkdir(tmp_dir)

# download files if not already present
for remote_file in tqdm_files:
    file_location = os.path.join(tmp_dir, remote_file)

    if not os.path.exists(file_location):
        ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov',)
        ftp.login()
        ftp.cwd("pubchem/Compound/CURRENT-Full/SDF")
        ftp.retrbinary("RETR " + remote_file, open(file_location, 'wb').write)

# process the downloaded files
pbar = tqdm(total=total_files)
all_molecules = []
total_skipped = []


def task_done(future):
    try:
        result = future.result()
        all_molecules.extend(result[0])
        total_skipped.append(result[1])

    except (ProcessExpired, TimeoutError) as e:
        print("File {} timed-out".format(e.args[0]))
    except Exception as e:
        raise e
    pbar.update()


with pebble.ProcessPool(max_workers=32) as pool:
    for downloaded_file in files:
        file_location = os.path.join(tmp_dir, downloaded_file)
        f = pool.schedule(process_sdf_file, args=(file_location, ),
                          timeout=480)
        f.add_done_callback(task_done)

pbar.close()

total_skipped = sum(total_skipped)

print("total skipped: {}".format(total_skipped))
print("total saved: {}".format(len(all_molecules)))

print("writing to database")

pubchem = db.pubchem
pubchem.insert_many(all_molecules)

print("done")
