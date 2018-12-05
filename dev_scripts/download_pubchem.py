"""
This module implements a script for extracting molecule names, SMILES and xyz
files from the Pubchem database.

It uses multiprocessing to process many files simultaneously.

WARNING: While in progress the disk usage maybe fairly large (100's GB).

Writes to MongograntStore. The data will be inserted into the "collection"
collection in the "mp_pubchem" database in the mongo database.
"""

import ftplib
import pebble
import os

import pybel
from pebble import ProcessExpired

from tqdm import tqdm
from concurrent.futures import TimeoutError
from maggma.advanced_stores import MongograntStore

tmp_dir = "tmp"

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

    coll.insert_many(pubchem_molecules)

    os.remove(filename)
    return len(pubchem_molecules), skipped


def task_done(future):
    try:
        result = future.result()
        total_completed.append(result[0])
        total_skipped.append(result[1])

    except (ProcessExpired, TimeoutError) as e:
        print("File {} timed-out".format(e.args[0]))
    except Exception as e:
        raise e

    pbar.update()


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

mp_pubchem = MongograntStore("rw:knowhere.lbl.gov", "mp_pubchem",
                             key="pubchem_id")
mp_pubchem.connect()
coll = mp_pubchem.collection

# process the downloaded files
pbar = tqdm(total=total_files)
total_completed = []
total_skipped = []

with pebble.ProcessPool(max_workers=32) as pool:
    for downloaded_file in files:
        file_location = os.path.join(tmp_dir, downloaded_file)
        f = pool.schedule(process_sdf_file, args=(file_location, ),
                          timeout=480)
        f.add_done_callback(task_done)

pbar.close()

total_skipped = sum(total_skipped)
total_completed = sum(total_completed)

print("total skipped: {}".format(total_skipped))
print("total saved: {}".format(total_completed))
