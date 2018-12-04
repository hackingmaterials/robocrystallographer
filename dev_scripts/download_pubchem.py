"""
This module implements a script for extracting molecule names and smiles
from the pubchem database.

It uses multiprocessing to download many files simulatenously.

WARNING: While in progress the disk usage maybe fairly large (several GB).
"""

import ftplib
import pebble
import os

import pybel
from pebble import ProcessExpired

from tqdm import tqdm
from monty.serialization import dumpfn
from concurrent.futures import TimeoutError

keys = ['PUBCHEM_OPENEYE_CAN_SMILES', 'PUBCHEM_OPENEYE_ISO_SMILES',
        'PUBCHEM_IUPAC_CAS_NAME', 'PUBCHEM_IUPAC_NAME',
        'PUBCHEM_IUPAC_TRADITIONAL_NAME', 'PUBCHEM_IUPAC_SYSTEMATIC_NAME',
        'PUBCHEM_IUPAC_OPENEYE_NAME']

list_ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov',)
list_ftp.login()

list_ftp.cwd("pubchem/Compound/CURRENT-Full/SDF")
files = [f for f in list_ftp.nlst() if '.sdf.gz' in f]
total_files = len(files)

pbar = tqdm(total=total_files)
results = []


def process_sdf_file(filename):
    ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov',)
    ftp.login()
    ftp.cwd("pubchem/Compound/CURRENT-Full/SDF")
    ftp.retrbinary("RETR " + filename, open(filename, 'wb').write)

    skipped = 0
    pubchem_molecules = {}
    for i, mol in enumerate(pybel.readfile('sdf', filename)):
        try:
            pubchem_id = mol.data['PUBCHEM_COMPOUND_CID']
            xyz = mol.write(format="xyz")

            data = {'xyz': xyz}
            for key in keys:
                if key in mol.data:
                    data[key] = mol.data[key]

            pubchem_molecules[pubchem_id] = data

        except KeyError:
            skipped += 1

    os.remove(filename)
    return pubchem_molecules, skipped


def task_done(future):
    try:
        result = future.result()
        results.append(result)
    except (ProcessExpired, TimeoutError) as e:
        print("File {} timed-out".format(e.args[1]))
    except Exception as e:
        raise e
    pbar.update()


with pebble.ProcessPool(max_workers=32) as pool:
    for remote_file_path in files:
        f = pool.schedule(process_sdf_file, args=(remote_file_path, ),
                          timeout=60)
        f.add_done_callback(task_done)

pbar.close()

molecules = {}
total_skipped = 0

for r_pubchem_molecules, r_skipped in results:
    molecules.update(r_pubchem_molecules)
    total_skipped += r_skipped

print("total skipped: {}".format(total_skipped))
print("total saved: {}".format(len(molecules.keys())))

dumpfn(molecules, "pubchem.json.gz")
