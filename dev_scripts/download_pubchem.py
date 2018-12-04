"""
This module implements a script for extracting molecule names and smiles
from the pubchem database.

It uses multiprocessing to download many files simulatenously.

WARNING: While in progress the disk usage maybe fairly large (several GB).
"""

import ftplib
import multiprocessing
import os

import pybel

from tqdm import tqdm
from monty.serialization import dumpfn

key_name_map = {'PUBCHEM_IUPAC_OPENEYE_NAME': 'openeye',
                'PUBCHEM_IUPAC_NAME': 'iupac',
                'PUBCHEM_IUPAC_CAS_NAME': 'cas',
                'PUBCHEM_IUPAC_SYSTEMATIC_NAME': 'iupac',
                'PUBCHEM_IUPAC_TRADITIONAL_NAME': 'traditional'}

n_proc = multiprocessing.cpu_count()


def process_sdf_file(filename):
    skipped = 0
    smiles_iso_to_can = {}
    smiles_can_to_name = {}

    for i, mol in enumerate(pybel.readfile('sdf', filename)):
        try:
            c_smiles_can = mol.data['PUBCHEM_OPENEYE_CAN_SMILES']
            c_smiles_iso = mol.data['PUBCHEM_OPENEYE_ISO_SMILES']

            names = {}
            for key in ['PUBCHEM_IUPAC_CAS_NAME',
                        'PUBCHEM_IUPAC_NAME',
                        'PUBCHEM_IUPAC_TRADITIONAL_NAME']:
                if key in mol.data:
                    names[key_name_map[key]] = mol.data[key]

            if not names:
                raise KeyError("name")

            smiles_iso_to_can[c_smiles_iso] = c_smiles_can
            smiles_can_to_name[c_smiles_can] = names

        except KeyError:
            # skip if no smiles or no names
            skipped += 1

    return smiles_iso_to_can, smiles_can_to_name, skipped


def download_and_process(filename):
    ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov',)
    ftp.login()
    ftp.cwd("pubchem/Compound/CURRENT-Full/SDF")
    ftp.retrbinary("RETR " + filename, open(filename, 'wb').write)

    data = process_sdf_file(filename)
    os.remove(filename)

    return data


list_ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov',)
list_ftp.login()

list_ftp.cwd("pubchem/Compound/CURRENT-Full/SDF")
files = [f for f in list_ftp.nlst() if '.sdf.gz' in f]
total_files = len(files)

pbar = tqdm(total=total_files)
results = []


def update(val):
    results.append(val)
    pbar.update()


pool = multiprocessing.Pool(n_proc)
for file in files:
    pool.apply_async(download_and_process, args=(file, ), callback=update)

pool.close()
pool.join()
pbar.close()

smiles_iso_to_can_full = {}
smiles_can_to_name_full = {}
total_skipped = 0

for smiles_iso, smiles_can, skip in results:
    smiles_iso_to_can_full.update(smiles_iso)
    smiles_can_to_name_full.update(smiles_can)
    total_skipped += skip

print("total skipped: {}".format(total_skipped))
print("total saved: {}".format(len(smiles_can_to_name_full.keys())))

smiles_can_to_name_full['smiles_map'] = smiles_iso_to_can_full

dumpfn(smiles_can_to_name_full, "smiles_to_name.json.gz")
