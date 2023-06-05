"""This module implements a script for extracting molecule names, SMILES and xyz
files from the Pubchem database.

It uses multiprocessing to process many files simultaneously.

WARNING: While in progress the disk usage maybe fairly large (~100 GB).

Writes to MongograntStore. The data will be inserted into the "collection"
collection in the "mp_pubchem" database in the mongo database. Contact
Shyam D for access to this collection. Alternatively, the script can easily be
reconfigured using alternative maggma stores. E.g. to write to a local MongoDB
database.
"""
from __future__ import annotations

import ftplib
import glob
import os
from concurrent.futures import TimeoutError

import pebble
import pybel
from maggma.advanced_stores import MongograntStore
from pebble import ProcessExpired
from tqdm import tqdm

tmp_dir = "tmp"

keys = [
    "PUBCHEM_OPENEYE_CAN_SMILES",
    "PUBCHEM_OPENEYE_ISO_SMILES",
    "PUBCHEM_IUPAC_CAS_NAME",
    "PUBCHEM_IUPAC_NAME",
    "PUBCHEM_IUPAC_TRADITIONAL_NAME",
    "PUBCHEM_IUPAC_SYSTEMATIC_NAME",
    "PUBCHEM_IUPAC_OPENEYE_NAME",
]

key_map = {
    "PUBCHEM_OPENEYE_CAN_SMILES": "smiles_can",
    "PUBCHEM_OPENEYE_ISO_SMILES": "smiles_iso",
    "PUBCHEM_IUPAC_CAS_NAME": "name_cas",
    "PUBCHEM_IUPAC_NAME": "name_iupac",
    "PUBCHEM_IUPAC_TRADITIONAL_NAME": "name_traditional",
    "PUBCHEM_IUPAC_SYSTEMATIC_NAME": "name_systematic",
    "PUBCHEM_IUPAC_OPENEYE_NAME": "name_openeye",
}


def process_sdf_file(filename):
    mp_pubchem = MongograntStore(
        "rw:knowhere.lbl.gov/mp_pubchem", "mp_pubchem", key="pubchem_id"
    )
    mp_pubchem.connect()
    coll = mp_pubchem.collection

    skipped = 0
    pubchem_molecules = []
    for _i, mol in enumerate(pybel.readfile("sdf", filename)):
        try:
            pubchem_id = int(mol.data["PUBCHEM_COMPOUND_CID"])
            xyz = mol.write(format="xyz")

            data = {"pubchem_id": pubchem_id, "xyz": xyz}
            for key in keys:
                if key in mol.data:
                    data[key_map[key]] = mol.data[key]

            pubchem_molecules.append(data)

        except KeyError:
            skipped += 1

    coll.insert_many(pubchem_molecules)

    os.rename(filename, filename + ".processed")
    return len(pubchem_molecules), skipped


def task_done(future):
    try:
        result = future.result()
        total_completed.append(result[0])
        total_skipped.append(result[1])

    except (ProcessExpired, TimeoutError) as e:
        print(f"File {e.args[0]} timed-out")
    except Exception as e:
        raise e

    pbar.update()


def download_done(future):
    future.result()
    download_pbar.update()


def download_file(remote_file):
    """First downloads the file to filename.tmp then moves to filename after."""
    ftp = ftplib.FTP(
        "ftp.ncbi.nlm.nih.gov",
    )
    ftp.login()

    ftp.cwd("pubchem/Compound/CURRENT-Full/SDF")
    file_location = os.path.join(tmp_dir, remote_file)
    ftp.retrbinary(
        "RETR " + remote_file, open(file_location + ".tmp", "wb").write  # noqa: SIM115
    )

    os.rename(file_location + ".tmp", file_location)
    return True


ftp = ftplib.FTP(
    "ftp.ncbi.nlm.nih.gov",
)
ftp.login()

ftp.cwd("pubchem/Compound/CURRENT-Full/SDF")
files = [f for f in ftp.nlst() if ".sdf.gz" in f]
total_files = len(files)

download_pbar = tqdm(total=total_files, desc="download")

if not os.path.exists(tmp_dir):
    os.mkdir(tmp_dir)

# download files if not already present
with pebble.ProcessPool(max_workers=10) as pool:
    for remote_file in files:
        if not os.path.exists(
            os.path.join(tmp_dir, remote_file)
        ) and not os.path.exists(os.path.join(tmp_dir, remote_file + ".processed")):
            f = pool.schedule(download_file, args=(remote_file,))
            f.add_done_callback(download_done)
        else:
            download_pbar.update()

download_pbar.close()

# process downloaded files
files = glob.glob(os.path.join(tmp_dir, "*.sdf.gz"))
total_files = len(files)

# process the downloaded files
pbar = tqdm(total=total_files, desc="processing")
total_completed = []
total_skipped = []

with pebble.ProcessPool() as pool:
    for downloaded_file in files:
        f = pool.schedule(process_sdf_file, args=(downloaded_file,), timeout=480)
        f.add_done_callback(task_done)

pbar.close()

total_skipped = sum(total_skipped)
total_completed = sum(total_completed)

print(f"total skipped: {total_skipped}")
print(f"total saved: {total_completed}")
