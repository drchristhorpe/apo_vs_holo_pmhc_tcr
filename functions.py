import os
from biopandas.pdb import PandasPdb
import math


def apo_holo_folder(apo:bool=False):
    if apo:
        return "apo_pmhc"
    else:
        return "holo_pmhc_tcr"


def fetch_pdb_file(pdb_code:str, domain:str, apo:bool=False):
    structure_filename = f"{pdb_code}_1_{domain}.pdb"

    filename = f"structures/{apo_holo_folder(apo)}/raw/{structure_filename}"

    url = f"https://coordinates.histo.fyi/structures/downloads/class_i/without_solvent/{structure_filename}"
    if not os.path.exists(filename):
        print (f"Fetching {pdb_code}...")
        os.system(f"wget -O {filename} {url}")
    else:
        if os.path.getsize(filename) == 0:
            print (f"Fetching {pdb_code}...")
            os.system(f"wget -O {filename} {url}")
        print (f"{pdb_code} already fetched.")


def load_pdb_file(pdb_code:str, domain:str, apo:bool=False):
    filename = f"structures/{apo_holo_folder(apo)}/clean/{pdb_code}_1_{domain}.pdb"
    ppdb = PandasPdb().read_pdb(filename)
    df = ppdb.df['ATOM']
    return df


def compare_apo_to_holo(apo_pdb_code, holo_pdb_code, domain, length, selection='all'):
    apo_df = load_pdb_file(apo_pdb_code, domain, apo=True)
    holo_df = load_pdb_file(holo_pdb_code, domain, apo=False)

    errors = []
    success = False
    position_rmsds = None

    position_rmsds = {}

    for residue_number in range(0, length):
        residue_number += 1
        position_rmsds[residue_number] = {}
        try:
            apo_residue = apo_df[apo_df["residue_number"] == residue_number]
            holo_residue = holo_df[holo_df["residue_number"] == residue_number]

            if apo_residue.shape == holo_residue.shape:
                if selection == 'backbone':
                    residue_rmsd = PandasPdb.rmsd(apo_residue, holo_residue, s='main chain', invert=False)
                elif selection == 'sidechain':
                    residue_rmsd = PandasPdb.rmsd(apo_residue, holo_residue, s='main chain', invert=True)
                else:
                    residue_rmsd = PandasPdb.rmsd(apo_residue, holo_residue)
                position_rmsds[residue_number] = residue_rmsd
                if math.isnan(residue_rmsd):
                   position_rmsds[residue_number] = None 
            else:
                errors.append(f"P{residue_number} cannot be compared due to a dataframe shape mismatch. Apo {apo_pdb_code} vs holo {holo_pdb_code}")
                position_rmsds[residue_number] = None
        except:
            errors.append(f"P{residue_number} cannot be compared due to a missing residue. Apo {apo_pdb_code} vs holo {holo_pdb_code}")
            position_rmsds[residue_number] = None

    if len(errors) > 0:
        success = False
    else:
        success = True
    return position_rmsds, success, errors
