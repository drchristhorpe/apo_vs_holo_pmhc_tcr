from Bio.PDB import PDBParser, PDBIO, Select
import os
import json

errors = []

with open('data/exclude.json', 'r') as file:
    exclude = json.load(file)


class NotDisordered(Select):
    def accept_atom(self, atom):
        return not atom.is_disordered() or atom.get_altloc() == "A"


def process_file(filename:str, folder:str):
    input_filepath = f"structures/{folder}/raw/{filename}" 
    output_filepath = f"structures/{folder}/clean/{filename}"
    pdb_code = filename.split('_')[0]

    if pdb_code not in exclude:
        parser = PDBParser()
        print (f"Processing {filename}")
        structure = parser.get_structure(pdb_code, input_filepath)

        io = PDBIO()
        io.set_structure(structure)
        io.save(output_filepath, select=NotDisordered())
    else:
        print (f"Skipping {filename}")


for folder in ['apo_pmhc', 'holo_pmhc_tcr']:
    folderpath = f"structures/{folder}/raw"
    for filename in os.listdir(folderpath):
        try:
            process_file(filename, folder)
        except:
            errors.append({'folder':folder,'filemane':filename})

print (errors)

with open('data/altloc_errors.json', 'w') as file:
    json.dump(errors, file, indent=4)



