import requests
import json

import constants
import functions


apo_and_holo = []

apo_structure_count = 0

with open('data/structures.json', 'r') as file:
    structures = json.load(file)

with open('data/stats.json', 'r') as file:
    stats = json.load(file)

for structure in structures:
    peptide = structures[structure]['peptide']

    if len(peptide) >= 8:
        peptide_api_url = f"{constants.BASE_API_URL}/peptide_sequences/{peptide.lower()}"

        r = requests.get(peptide_api_url)

        members = r.json()['set']['members']

        for member in members:
            if member['complex_type'] == 'class_i_with_peptide':
                structures[structure]['apo'].append({'pdb_code':member['pdb_code'], 'resolution':member['resolution']})

                print (member['pdb_code'])
                # Download the PDB file
                functions.fetch_pdb_file(member['pdb_code'], 'peptide', apo=True)
                functions.fetch_pdb_file(member['pdb_code'], 'abd', apo=True)

                apo_structure_count += 1


            if structure not in apo_and_holo:
                apo_and_holo.append(structure)




with open('data/structures.json', 'w') as file:
    json.dump(structures, file, indent=4)

print (f"{len(apo_and_holo)} complexes have apo and holo structures available.")

print (f"{apo_structure_count} apo structures downloaded.")

stats['apo'] = {
    'complexes': len(apo_and_holo),
    'structures': apo_structure_count
}

with open('data/stats.json', 'w') as file:
    json.dump(stats, file, indent=4)