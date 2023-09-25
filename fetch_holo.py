import requests
import json

import constants
import functions

holo_structures = {}

holo_api_url = f"{constants.BASE_API_URL}/complex_types/class_i_with_peptide_and_alpha_beta_tcr"

r = requests.get(holo_api_url)
data = r.json()['set']  
page_numbers = data['pagination']['pages']

holo_structure_count = 0
complex_count = 0


for page_number in page_numbers:
    
    holo_page_api_url = f"{holo_api_url}?page_number={page_number}"
    r = requests.get(holo_page_api_url)
    members = r.json()['set']['members']

    for member in members:

        holo_structure_count += 1
        pdb_code = member['pdb_code']
        allele = member['allele']['alpha']['slug']
        peptide = member['assigned_chains']['peptide']['sequence']

        print (pdb_code)

        complex_key = f"{allele}_{peptide.lower()}"
        if not complex_key in holo_structures:
            complex_count += 1
            holo_structures[complex_key] = {
            'allele': allele,
            'peptide': peptide,
            'holo': [],
            'apo': []
        }
        holo_structures[complex_key]['holo'].append({'pdb_code':pdb_code, 'resolution': member['resolution']})
        
        # Download the PDB file
        functions.fetch_pdb_file(pdb_code, 'peptide', apo=False)
        functions.fetch_pdb_file(pdb_code, 'abd', apo=False)


with open('data/structures.json', 'w') as file:
    json.dump(holo_structures, file, indent=4)

print (f"{complex_count} complexes have holo structures available.")
print (f"{holo_structure_count} holo structures downloaded.")

stats = {
    'holo': {
        'complexes': complex_count,
        'structures': holo_structure_count
    }
}

with open('data/stats.json', 'w') as file:
    json.dump(stats, file, indent=4)








