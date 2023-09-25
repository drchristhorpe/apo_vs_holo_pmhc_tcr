import json

import functions


with open('data/structures.json', 'r') as file:
    structures = json.load(file)

with open('data/exclude.json', 'r') as file:
    exclude = json.load(file)

peptide_rmsds = {}
abd_rmsds = {}
peptide_errors = {}
abd_errors = {}

both_success_count = 0

for complex in structures:
    if len(structures[complex]['apo']) > 0:
        peptide_rmsds[complex] = {}
        abd_rmsds[complex] = {}

        print (complex)
        peptide = structures[complex]['peptide']
        allele = structures[complex]['allele']
        
        print (allele)
        print (peptide)

        peptide_success = False
        abd_success = False

        for holo_structure in structures[complex]['holo']:
            holo_pdb_code = holo_structure['pdb_code']
            if holo_pdb_code in exclude:
                break
            else:
                for apo_structure in structures[complex]['apo']:
                    apo_pdb_code = apo_structure['pdb_code']
                    if apo_pdb_code in exclude:
                        break
                    else:
                        print (f"Compare {holo_pdb_code} vs {apo_pdb_code}")
                        position_rmsds, success, errors = functions.compare_apo_to_holo(apo_pdb_code, holo_pdb_code, 'peptide', len(peptide))

                        if success:
                            print (position_rmsds)
                            comparison_key = f"{apo_pdb_code}_vs_{holo_pdb_code}"
                            peptide_rmsds[complex][comparison_key] = position_rmsds
                            peptide_success = True
                        else:
                            print (errors)
                            if not complex in peptide_errors:
                                peptide_errors[complex] = []
                            peptide_errors[complex].append(
                                {
                                    'errors': errors,
                                    'apo': apo_pdb_code,
                                    'holo': holo_pdb_code
                                }
                            )
                        
                        position_rmsds, success, errors = functions.compare_apo_to_holo(apo_pdb_code, holo_pdb_code, 'abd', 180)

                        if success:
                            print (position_rmsds)
                            comparison_key = f"{apo_pdb_code}_vs_{holo_pdb_code}"
                            abd_rmsds[complex][comparison_key] = position_rmsds
                            abd_success = True
                        else:
                            print (errors)
                            if not complex in abd_errors:
                                abd_errors[complex] = []
                            abd_errors[complex].append(
                                {
                                    'errors': errors,
                                    'apo': apo_pdb_code,
                                    'holo': holo_pdb_code
                                }
                            )
        print ('')
        if peptide_success and abd_success:
            both_success_count += 1


print (both_success_count)

with open('data/peptide_errors.json', 'w') as file:
    json.dump(peptide_errors, file, indent=4)

with open('data/abd_errors.json', 'w') as file:
    json.dump(abd_errors, file, indent=4)

with open('data/peptide_rmsds.json', 'w') as file:
    json.dump(peptide_rmsds, file, indent=4)

with open('data/abd_rmsds.json', 'w') as file:
    json.dump(abd_rmsds, file, indent=4)

