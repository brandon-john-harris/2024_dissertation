import os.path
from rdkit import Chem
from rdkit.Chem import rdMolAlign
import pandas as pd
from itertools import product
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


def filter_files(file_names, search_string, output_folder):
    for file_name in file_names:
        with open(file_name, 'r') as file:
            filtered_lines = [line for line in file if line.startswith(search_string)]
            base_name = os.path.basename(file_name)
            output_file_path = os.path.join(output_folder, base_name)

            with open(output_file_path, 'w') as file:
                file.writelines(filtered_lines)


def match_coords(reference, test, reference_name, chai1_ligand_name, save_output_path):
    pdb1 = Chem.MolFromPDBFile(reference, sanitize=True, removeHs=True)
    pdb2 = Chem.MolFromPDBFile(test, sanitize=False, removeHs=True)
    Chem.MolToMolFile(pdb1, f'{save_output_path}/sanitized_{reference_name}.mol')
    Chem.MolToMolFile(pdb2, f'{save_output_path}/{chai1_ligand_name}.mol')

    mol1 = Chem.MolFromMolFile(f'{save_output_path}/sanitized_{reference_name}.mol', sanitize=True, removeHs=True)
    mol2 = Chem.MolFromMolFile(f'{save_output_path}/{chai1_ligand_name}.mol', sanitize=True, removeHs=True)
    rmsd_match_value = rdMolAlign.CalcRMS(mol1, mol2)
    print('RMSD:', rmsd_match_value)
    return rmsd_match_value

root_path = '/Users/user/Data/'  # path to data, adjust as needed
pdb_code = '6j8i'  # adjust as needed
test_pdb_path = f'{root_path}/{pdb_code}/converted_pdb'

# Create a permutation of possible file names Chai-1 can produce from ranking the models
test_file_permutations = []
for i in range(0, 5):
    for j in range(0, 5):
        test_file_permutations.append(f'{test_pdb_path}/LIG_pred.model_idx_{i}.rank_{j}.pdb')

# Determine which Chai-1 files exist
test_files = []
for file in test_file_permutations:
    if os.path.exists(file):
        test_files.append(file)

# Reference pdb files for RMSD comparison
ref_pdb_path = f'{root_path}/{pdb_code}/reference'  # adjust as needed
reference_files = [f'{ref_pdb_path}/reference_structure.pdb']  # adjust as needed

# Create the DataFrame to save RMSD results
rmsd_results = pd.DataFrame(columns=['PDB Code', 'Reference File', 'Chai-1 File', 'RMSD'])

# Remove everything but the HETATM information in the Chai-1 files, since it can sometimes save as multiple models
save_path = f'{root_path}/{pdb_code}/converted_rdkit_mol'
os.makedirs(save_path, exist_ok=True)
search_string = 'HETATM'
os.makedirs(test_pdb_path, exist_ok=True)
filter_files(test_files, search_string, test_pdb_path)

# Do the same for the reference files
filter_files(reference_files, search_string, ref_pdb_path)

# Create the Reference:Chai-1 file pairs for RMSD calculation
file_pairs = product(reference_files, test_files)

# Calculate best symmetry-corrected RMSD for each Reference:Chai-1 pair (not re-aligning coordinates)
for reference_file, test_file in file_pairs:
    reference_ligand_name = os.path.splitext(os.path.basename(reference_file))[0]
    test_ligand_name = os.path.splitext(os.path.basename(test_file))[0]
    rmsd_value = reordered_coords = match_coords(reference_file, test_file, reference_ligand_name, test_ligand_name,
                                                 save_path)
    rmsd_results = rmsd_results.append(
        {'PDB Code': pdb_code, 'Reference File': f'{reference_ligand_name}', 'Chai-1 File': f'{test_ligand_name}',
         'RMSD': rmsd_value if rmsd_value is not None else 99999}, ignore_index=True)

# Save CSV file containing RMSD information
save_data_path = f'{root_path}/{pdb_code}/processed_data'
os.makedirs(save_data_path, exist_ok=True)
rmsd_results.to_csv(f'{save_data_path}/rmsd_all_combinations.csv', index=False)
exit()
