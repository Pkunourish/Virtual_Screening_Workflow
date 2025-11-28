import os
import pandas as pd
from openbabel import pybel
from freesasa import Structure, calc

# Process: Calculate F104 RSA in holo structure.
# Note that ligand poses are required to generate complex structures before running this code. Not provided.
# Configuration
PROTEIN_STRUC = "./1443.pdb" # apo structure in docking
LIGANDS_DIR = "./ligand" # ligand poses is located at ./ligand/
COMPLEX_DIR = "./complex" # complex structures is stored at ./complex/
CHAIN_ID = 'A' 
RESIDUE_NUM = '1'
output_file = "./rsa_1443_select.csv"

df_map = pd.read_csv('./pharm_1443_select.csv')

def build_complex(protein_file, ligand_file, output_file):
    # Combine protein and ligand into a complex structure.
    # Using ligand poses from vina results. .pdbqt file.
    try:
        protein = next(pybel.readfile("pdb", protein_file))
        ligand = next(pybel.readfile("pdbqt", ligand_file))

        complex_mol = pybel.ob.OBMol()
        complex_mol += protein.OBMol
        complex_mol += ligand.OBMol
        complex_pybel = pybel.Molecule(complex_mol)

        complex_pybel.write("pdb", output_file, overwrite=True)
        return True
    except Exception as e:
        print(f"Error building complex for {ligand_file}: {e}")
        return False

# Reading ligand pose files.    
ligand_files = [f for f in os.listdir(LIGANDS_DIR) if f.endswith(('.pdbqt'))]
print(f"Found {len(ligand_files)} ligands.")

# Calculate RSA for each complex structure using freesasa.
results = []
failed = 0
success = 0
for i, ligand_file in enumerate(ligand_files):
    ligand_path = os.path.join(LIGANDS_DIR, ligand_file)
    ligand_name = os.path.splitext(ligand_file)[0]
    complex_pdb = os.path.join(COMPLEX_DIR, f"{ligand_name}_complex.pdb")

    print(f"processing ligand {i+1}/{len(ligand_files)}: {ligand_name}")
    build_complex(PROTEIN_STRUC, ligand_path, complex_pdb)
    
    try:
        structure = Structure(complex_pdb, options={"hetatm":True})
    except Exception as e:
        print(f"error reading complex structure for {ligand_name}: {e}")
        continue
    
    result = calc(structure)
    if result.residueAreas()[CHAIN_ID][RESIDUE_NUM].relativeTotal is None:
        print(f"failed to calculate RSA for {ligand_name}")
        failed = failed + 1
    else:
        results.append({
        'Ligand': ligand_name,
        'Total_RSA' : result.residueAreas()[CHAIN_ID][RESIDUE_NUM].relativeTotal,
        'Side_chain_RSA' : result.residueAreas()[CHAIN_ID][RESIDUE_NUM].relativeSideChain
        })
        success = success + 1
print(f"RSA calculation completed: {len(ligand_files)} total, {success} successful, {failed} failed.")

results_df = pd.DataFrame(results)
results_df = results_df.sort_values(by = 'Side_chain_RSA', ascending=True)
results_df.columns = ['title_vina', 'Total_RSA', 'Side_chain_RSA']
results_df['title_vina'] = results_df['title_vina'] + ".pdbqt"

# SASA results is added to the map file.
df_map = pd.read_csv('./pharm_1443_select.csv')
rsa_select = df_map.merge(results_df.copy(), 
    on = 'title_vina',
    how='left',
    suffixes=('', '_dup'))
rsa_select = rsa_select[rsa_select['Side_chain_RSA']<=0.07]
rsa_select = rsa_select[rsa_select['Total_RSA']<=0.35]
rsa_select.to_csv(output_file, index=False)

# Optional: extract molecule names for binding pose extractions. dos format.
rsa_select['title_vina'].copy().to_csv('./1443_sasa_vinalist.txt', index=False, header=False)
rsa_select['title_glide'].copy().to_csv('./1443_sasa_glidelist.txt', index=False, header=False)