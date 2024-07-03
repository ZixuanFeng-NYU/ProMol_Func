import os
import pandas as pd

# Load the DataFrame and fill NaN values
df_LIT_PCBA_pdb_code = pd.read_csv("LIT-PCBA_pdb_code.csv")
print(df_LIT_PCBA_pdb_code)

# Initialize lists
PDB_id, Sequence, SMILES, Target_Name, Class = [], [], [], [], []

# Get the list of FASTA files
folder_directory = "LIT-PCBA"
target_list = os.listdir(folder_directory)

# Process each target folder
for target in target_list:
    target_folder = os.path.join(folder_directory, target)
    if os.path.isdir(target_folder):
        pdb_id=df_LIT_PCBA_pdb_code.loc[df_LIT_PCBA_pdb_code['Target Name'] == target, 'PDB code'].values[0]
        print(target, pdb_id)
        # Read sequence from FASTA file
        fasta_path = os.path.join(target_folder,f"{pdb_id}.fasta")
        with open(fasta_path, "r") as f:
            lines = f.readlines()
            seq = ''.join(line.strip() for line in lines if not line.startswith(">"))

            # Read SMILES strings from corresponding ligand file
            ligand_file = os.path.join(target_folder, "actives.smi")
            print(ligand_file)
            with open(ligand_file, "r") as file_2:
                smiles_list = file_2.read().strip().splitlines()
                if not smiles_list:
                    print(f"Warning: {ligand_file} is empty.")
                for smi in smiles_list:
                    smi = smi.split("\t")[0]  # Take the first part if the line contains more than one part
                    PDB_id.append(pdb_id)
                    Sequence.append(seq)
                    SMILES.append(smi)
                    Target_Name.append(target)
                    Class.append(1)

            inactives_file = os.path.join(target_folder, "inactives.smi")
            print(inactives_file)
            with open(inactives_file, "r") as file_3:
                smiles_list = file_3.read().strip().splitlines()
                if not smiles_list:
                    print(f"Warning: {inactives_file} is empty.")
                for smi in smiles_list:
                    smi = smi.split("\t")[0]  # Take the first part if the line contains more than one part
                    PDB_id.append(pdb_id)
                    Sequence.append(seq)
                    SMILES.append(smi)
                    Target_Name.append(target)
                    Class.append(0)

# Create a DataFrame
data = pd.DataFrame({
    "smiles": SMILES,
    "Class": Class,
    "pro_id": PDB_id,
    "Sequence": Sequence,
    "Target_Name": Target_Name
})

# Save DataFrame to CSV
output_csv_path = "LIT-PCBA_protein_ligand_data.csv"
data.to_csv(output_csv_path, index=False)

for item in set(data['Target_Name']):
    df_item=data[data['Target_Name']==item]
    df_item.to_csv(f"LIT_PCBA_EF/data_per_target/{item}_protein_ligands.csv",index=False)

