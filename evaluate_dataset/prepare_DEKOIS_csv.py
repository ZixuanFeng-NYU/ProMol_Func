import pandas as pd
import os

# Load the DataFrame and fill NaN values
df_pdb_code = pd.read_csv("DEKOIS2.0_library_targets_pdb_code.csv")
df_pdb_code['Target Name'] = df_pdb_code['Target Name'].fillna('NA')

# Initialize lists
PDB_id, Sequence, SMILES,Class, Target_Name = [], [], [], [],[]

# Get the list of FASTA files
fasta_directory = "DEKOIS2.0_library_protein"
file_list = os.listdir(fasta_directory)

# Process each FASTA file
for file_ in file_list:
    if file_.endswith(".fasta"):
        pdb_id = file_.split(".")[0]

        # Retrieve the target name and replace NaN with 'NA'
        target_name = df_pdb_code.loc[df_pdb_code['PDB code'] == pdb_id, 'Target Name'].values
        target_name = 'NA' if len(target_name) == 0 or pd.isna(target_name[0]) else target_name[0]

        # Read sequence from FASTA file
        fasta_path = os.path.join(fasta_directory, file_)
        with open(fasta_path, "r") as f:
            lines = f.readlines()
            seq = ''.join(line.strip() for line in lines if not line.startswith(">"))

        # Read SMILES strings from corresponding ligand file
        ligand_file = os.path.join("DEKOIS2.0_library", "ligands", f"{target_name}_ligand.smi")
        if os.path.exists(ligand_file):
            with open(ligand_file, "r") as file_2:
                smiles_list = file_2.read().strip().splitlines()
                if not smiles_list:
                    print(f"Warning: {ligand_file} is empty.")
                for smi in smiles_list:
                    smi = smi.split("\t")[0]  # Take the first part if the line contains more than one part
                    PDB_id.append(pdb_id)
                    Sequence.append(seq)
                    SMILES.append(smi)
                    Class.append(1)
                    Target_Name.append(target_name)
        else:
            # If ligand file does not exist, add one row with None for SMILES
            PDB_id.append(pdb_id)
            Sequence.append(seq)
            SMILES.append(None)
            Class.append(None)
            Target_Name.append(target_name)

        # Read SMILES strings from corresponding ligand file
        decoy_file = os.path.join("DEKOIS2.0_library", "decoys", f"{target_name}_Celling-v1.12_decoyset.sdf.smi")
        if os.path.exists(decoy_file):
            with open(decoy_file, "r") as file_2:
                smiles_list = file_2.read().strip().splitlines()
                if not smiles_list:
                    print(f"Warning: {ligand_file} is empty.")
                for smi in smiles_list:
                    smi = smi.split("\t")[0]  # Take the first part if the line contains more than one part
                    PDB_id.append(pdb_id)
                    Sequence.append(seq)
                    SMILES.append(smi)
                    Class.append(0)
                    Target_Name.append(target_name)
        else:
            # If ligand file does not exist, add one row with None for SMILES
            PDB_id.append(pdb_id)
            Sequence.append(seq)
            SMILES.append(None)
            Class.append(None)
            Target_Name.append(target_name)

# Create a DataFrame
data = pd.DataFrame({
    "PDB_id": PDB_id,
    "Sequence": Sequence,
    "SMILES": SMILES,
    "Class":Class,
    "Target_Name": Target_Name
})

# Save DataFrame to CSV
output_csv_path = "DEKOIS2.0_library_protein_ligand_data.csv"
data.to_csv(output_csv_path, index=False)

