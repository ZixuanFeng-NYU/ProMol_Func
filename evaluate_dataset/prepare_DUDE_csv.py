import pandas as pd
import os

# Load the DataFrame and fill NaN values
df_DUDE_pdb_code = pd.read_csv("DUDE_targets_pdb_code.csv")
df_DUDE_pdb_code['Target Name'] = df_DUDE_pdb_code['Target Name'].str.lower()
print(df_DUDE_pdb_code)

# Initialize lists
PDB_id, Sequence, SMILES, Target_Name, Class = [], [], [], [],[]

# Get the list of FASTA files
fasta_directory = "DUD-E/DUD-E_targets_pdb/"
file_list = os.listdir(fasta_directory)

# Process each FASTA file
for file_ in file_list:
    if file_.endswith(".fasta"):
        pdb_id = file_.split(".")[0]

        # Retrieve the target name and replace NaN with 'NA'
        target_name = df_DUDE_pdb_code.loc[df_DUDE_pdb_code['PDB code'] == pdb_id, 'Target Name'].values[0]
        print(target_name)
        # Read sequence from FASTA file
        fasta_path = os.path.join(fasta_directory, file_)
        with open(fasta_path, "r") as f:
            lines = f.readlines()
            seq = ''.join(line.strip() for line in lines if not line.startswith(">"))

        # Read SMILES strings from corresponding ligand file
        ligand_file = os.path.join("DUD-E", "all",f"{target_name}", "actives_final.ism")
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
                Target_Name.append(target_name)
                Class.append(1)
        
        decoys_file=os.path.join("DUD-E","all",f"{target_name}","decoys_final.ism")
        print(decoys_file)
        with open(decoys_file,"r") as file_3:
            smiles_list=file_3.read().strip().splitlines()
            if not smiles_list:
                print(f"Warning: {ligand_file} is empty.")
            for smi in smiles_list:
                smi = smi.split("\t")[0]  # Take the first part if the line contains more than one part
                PDB_id.append(pdb_id)
                Sequence.append(seq)
                SMILES.append(smi)
                Target_Name.append(target_name)
                Class.append(0)


# Create a DataFrame
data = pd.DataFrame({
    "smiles": SMILES,
    "Class": Class,
    "pro_id": PDB_id,
    "Sequence": Sequence,
    "Target_Name": Target_Name
})
print("all data:",data)
# Save DataFrame to CSV
output_csv_path = "DUDE_protein_ligand_data.csv"
data.to_csv(output_csv_path, index=False)

for item in set(data['Target_Name']):
    print("item:",item)
    df_item=data[data['Target_Name']==item]
    print("df_item",df_item)
    df_item.to_csv(f"DUDE_EF/data_per_target/{item}_protein_ligands.csv",index=False)


