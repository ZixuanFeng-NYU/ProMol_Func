## convert small molecule sdf to smiles
import os
import pandas as pd

target_list=os.listdir("DUD-E/all")
for target in target_list:
    files=os.listdir("DUD-E/all/"+target)
    for file_ in files:
        if file_.endswith(".gz"):
            os.system(f"gzip -d DUD-E/all/{target}/{file_}")

from rdkit import Chem
import requests

target_list = os.listdir("DUD-E/all")

import subprocess
input_dir = "DUD-E/DUD-E_targets_pdb"

# Ensure the output directory exists
if not os.path.exists(input_dir):
    os.makedirs(input_dir)

df_DUDE=pd.read_csv("DUDE_targets_pdb_code.csv")

def download_pdb(pdb_id, save_dir="."):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)

    if response.status_code == 200:
        file_path = f"{save_dir}/{pdb_id}.pdb"
        with open(file_path, "w") as f:
            f.write(response.text)
        print(f"PDB file saved: {file_path}")
    else:
        print(f"Error: PDB ID {pdb_id} not found.")
for id_ in df_DUDE['PDB code']:
    download_pdb(id_,input_dir)

# Get the list of protein directories
list_of_protein = os.listdir(input_dir)

# Iterate through each protein directory
for target in list_of_protein:
    pdb_id=target.split(".")[0]
    pdb_file = os.path.join(input_dir, pdb_id+".pdb")
    output_fasta = os.path.join(input_dir, pdb_id+".fasta")

    # Check if the PDB file exists
    if os.path.exists(pdb_file):
        try:
            # Construct the command string
            command = f"python pdb2fasta/pdb2fasta.py {pdb_file} > {output_fasta}"

            # Run the command
            subprocess.run(command, shell=True, check=True)
            print(f"Successfully converted {pdb_file} to {output_fasta}")
        except subprocess.CalledProcessError as e:
            print(f"Error occurred while processing {pdb_file}: {e}")
    else:
        print(f"PDB file not found: {pdb_file}")
## save longest chain

# Directory containing the FASTA files
input_dir = "DUD-E/DUD-E_targets_pdb/"

# Get the list of FASTA files
list_of_fasta = os.listdir(input_dir)

# Iterate through each FASTA file
for file_ in list_of_fasta:
    if file_.endswith(".fasta"):
        fasta_path = os.path.join(input_dir, file_)

        with open(fasta_path, "r") as f:
            # Read all lines in the file
            lines = f.readlines()

        # Variables to track the current sequence and the longest sequence
        current_header = None
        current_sequence = []
        longest_header = None
        longest_sequence = ""
        for line in lines:
            if line.startswith('>'):
                # Save the previous sequence if it was longer than the current longest
                if current_header and len(''.join(current_sequence)) > len(longest_sequence):
                    longest_header = current_header
                    longest_sequence = ''.join(current_sequence)

                # Start a new sequence
                current_header = line.strip()
                current_sequence = []
            else:
                # Add line to the current sequence
                current_sequence.append(line.strip())

        # Final check for the last sequence in the file
        if current_header and len(''.join(current_sequence)) > len(longest_sequence):
            longest_header = current_header
            longest_sequence = ''.join(current_sequence)

        # Write the longest sequence back to the file
        with open(fasta_path, "w") as f_out:
            if longest_header and longest_sequence:
                f_out.write(f"{longest_header}\n")
                f_out.write(f"{longest_sequence}\n")

        # Print the result
        print(f"Processed {file_}: longest sequence length {len(longest_sequence)}")
## prepare protein ligand csv

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
