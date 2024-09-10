## convert pdb to fasta file
import os
import subprocess

input_dir = "LIT-PCBA"
output_dir = "LIT-PCBA"

# Ensure the output directory exists
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Get the list of protein directories
list_of_folders = os.listdir(input_dir)

# Iterate through each protein directory
for target in list_of_folders:
    target_path = os.path.join(input_dir, target)

    if os.path.isdir(target_path):
        # Iterate through each file in the protein directory
        for file in os.listdir(target_path):
            if file.endswith(".pdb"):
                pdb_id = os.path.splitext(file)[0]
                pdb_file = os.path.join(target_path, file)
                output_fasta = os.path.join(output_dir, target, f"{pdb_id}.fasta")

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

## save longest sequence, monomer if multiple chain
import os

# Directory containing the FASTA files
input_dir = "LIT-PCBA"

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

## prepare protein ligand data

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
