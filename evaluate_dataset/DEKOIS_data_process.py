import rarfile
with rarfile.RarFile('DEKOIS2.0_library.rar', 'r') as rar:
    rar.extractall()

file_list=rar.namelist()
for file_name in file_list:
    print(file_name)
import gzip
import shutil
import os

extracted_dir = '.'
for file_name in file_list:
    if file_name.endswith('.gz'):
        with gzip.open(os.path.join(extracted_dir, file_name), 'rb') as f_in:
            uncompressed_file_name = os.path.splitext(file_name)[0]
            with open(os.path.join(extracted_dir, uncompressed_file_name), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
print("Files have been successfully unzipped.")

target_list=os.listdir("DEKOIS2.0_library/decoys/")
for target in target_list:
    if target.endswith(".sdf"):
        os.system(f"obabel -isdf DEKOIS2.0_library/decoys/{target} -osmi -O DEKOIS2.0_library/decoys/{target}.smi")


## convert pdb to fasta 
import os
import subprocess

input_dir = "DEKOIS2.0_library_protein"
output_dir = "DEKOIS2.0_library_protein"

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Get the list of protein directories
list_of_protein = os.listdir(input_dir)

# Iterate through each protein directory
for file_pdb in list_of_protein:
    if file_pdb.endswith(".pdb"):
        pdb_id = file_pdb.split(".")[0]
        pdb_file = os.path.join(input_dir, f"{pdb_id}.pdb")
        output_fasta = os.path.join(output_dir, f"{pdb_id}.fasta")

        # Check if the pdb file exists
        if os.path.exists(pdb_file):
            # Construct the command string
            command = f"python pdb2fasta/pdb2fasta.py {pdb_file} > {output_fasta}"

            try:
                # Run the command
                subprocess.run(command, shell=True, check=True)
                print(f"Successfully processed {pdb_file}")
            except subprocess.CalledProcessError as e:
                print(f"Error processing {pdb_file}: {e}")
        else:
            print(f"PDB file does not exist: {pdb_file}")

## save longest chain

# Directory containing the FASTA files
input_dir = "DEKOIS2.0_library_protein"

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
## prepare DEKOIS protein ligand csv

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
