## convert pdb to fasta file
import os
import subprocess

input_dir = "CASF-2016/coreset"
output_dir = "CASF-2016_protein_fasta"

# Ensure the output directory exists
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Get the list of protein directories
list_of_protein = os.listdir(input_dir)

# Iterate through each protein directory
for pdb_id in list_of_protein:
    pdb_file = f"{input_dir}/{pdb_id}/{pdb_id}_protein.pdb"
    output_fasta = f"{output_dir}/{pdb_id}_protein.fasta"

    # Construct the command string
    command = f"python pdb2fasta/pdb2fasta.py {pdb_file} > {output_fasta}"

    # Run the command
    subprocess.run(command, shell=True, check=True)


## save longest sequence if multiple chain

# Directory containing the FASTA files
input_dir = "CASF-2016_protein_fasta"

# Get the list of FASTA files
list_of_fasta = os.listdir(input_dir)

# Iterate through each FASTA file
for fasta in list_of_fasta:
    fasta_path = os.path.join(input_dir, fasta)

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
    print(f"Processed {fasta}: longest sequence length {len(longest_sequence)}")


## generate protein_ligand data

# Initialize lists
PDB_id, Sequence, SMILES = [], [], []

# Get list of FASTA files
fasta_list = os.listdir("CASF-2016_protein_fasta")

# Process each FASTA file
for fasta in fasta_list:
    pdb_id = fasta.split("_")[0]
    PDB_id.append(pdb_id)

    # Read sequence from FASTA file
    with open(f"CASF-2016_protein_fasta/{fasta}", "r") as f:
        lines = f.readlines()
        seq = ''.join(line.strip() for line in lines if not line.startswith(">"))
        Sequence.append(seq)

    # Read SMILES string from corresponding ligand file
    ligand_file = f"CASF-2016/coreset/{pdb_id}/{pdb_id}_ligand.smi"
    if os.path.exists(ligand_file):
        with open(ligand_file, "r") as file_2:
            smi = file_2.read().strip().split()[0]  # Read and split the line, then take the first part
            SMILES.append(smi)
    else:
        SMILES.append(None)

# Create a DataFrame
data = pd.DataFrame({
    "PDB_id": PDB_id,
    "Sequence": Sequence,
    "SMILES": SMILES
})

# Save DataFrame to CSV
data.to_csv("CASF-2016_protein_ligand_data.csv", index=False)
