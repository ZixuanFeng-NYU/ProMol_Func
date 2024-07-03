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

