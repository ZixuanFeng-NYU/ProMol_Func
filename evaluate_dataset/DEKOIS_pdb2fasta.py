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

