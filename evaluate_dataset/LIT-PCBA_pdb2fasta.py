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

