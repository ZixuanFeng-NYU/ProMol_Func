import os
import pandas as pd
from rdkit import Chem

# Paths to ligand and decoy directories
ligand_dir = "DEKOIS2.0/ligands/"
decoy_dir = "DEKOIS2.0/decoys/"

# Process ligand files
ligand_files = os.listdir(ligand_dir)
for ligand_f in ligand_files:
    if not ligand_f.endswith(".sdf"):  # Skip non-SDF files
        continue

    target = ligand_f.split(".")[0]
    sdf_path = os.path.join(ligand_dir, ligand_f)  # Full file path

    # Load SDF molecules
    ligand_supplier = Chem.SDMolSupplier(sdf_path)

    Target, ligand_smiles, Class = [], [], []
    for mol in ligand_supplier:
        if mol is not None:
            smiles = Chem.MolToSmiles(mol)
            ligand_smiles.append(smiles)
            Target.append(target)
            Class.append(1)

    # Save results to CSV
    df = pd.DataFrame({'pro_id': Target, 'smiles': ligand_smiles, 'Class': Class})
    df.to_csv(os.path.join(ligand_dir, f"{target.lower()}_ligands.csv"), index=False)

# Process decoy files
decoy_files = os.listdir(decoy_dir)
for decoy_f in decoy_files:
    if not decoy_f.endswith(".sdf"):  # Skip non-SDF files
        continue

    target = decoy_f.split("_")[0]
    sdf_path = os.path.join(decoy_dir, decoy_f)  # Full file path

    # Load SDF molecules
    decoy_supplier = Chem.SDMolSupplier(sdf_path)

    Target, decoy_smiles, Class = [], [], []
    for mol in decoy_supplier:
        if mol is not None:
            smiles = Chem.MolToSmiles(mol)
            decoy_smiles.append(smiles)
            Target.append(target)
            Class.append(0)

    # Save results to CSV
    df = pd.DataFrame({'pro_id': Target, 'smiles': decoy_smiles, 'Class': Class})
    df.to_csv(os.path.join(decoy_dir, f"{target.lower()}_decoys.csv"), index=False)
# Ensure the output directory exists
os.makedirs("DEKOIS2.0_library_protein", exist_ok=True)
import subprocess
folder_list=os.listdir("DEKOIS2")
for folder in folder_list:
    pdb_file=f"DEKOIS2/{folder}/protein/{folder}_protein.pdb"
    output_fasta=f"DEKOIS2.0_library_protein/{folder}_protein.fasta"
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
# Directory containing the FASTA files
input_dir = "DEKOIS2.0_library_protein"

# Get the list of FASTA files
list_of_fasta = os.listdir(input_dir)

# Iterate through each FASTA file
for file_ in list_of_fasta:
    if file_.endswith(".fasta"):
        fasta_path = os.path.join("DEKOIS2.0_library_protein", file_)

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
combined_data=pd.DataFrame()
fasta_dir=os.listdir("DEKOIS2.0_library_protein")
for file_ in fasta_dir:
    if file_.endswith(".fasta"):
        target=file_.split("_")[0]
        fasta_path = os.path.join("DEKOIS2.0_library_protein", file_)
        with open(fasta_path, "r") as f:
            lines = f.readlines()
            seq = ''.join(line.strip() for line in lines if not line.startswith(">"))
        ligand_file=os.path.join(ligand_dir,f"{target}_ligands.csv")
        df_ligand=pd.read_csv(ligand_file)
        decoy_file=os.path.join(decoy_dir,f"{target}_decoys.csv")
        df_decoy=pd.read_csv(decoy_file)
        df_all=pd.concat([df_ligand,df_decoy])
        df_all['Sequence']=seq
        df_all=df_all[['smiles','Class','pro_id','Sequence']]
        combined_data = pd.concat([combined_data, df_all], ignore_index=True)


# Save DataFrame to CSV
combined_data = combined_data.fillna("NA")
combined_data.to_csv("DEKOIS2.0_protein_ligand_decoys.csv",index=False)
