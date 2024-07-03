import pandas as pd
import os

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

