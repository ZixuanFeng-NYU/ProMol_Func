import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import TanimotoSimilarity

# Load the data
df_BindingDB = pd.read_csv("/vast/zf2012/05-16-2024_BindingDB_data/BindingDB_IC50_lower_than_5uM.csv")
df_DUDE = pd.read_csv("DUDE_protein_ligand_data.csv")
print(df_BindingDB)
# Find duplicated proteins based on sequence
duplicated_protein = df_BindingDB[df_BindingDB["BindingDB Target Chain Sequence"].isin(df_DUDE['Sequence'])]
print(duplicated_protein)
duplicated_protein_2 = df_DUDE[df_DUDE["Sequence"].isin(df_BindingDB["BindingDB Target Chain Sequence"])]
print(duplicated_protein_2)
# Function to calculate similarity between SMILES strings
def is_identical_smiles(smiles1, smiles2):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    if mol1 is None or mol2 is None:
        return False
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2)
    #print(TanimotoSimilarity(fp1, fp2))

    return TanimotoSimilarity(fp1, fp2)==1.0 

# Further filter to find duplicated proteins with duplicated SMILES
duplicated_protein_duplicated_smiles = duplicated_protein[duplicated_protein.apply(
    lambda row: any(is_identical_smiles(row['Ligand SMILES'], DUDE_smiles) for DUDE_smiles in df_DUDE[df_DUDE['Sequence'] == row['BindingDB Target Chain Sequence']]['SMILES']),
    axis=1
)]

# Print the result
print(duplicated_protein_duplicated_smiles)
print(duplicated_protein_duplicated_smiles.columns)
print(duplicated_protein_duplicated_smiles["UniProt (SwissProt) Primary ID of Target Chain"])
# Remove duplicated_protein_duplicated_smiles from df_BindingDB
df_BindingDB_removed_DUDE = df_BindingDB[~df_BindingDB.index.isin(duplicated_protein_duplicated_smiles.index)]

# Print the filtered DataFrame
df_BindingDB_removed_DUDE.to_csv("BindingDB_removed_DUDE.csv",index=False)

