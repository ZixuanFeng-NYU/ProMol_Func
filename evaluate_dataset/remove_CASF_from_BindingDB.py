import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import TanimotoSimilarity
import numpy as np
# Load the data
df_BindingDB = pd.read_csv("/vast/zf2012/05-16-2024_BindingDB_data/BindingDB_IC50_lower_than_5uM.csv")
df_CASF = pd.read_csv("CASF-2016_protein_ligand_data.csv")

# Find duplicated proteins based on sequence
duplicated_protein = df_BindingDB[df_BindingDB["BindingDB Target Chain Sequence"].isin(df_CASF['Sequence'])]

# Function to calculate similarity between SMILES strings
def is_identical_smiles(smiles1, smiles2):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    if mol1 is None or mol2 is None:
        return False
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2)
    fp1_array = np.array(fp1)
    length_fp1 = len(fp1_array)
    print(length_fp1)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2)
    print(TanimotoSimilarity(fp1, fp2))

    return TanimotoSimilarity(fp1, fp2)==1.0 

# Further filter to find duplicated proteins with duplicated SMILES
duplicated_protein_duplicated_smiles = duplicated_protein[duplicated_protein.apply(
    lambda row: any(is_identical_smiles(row['Ligand SMILES'], casf_smiles) for casf_smiles in df_CASF[df_CASF['Sequence'] == row['BindingDB Target Chain Sequence']]['SMILES']),
    axis=1
)]

# Print the result
print(duplicated_protein_duplicated_smiles)


