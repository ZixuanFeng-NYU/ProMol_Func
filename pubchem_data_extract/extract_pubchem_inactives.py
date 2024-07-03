import pandas as pd

df_pubchem_inactives = pd.read_csv("/vast/zf2012/05-13-2024_ProMol_Func_general_model/pubchem_inactives.csv")
df_LIT_PCBA=pd.read_csv("testdataset/LIT-PCBA_protein_ligand_data.csv")
df_LIT_PCBA_inactives=df_LIT_PCBA[df_LIT_PCBA['Class']==0.0]

# Find duplicated proteins based on sequence
duplicated_protein = df_pubchem_inactives[df_pubchem_inactives["Sequence"].isin(df_LIT_PCBA_inactives['Sequence'])]
print(duplicated_protein)

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
    lambda row: any(is_identical_smiles(row['smiles'], smiles) for smiles in df_LIT_PCBA_inactive[df_LIT_PCBA_inactive['Sequence'] == row['Sequence']]['SMILES']),
    axis=1
)]

# Print the result
print(duplicated_protein_duplicated_smiles)

# Remove duplicated_protein_duplicated_smiles from df_BindingDB
df_pubchem_inactives_removed_LIT_PCBA = df_pubchem_inactives[~df_pubchem_inactives.index.isin(duplicated_protein_duplicated_smiles.index)]

# Group by 'pro_id'
grouped = df_pubchem_inactives_removed_LIT_PCBA.groupby('pro_id')

# Define the number of samples to extract per 'pro_id'
samples_per_pro_id = 2000  # Adjust as needed

# Calculate the total number of samples
total_samples = samples_per_pro_id * len(grouped)

# Initialize an empty list to store sampled DataFrames
samples = []

# Sample fixed number of rows from each group
for _, group in grouped:
    if len(group) >= samples_per_pro_id:
        sampled_group = group.sample(n=samples_per_pro_id, replace=False)
    else:
        sampled_group = group.sample(n=len(group), replace=True)  # Replace=True to allow sampling with replacement if the group size is smaller
    samples.append(sampled_group)

# Concatenate the sampled subsets
sampled_df = pd.concat(samples)

# Sample remaining rows from the entire DataFrame
remaining_samples = total_samples - len(sampled_df)
if remaining_samples > 0:
    remaining_df = df_pubchem_inactives.drop(sampled_df.index)
    additional_samples = remaining_df.sample(n=remaining_samples, replace=False)
    sampled_df = pd.concat([sampled_df, additional_samples])

# Shuffle the DataFrame to ensure randomness
sampled_df = sampled_df.sample(frac=1).reset_index(drop=True)

# Save the sampled DataFrame to a CSV file
sampled_df.to_csv("/vast/zf2012/05-13-2024_ProMol_Func_general_model/sampled_pubchem_inactives.csv", index=False)

