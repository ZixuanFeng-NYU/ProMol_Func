import os
import pandas as pd



# Get the list of files in the directory
list_of_files = os.listdir("pubchem_bioassays_with_protein_targets")
assay_type_df=pd.read_csv("pubchem_bioassays_with_protein_target.csv")
assay_type_df['AID']=assay_type_df['AID'].astype(int)
# Iterate through each file
for file_name in list_of_files:
    print(file_name)
    aid = file_name.split("_")[1].split(".")[0]  # Extract the assay ID from the file name
    print(aid)
    protein_id = assay_type_df.loc[assay_type_df['AID'] == int(aid), 'Protein_id'].iloc[0]
    print(protein_id)
    # Read the CSV file into a DataFrame
    df = pd.read_csv(f"pubchem_bioassays_with_protein_targets/{file_name}")

    # Select only the relevant columns and rename them
    df = df[['PUBCHEM_EXT_DATASOURCE_SMILES', 'PUBCHEM_ACTIVITY_OUTCOME']]
    df = df.rename(columns={'PUBCHEM_EXT_DATASOURCE_SMILES': 'smiles'})

    # Convert the 'PUBCHEM_ACTIVITY_OUTCOME' column to binary labels (Active: 1, Inactive: 0)
    df['Class'] = df['PUBCHEM_ACTIVITY_OUTCOME'].map({'Active': 1, 'Inactive': 0})

    # Drop rows with missing values
    df = df.dropna()

    # Add the assay ID as a new column
    df['pro_id'] = protein_id
    df=df[['smiles','Class','pro_id']]
    df.to_csv(f"edit_pubchem_bioassay_csv/edited_{file_name}")

