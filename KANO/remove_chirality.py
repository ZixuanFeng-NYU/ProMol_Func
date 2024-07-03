import pandas as pd
import sys

def remove_chirality(df_path):
    df = pd.read_csv(df_path)
    df_name = df_path.split("/")[-1]  # Corrected index to get the file name
    SMILES = []
    for smiles in df['smiles']:
        smiles = smiles.replace("@", "")
        SMILES.append(smiles)
    df['smiles'] = SMILES
    df.to_csv(f"data_rm_chirality/{df_name}", index=False)  # Added index=False to avoid writing row indices

if __name__ == "__main__":
    df_path = sys.argv[1]
    remove_chirality(df_path)



