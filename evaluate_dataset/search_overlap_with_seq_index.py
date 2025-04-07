import pandas as pd
from rdkit import Chem

# Load data
df_lit = pd.read_csv("lit-pcba-protein.csv")
df_lit_prediction = pd.read_csv("../evaluate_dataset/CASF-2016_organized_data_ProMol_Func_prediction.csv")

df_dekois = pd.read_csv("DEKOIS_81targets.csv")
df_dekois_prediction = pd.read_csv("../evaluate_dataset/DEKOIS2_ProMol_Func_prediction.csv")

df_dude = pd.read_csv("DUDE_102targets.csv")
df_dude_prediction = pd.read_csv("../evaluate_dataset/DUDE_ProMol_Func_prediction.csv")

df_casf2016 = pd.read_csv("CASF2016_57_Seq.csv")
df_casf2016_prediction = pd.read_csv("../evaluate_dataset/CASF-2016_organized_data_ProMol_Func_prediction.csv")

df_general = pd.read_csv("../KANO/data/ProMol_Func_general_model_data_06102024_add_decoys_0331version.csv")

# Canonical SMILES function
def canonical_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Chem.MolToSmiles(mol, canonical=True)
    else:
        return None

# Generalized overlap checking function
def check_overlap(df_query, df_query_preds, label):
    for seq in df_query['Sequence'].unique():
        overlap = df_general[df_general['Sequence'] == seq].copy()
        if not overlap.empty:
            pro_ids = df_query[df_query['Sequence'] == seq]['pro_id'].unique()
            overlap['canonical'] = [canonical_smiles(s) for s in overlap['smiles']]
            df_query_seq = df_query_preds[df_query_preds['pro_id'].isin(pro_ids)].copy()
            df_query_seq['canonical'] = [canonical_smiles(s) for s in df_query_seq['smiles']]
            df_overlap = df_query_seq[df_query_seq['canonical'].isin(overlap['canonical'])]
            if not df_overlap.empty:
                print(f"\n=== Overlap found in {label} for sequence: {seq} ===")
                print(df_overlap)

# Run overlap checks
check_overlap(df_lit, df_lit_prediction, 'LIT-PCBA')
check_overlap(df_dekois, df_dekois_prediction, 'DEKOIS')
check_overlap(df_dude, df_dude_prediction, 'DUDE')
check_overlap(df_casf2016, df_casf2016_prediction, 'CASF2016')

