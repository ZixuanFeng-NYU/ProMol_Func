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

# Function to calculate EF for each pro_id
def calculate_ef(df, percentages=[0.005, 0.01, 0.05, 0.10]):
    ef_results = []

    # Process each pro_id separately
    for pro_id, group in df.groupby('pro_id'):
        group = group.sort_values(by='preds', ascending=False)  # Sort by prediction scores

        total_compounds = len(group)  # Total compounds for this pro_id
        total_actives = (group['Targets'] == 1).sum()  # Total actives
        total_decoys = (group['Targets'] == 0).sum()  # Total decoys

        ef_dict = {
            'pro_id': pro_id,
            'total_actives': total_actives,
            'total_decoys': total_decoys,
            'total_compounds': total_compounds
        }

        # Calculate EF for each percentage
        for pct in percentages:
            num_selected = max(1, math.ceil(total_compounds * pct))  # Ensure at least 1 compound
            top_selected = group.head(num_selected)  # Get top x% compounds
            actives_in_top = (top_selected['Targets'] == 1).sum()  # Count actives in top x%

            # Compute EF (handle case when total_actives is zero)
            ef = ((actives_in_top / num_selected) / (total_actives / total_compounds)) if total_actives > 0 else 0
            ef_dict[f'EF{int(pct*100)}%'] = round(ef, 4)  # Round for better readability

        ef_results.append(ef_dict)

    ef_df = pd.DataFrame(ef_results)

    # Add a final row with mean EF values
    mean_values = ef_df.mean(numeric_only=True).to_dict()
    mean_values['pro_id'] = 'Mean'  # Label the mean row
    ef_df = pd.concat([ef_df, pd.DataFrame([mean_values])], ignore_index=True)

    return ef_df

# Generalized overlap checking function
def check_overlap(df_query, df_query_preds, label):
    df_overlapped_all=pd.DataFrame()
    df_query_preds['canonical'] = [canonical_smiles(s) for s in df_query_preds['smiles']]
    for seq in df_query['Sequence'].unique():
        overlap = df_general[df_general['Sequence'] == seq].copy()
        if not overlap.empty:
            pro_ids = df_query[df_query['Sequence'] == seq]['pro_id'].unique()
            #print(pro_ids)
            overlap['canonical'] = [canonical_smiles(s) for s in overlap['smiles']]
            df_query_seq = df_query_preds[df_query_preds['pro_id'].isin(pro_ids)].copy()
            df_query_seq['canonical'] = [canonical_smiles(s) for s in df_query_seq['smiles']]
            #print(df_query_seq)
            df_overlap = df_query_seq[df_query_seq['canonical'].isin(overlap['canonical'])]
            if not df_overlap.empty:
                print(f"\n=== Overlap found in {label} for protein: {pro_ids} ===")
                print(df_overlap)
                df_overlapped_all = pd.concat([df_overlapped_all, df_overlap], ignore_index=True)
                print(len(df_overlap), "compounds overlapped")
    if not df_overlapped_all.empty:
        df_query_preds_cleaned = df_query_preds[~df_query_preds[['canonical', 'pro_id']].apply(tuple, axis=1).isin(df_overlapped_all[['canonical', 'pro_id']].apply(tuple, axis=1))]
        print(df_query_preds_cleaned)
        print(df_query_preds)
        ef_df=calculate_ef(df_query_preds_cleaned)
        print(ef_df)



# Run overlap checks
check_overlap(df_lit, df_lit_prediction, 'LIT-PCBA')
check_overlap(df_dekois, df_dekois_prediction, 'DEKOIS')
check_overlap(df_dude, df_dude_prediction, 'DUDE')
check_overlap(df_casf2016, df_casf2016_prediction, 'CASF2016')

