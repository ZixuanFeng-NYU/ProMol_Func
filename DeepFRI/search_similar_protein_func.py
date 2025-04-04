import json
import os
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import sys

# Load all protein functions
def load_protein_functions(pro_func_list, directory):
    pro_func_dict = {}
    for pro_id in pro_func_list:
        # Ensure the filename does not already include the suffix
        if not pro_id.endswith("_MF_pred_scores.json"):
            pro_id = pro_id + "_MF_pred_scores.json"
        file_path = os.path.join(directory, pro_id)
        try:
            with open(file_path, "r") as file:
                data = json.load(file)
                pro_func_dict[pro_id] = data['Y_hat']
        except Exception as e:
            print("Error loading {}: {}".format(file_path, e))
    return pro_func_dict
def search_similar_protein_functions(target_pro_id, directory):
    similari_protein_list,similarity_list = [],[]
    pro_func_list = os.listdir(directory)

    # Load all protein functions
    pro_func_dict = load_protein_functions(pro_func_list, directory)

    # Ensure the target filename does not already include the suffix
    if not target_pro_id.endswith("_MF_pred_scores.json"):
        target_pro_id = target_pro_id + "_MF_pred_scores.json"
    target_file_path = os.path.join(directory, target_pro_id)
    try:
        with open(target_file_path, "r") as file:
            target_data = json.load(file)
            pro_func_target = target_data['Y_hat']
    except Exception as e:
        print("Error loading {}: {}".format(target_file_path, e))
        return

    # Calculate similarity
    for pro_id, pro_func in pro_func_dict.items():
        if pro_id != target_pro_id:
            similarity_matrix = cosine_similarity(pro_func_target, pro_func)
            similarity = similarity_matrix[0][0]
            if round(similarity, 2) >= 0.80:
                print("{} has similar pro_func with {}".format(pro_id, target_pro_id))
                similari_protein_list.append(pro_id)
                similarity_list.append(similarity)

    # Save similar proteins to CSV
    print("similari_protein_list",similari_protein_list)
    df = pd.DataFrame(similari_protein_list, columns=['Protein_ID'])
    df['similarity']=similarity_list
    df.to_csv(f"{target_pro_id}_similar_protein.csv", index=False)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python script.py <target_protein_id>")
        sys.exit(1)

    target_pro_id = sys.argv[1]
    directory = "DeepFRI_outputs"
    search_similar_protein_functions(target_pro_id, directory)
    df_similar=pd.read_csv(f"{target_pro_id}_MF_pred_scores.json_similar_protein.csv")
    simlar_protein_funcs=[]
    for id_ in df_similar['Protein_ID']:
        simlar_protein_funcs.append(id_.split("_")[0])
    print("simlar_protein_funcs id:",simlar_protein_funcs)
    df_general=pd.read_csv("../KANO/data/ProMol_Func_general_model_data_06102024_add_decoys_0331version.csv")
    df_general['pro_id']=df_general['pro_id'].astype(str)
    df_similar_found=df_general[df_general['pro_id'].isin(simlar_protein_funcs)]
    print(df_similar_found)
    df_similar_found.to_csv(f"{target_pro_id}_similar_protein_data.csv", index=False)


