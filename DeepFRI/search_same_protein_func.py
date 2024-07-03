import json
from sklearn.metrics.pairwise import cosine_similarity
import os
import pandas as pd

HSP70_related_protein=[]
pro_func_list=os.listdir("/vast/zf2012/03-13-2024_KANO_DeepFRI/DeepFRI/BindingDB_outputs")
with open("/vast/zf2012/03-13-2024_KANO_DeepFRI/DeepFRI/BindingDB_outputs/C3TRK2_MF_pred_scores.json","r") as file:
    data=json.load(file)
    pro_func_Ecoli=data['Y_hat']
    print(pro_func_Ecoli)
    for pro_id in pro_func_list:
        with open (f"/vast/zf2012/03-13-2024_KANO_DeepFRI/DeepFRI/BindingDB_outputs/{pro_id}","r") as file:
            data=json.load(file)
            pro_func=data['Y_hat']
            similarity_matrix = cosine_similarity(pro_func_Ecoli, pro_func)
            similarity=similarity_matrix[0][0]
            if round(similarity,2) >=0.80:
                print(pro_id, "have similar pro_func with Ecoli_DnaK")
                HSP70_related_protein.append(pro_id.split("_")[0])


pro_func_list=os.listdir("/vast/zf2012/03-13-2024_KANO_DeepFRI/DeepFRI/BindingDB_outputs")
with open("/vast/zf2012/03-13-2024_KANO_DeepFRI/DeepFRI/BindingDB_outputs/P9WMJ9_MF_pred_scores.json","r") as file:
    data=json.load(file)
    pro_func_Mtb=data['Y_hat']
    print(pro_func_Mtb)
    for pro_id in pro_func_list:
        with open (f"/vast/zf2012/03-13-2024_KANO_DeepFRI/DeepFRI/BindingDB_outputs/{pro_id}","r") as file:
            data=json.load(file)
            pro_func=data['Y_hat']
            similarity_matrix = cosine_similarity(pro_func_Mtb, pro_func)
            similarity=similarity_matrix[0][0]
            if round(similarity,2) >=0.80:
                print(pro_id, "have similar pro_func with Mtb_DnaK")
                HSP70_related_protein.append(pro_id.split("_")[0])

pro_func_list=os.listdir("/vast/zf2012/03-13-2024_KANO_DeepFRI/DeepFRI/BindingDB_outputs")
with open("/vast/zf2012/03-13-2024_KANO_DeepFRI/DeepFRI/BindingDB_outputs/P0DMV8_MF_pred_scores.json","r") as file:
    data=json.load(file)
    pro_func_Human_HSP70=data['Y_hat']
    print(pro_func_Human_HSP70)
    for pro_id in pro_func_list:
        with open (f"/vast/zf2012/03-13-2024_KANO_DeepFRI/DeepFRI/BindingDB_outputs/{pro_id}","r") as file:
            data=json.load(file)
            pro_func=data['Y_hat']
            similarity_matrix = cosine_similarity(pro_func_Human_HSP70, pro_func)
            similarity=similarity_matrix[0][0]
            if round(similarity,2) >=0.80:
                print(pro_id, "have similar pro_func with Human HSP70")
                HSP70_related_protein.append(pro_id.split("_")[0])
print(HSP70_related_protein)

df_refined_BindingDB=pd.read_csv("/vast/zf2012/01-31-2024_KANO_ESM2/esm/BindingDB_IC50.csv")
df_related_protein=df_refined_BindingDB[df_refined_BindingDB['UniProt (SwissProt) Primary ID of Target Chain'].isin(HSP70_related_protein)]
df_related_protein.to_csv("Mix_HSP70_related_protein_in_BindingDB.csv",index=False)

print("pro_func cosine similarity between MtbDnaK and Ecoli DnaK",cosine_similarity(pro_func_Ecoli, pro_func_Mtb))
print("pro_func cosine similarity between Human HSP70 and Ecoli DnaK",cosine_similarity(pro_func_Human_HSP70, pro_func_Ecoli))
print("pro_func cosine similarity between Human HSP70 and Mtb DnaK",cosine_similarity(pro_func_Human_HSP70, pro_func_Mtb))

GBA1_related_protein=[]
pro_func_list=os.listdir("/vast/zf2012/03-13-2024_KANO_DeepFRI/DeepFRI/BindingDB_outputs")
with open("/vast/zf2012/03-13-2024_KANO_DeepFRI/DeepFRI/BindingDB_outputs/P04062_MF_pred_scores.json","r") as file:
    data=json.load(file)
    pro_func_Human_GBA1=data['Y_hat']
    print(pro_func_Human_GBA1)
    for pro_id in pro_func_list:
        with open (f"/vast/zf2012/03-13-2024_KANO_DeepFRI/DeepFRI/BindingDB_outputs/{pro_id}","r") as file:
            data=json.load(file)
            pro_func=data['Y_hat']
            similarity_matrix = cosine_similarity(pro_func_Human_GBA1, pro_func)
            similarity=similarity_matrix[0][0]
            if round(similarity,2) >=0.80:
                print(pro_id, "have similar pro_func with Human GBA1")
                GBA1_related_protein.append(pro_id.split("_")[0])
print(GBA1_related_protein)

df_refined_BindingDB=pd.read_csv("/vast/zf2012/01-31-2024_KANO_ESM2/esm/BindingDB_IC50.csv")
df_related_protein=df_refined_BindingDB[df_refined_BindingDB['UniProt (SwissProt) Primary ID of Target Chain'].isin(GBA1_related_protein)]
df_related_protein.to_csv("Human_GBA1_related_protein_in_BindingDB.csv",index=False)

