#!/bin/bash
import pandas as pd
import json
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib.pyplot as plt

#df_refined=pd.read_csv("/vast/zf2012/02-21-2024_KANO_ESM2_attention/KANO/data/refined_BindingDB_data.csv")
#pro_func=[]
#for pro_id in df_refined['pro_id']:
#    with open (f"/vast/zf2012/03-13-2024_KANO_DeepFRI/DeepFRI/BindingDB_outputs/{pro_id}_MF_pred_scores.json","r") as file:
#        data=json.load(file)
#        pro=data['Y_hat'][0]
#        pro_func.append(pro)
#df_new=pd.DataFrame(pro_func)
#df_new['pro_id']=df_refined['pro_id']
#df_new=df_new.drop_duplicates()
#df_new.to_csv("refined_BindingDB_deepfri_pro_func.csv",index=False))
df_refined=pd.read_csv("refined_BindingDB_deepfri_pro_func.csv")
df_target=pd.read_csv("deepfri_pro_func_tSNE_target_values.csv")
df=pd.concat([df_refined, df_target],axis=0)
df=df.set_index("pro_id")
data_array = df.values
cos_sim_matrix = cosine_similarity(data_array)
cos_sim_matrix_list=list(cos_sim_matrix)
df['cos_sim_matrix_list']=cos_sim_matrix_list
df_GBA_cos_sim_matrix=df[df.index=='P04062']
df_FEN1_cos_sim_matrix=df[df.index=='P39748']
df_VDR_cos_sim_matrix=df[df.index=='P11473']
df_KAT2A_cos_sim_matrix=df[df.index=='Q92830']
# Plot histograms for FEN1 and GBA on the same figure
plt.figure()  # Create a new figure
plt.hist(df_FEN1_cos_sim_matrix['cos_sim_matrix_list'], bins=20, color='grey', label='FEN1', alpha=0.5)  
plt.hist(df_GBA_cos_sim_matrix['cos_sim_matrix_list'], bins=20, color='pink', label='GBA', alpha=0.5)
plt.xlabel('Cosine Similarity')
plt.ylabel('Frequency')
plt.title('Histogram of Cosine Similarity of Protein Functions (FEN1 and GBA)')
plt.legend()  
plt.savefig('Histogram_of_Cosine_Similarity_FEN1_GBA.png',dpi=800)  # Save the figure
plt.show()

plt.figure()  # Create a new figure
# Separate the plot for VDR and KAT2A
plt.hist(df_VDR_cos_sim_matrix['cos_sim_matrix_list'], bins=20, color='green', label='VDR', alpha=0.5)  
plt.hist(df_KAT2A_cos_sim_matrix['cos_sim_matrix_list'], bins=20, color='yellow', label='KAT2A', alpha=0.5)
plt.xlabel('Cosine Similarity')
plt.ylabel('Frequency')
plt.title('Histogram of Cosine Similarity of Protein Functions (VDR and KAT2A)')
plt.legend()  
plt.savefig('Histogram_of_Cosine_Similarity_VDR_KAT2A.png',dpi=800)  # Save the figure
plt.show()


df_MtbDnaK_cos_sim_matrix=df[df.index=='P9WMJ9']
df_HSP90_cos_sim_matrox=df[df.index=='P07900']
df_Ecoli_DnaK_cos_sim_matrix=df[df.index=='C3TRK2']
df_HSP70_cos_sim_matrix=df[df.index=='P0DMV8']

plt.figure()  # Create a new figure
# Separate the plot for VDR and KAT2A
plt.hist(df_MtbDnaK_cos_sim_matrix['cos_sim_matrix_list'], bins=20, color='purple', label='MtbDnaK', alpha=0.5)
plt.hist(df_KAT2A_cos_sim_matrix['cos_sim_matrix_list'], bins=20, color='orange', label='HSP90', alpha=0.5)
plt.xlabel('Cosine Similarity')
plt.ylabel('Frequency')
plt.title('Histogram of Cosine Similarity of Protein Functions (MtbDnaK and HSP90)')
plt.legend()
plt.savefig('Histogram_of_Cosine_Similarity_MtbDnaK_and_HSP90.png',dpi=800)  # Save the figure
plt.show()


plt.figure()  # Create a new figure
# Separate the plot for VDR and KAT2A
plt.hist(df_Ecoli_DnaK_cos_sim_matrix['cos_sim_matrix_list'], bins=20, color='yellow', label='Ecoli_DnaK', alpha=0.5)
plt.hist(df_HSP70_cos_sim_matrix['cos_sim_matrix_list'],bins=20, color='blue', label='Human_HSP70', alpha=0.5)
plt.xlabel('Cosine Similarity')
plt.ylabel('Frequency')
plt.title('Histogram of Cosine Similarity of HSP70')
plt.legend()
plt.savefig('Histogram_of_Cosine_Similarity_HSP70.png',dpi=800)  # Save the figure
plt.show()


