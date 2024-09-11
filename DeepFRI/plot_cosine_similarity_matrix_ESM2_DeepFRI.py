#!/bin/bash
import pandas as pd
import json
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
import numpy as np

df=pd.read_csv("ESM2_DeepFRI_uniprot_id.csv")
print(df)
import ast
ESM_list=[]
for esm_item in df['ESM2']:
    esm_list=ast.literal_eval(esm_item)
    ESM_list.append(esm_list)
df_ESM2=pd.DataFrame(ESM_list)
ESM2_array=df_ESM2.values
ESM2_cos_sim_matrix = cosine_similarity(ESM2_array)


DeepFRI_list=[]

for pro_func_item in df['DeepFRI']:
    pro_func_list=ast.literal_eval(pro_func_item)
    DeepFRI_list.append(pro_func_list)

df_DeepFRI=pd.DataFrame(DeepFRI_list)
DeepFRI_array=df_DeepFRI.values
DeepFRI_cos_sim_matrix = cosine_similarity(DeepFRI_array)
# Normalize ESM2_cos_sim_matrix
scaler = MinMaxScaler()
ESM2_cos_sim_matrix_normalized = scaler.fit_transform(ESM2_cos_sim_matrix)

# Normalize DeepFRI_cos_sim_matrix
DeepFRI_cos_sim_matrix_normalized = scaler.fit_transform(DeepFRI_cos_sim_matrix)
print(np.shape(ESM2_cos_sim_matrix_normalized))

plt.hist(ESM2_cos_sim_matrix_normalized[ESM2_cos_sim_matrix_normalized != 1],bins=100,alpha=0.4,density=True,log=False,color='blue',label='ESM2')
plt.hist(DeepFRI_cos_sim_matrix_normalized[DeepFRI_cos_sim_matrix_normalized != 1],bins=100,alpha=0.4,density=True,log=False,color='red',label='DeepFRI')
plt.xlabel("cosine similarity")
plt.ylabel("Count")
plt.legend()
plt.savefig("cosine similarity distribution",dpi=800)
plt.show()


