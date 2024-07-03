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

#Plotting ESM2_cos_sim_matrix
plt.figure(figsize=(8, 6))
plt.imshow(ESM2_cos_sim_matrix_normalized, cmap='viridis', interpolation='nearest')
plt.title('ESM2 Cosine Similarity Matrix')
plt.colorbar()
plt.savefig('ESM2_cos_sim_matrix.png')  # Save the figure as PNG
plt.show()

# Plotting DeepFRI_cos_sim_matrix
plt.figure(figsize=(8, 6))
plt.imshow(DeepFRI_cos_sim_matrix_normalized, cmap='viridis', interpolation='nearest')
plt.title('DeepFRI Cosine Similarity Matrix')
plt.colorbar()
plt.savefig('DeepFRI_cos_sim_matrix.png')  # Save the figure as PNG
plt.show()

# Reorder ESM2_cos_sim_matrix_normalized
ESM2_reordered_indices = np.argsort(-np.sum(ESM2_cos_sim_matrix_normalized, axis=1))
ESM2_cos_sim_matrix_reordered = ESM2_cos_sim_matrix_normalized[ESM2_reordered_indices]
ESM2_cos_sim_matrix_reordered = ESM2_cos_sim_matrix_reordered[:, ESM2_reordered_indices]

# Reorder DeepFRI_cos_sim_matrix_normalized
DeepFRI_reordered_indices = np.argsort(-np.sum(DeepFRI_cos_sim_matrix_normalized, axis=1))
DeepFRI_cos_sim_matrix_reordered = DeepFRI_cos_sim_matrix_normalized[DeepFRI_reordered_indices]
DeepFRI_cos_sim_matrix_reordered = DeepFRI_cos_sim_matrix_reordered[:, DeepFRI_reordered_indices]

# Plotting reordered ESM2_cos_sim_matrix
plt.figure(figsize=(8, 6))
plt.imshow(ESM2_cos_sim_matrix_reordered, cmap='viridis', interpolation='nearest')
plt.title('Reordered ESM2 Cosine Similarity Matrix')
plt.colorbar()
plt.savefig('Reordered_ESM2_cos_sim_matrix.png')  # Save the figure as PNG
plt.show()

# Plotting reordered DeepFRI_cos_sim_matrix
plt.figure(figsize=(8, 6))
plt.imshow(DeepFRI_cos_sim_matrix_reordered, cmap='viridis', interpolation='nearest')
plt.title('Reordered DeepFRI Cosine Similarity Matrix')
plt.colorbar()
plt.savefig('Reordered_DeepFRI_cos_sim_matrix.png')  # Save the figure as PNG
plt.show()
