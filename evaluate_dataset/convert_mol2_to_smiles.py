import os

pdb_id_list=os.listdir("CASF-2016/coreset")
for pdb_id in pdb_id_list:
    os.system(f"obabel -imol2 CASF-2016/coreset/{pdb_id}/{pdb_id}_ligand.mol2 -osmi -O CASF-2016/coreset/{pdb_id}/{pdb_id}_ligand.smi")
    
