# ProMol_Func: Protein functions talk to small Molecules functional groups

[ProMol_Func__general.pdf](https://github.com/user-attachments/files/19197139/ProMol_Func__general.pdf)




ProMol_Func bridges small molecules graph encoder model KANO and structure-based protein functions prediction model DeepFRI for protein-ligand binders classification/regression tasks. By default, the saved three models under "saved_models/0610data_5FFN_3models" is used. This is the same sets of models which were used to evluate general model performance on DUD-E,LIT-PCBA, DEKOIS2.0, and CASF-2016. Results can be found in paper


## Usage
### Install
We will set up the environment using Anaconda. Clone the current repo
```
git clone https://github.com/ZixuanFeng-NYU/ProMol_Func.git
```
### Protein Functions Prediction
A pro_sequence.csv file needs to be prepared with one column ['pro_id'] and one column ['Sequence']. DeepFRI Pretrained models can be downloaded from:
https://users.flatironinstitute.org/~renfrew/DeepFRI_data/trained_models.tar.gz (run DeepFRI on GPU). Uncompress tar.gz file into the DeepFRI directory (tar xvzf trained_models.tar.gz -C /path/to/DeepFRI).
```
cd DeepFRI
mkdir DeepFRI_outputs
python predict_protein_functions.py sample_protein.csv
```
### General model Prediction
```
cd KANO
mkdir general_model_prediciton
python model_prediction.py --gpu 0 --test_path examples/sample_data.csv --preds_path general_model_prediction  --checkpoint_dir ../saved_models/0610data_5FFN_3models/
```
### Finetune on Target-specific task
Original chaperone data we collected inlucdes 568578 pairs, here we only provide sample chaperone data (data used for validation in the paper) including ligands for HSP70, EcoliDnaK, HSP90, and some protein with close protein functions. Before finetuning, you need to predict protein functions of target proteins.
```
cd DeepFRI
python predict_protein_functions.py chaperone_protein_seq.csv
cd ../KANO
mkdir Target_specific_model
python model_train.py  --data_path ./examples/chaperone_sample_data.csv    --metric 'accuracy'        --dataset_type classification --epochs 40    --gpu 0    --batch_size 256    --ensemble_size 1    --num_runs 1  --seed 1  --init_lr 1e-4    --split_type 'scaffold_balanced'         --step 'functional_prompt'    --ffn_num_layers 5  --split_sizes 0.8 0.1 0.1  --exp_name finetune      --exp_id finetune         --checkpoint_path "../saved_models/0610data_5FFN_3models/ProMol_func_general_0610data_5FFN_1/run_0/model_0/model.pt"
```
### Acknowledgements
Thanks for the following released code bases:
[KANO](https://github.com/HICAI-ZJU/KANO), [DeepFRI](https://github.com/flatironinstitute/DeepFRI), [RDKit](https://github.com/rdkit/rdkit), [pdb2fasta](https://github.com/alexholehouse/pdb2fasta)

