# ProMol_Func: Protein functions talk to small Molecules functional groups
<img width="1738" height="936" alt="Overview-2" src="https://github.com/user-attachments/assets/79f432f7-cf06-44c3-932b-3a9bdfef5549" />

ProMol_Func bridges small molecules graph encoder model KANO and sequence-based protein functions prediction model DeepFRI for protein-ligand binding probability prediction task. By default, the saved three models under "saved_models/0610data_5FFN_3models" is used. This is the same sets of models which were used to evluate general model performance on DUD-E,LIT-PCBA, DEKOIS2.0. Results can be found in paper


## Usage
### Install
We will set up the environment using Anaconda. Clone the current repo
```
git clone https://github.com/ZixuanFeng-NYU/ProMol_Func.git
```
### Protein Functions Prediction
A pro_sequence.csv file needs to be prepared with one column ['pro_id'] and one column ['Sequence']. DeepFRI Pretrained models can be downloaded from:
https://users.flatironinstitute.org/~renfrew/DeepFRI_data/trained_models.tar.gz (run DeepFRI on GPU). Uncompress tar.gz file into the DeepFRI directory (tar xvzf trained_models.tar.gz -C /path/to/DeepFRI). (DeepFRI needs a separate conda env with python3.7)
```
cd DeepFRI
mkdir DeepFRI_outputs
conda activate deepfri_env
python predict_protein_functions.py sample_protein.csv
```
### General model Prediction
```
cd KANO
mkdir general_model_prediciton
python model_prediction.py --gpu 0 --test_path examples/sample_data.csv --preds_path general_model_prediction  --checkpoint_dir ../saved_models/0610data_5FFN_3models/
```
### Reproduction of evaluation results
download and decompress LIT-PCBA.zip from [https://doi.org/10.5281/zenodo.16825387](https://doi.org/10.5281/zenodo.16825387) and put it under evaluate_dataset
```
cd evaluate_dataset
python LIT-PCBA_data_process.py
cd ../DeepFRI
conda activate deepfri_env 
python predict_protein_functions.py lit-pcba-protein.csv
cd ../KANO
conda deactivate
python model_prediction.py --gpu 0 --test_path ../evaluate_dataset/LIT_PCBA_EF/data_per_target/ADRB2_protein_ligands.csv --preds_path ../evaluate_dataset/LIT_PCBA_EF/data_per_target --checkpoint_dir ../saved_models/0610data_5FFN_3models/
###After done with other targets
python LIT_PCBA_EF.py
```

### Retrain General model
ProMol_Func_general_model_data_06102024_add_decoys.csv is given in [https://doi.org/10.5281/zenodo.16825387](https://doi.org/10.5281/zenodo.16825387)
### Acknowledgements
Thanks for the following released code bases:
[KANO](https://github.com/HICAI-ZJU/KANO), [DeepFRI](https://github.com/flatironinstitute/DeepFRI), [RDKit](https://github.com/rdkit/rdkit), [pdb2fasta](https://github.com/alexholehouse/pdb2fasta)

