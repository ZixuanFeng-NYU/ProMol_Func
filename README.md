# ProMol_Func: A Structure-Free Deep Learning Model for Virtual Screening 
<img width="1738" height="936" alt="Overview-2" src="https://github.com/user-attachments/assets/79f432f7-cf06-44c3-932b-3a9bdfef5549" />

ProMol_Func bridges small molecules graph encoder model KANO and sequence-based protein functions prediction model DeepFRI for protein-ligand binding probability prediction task. By default, the saved three models under "saved_models/0610data_5FFN_3models" is used. This is the same sets of models which were used to evluate general model performance on DUD-E,LIT-PCBA, DEKOIS2.0. Results can be found in paper


## Usage
### Install
Clone the current repo
```
git clone https://github.com/ZixuanFeng-NYU/ProMol_Func.git
```
### Set up environment 
Users need to install dependency packages (python 3.7)
```
pip3 install torch torchvision
pip install rdkit
pip install numpy==1.20.3
pip install gensim==4.2.0
pip install nltk==3.4.5
conda install -c conda-forge "owlready2>=0.25,<0.26" 
pip install Owlready2==0.37  
pip install torch-scatter==2.0.9 -f https://data.pyg.org/whl/torch-1.13.1+cu117.html
pip install tensorflow-gpu==2.3.1
pip install biopython==1.76
pip install scikit-learn==0.23.1
pip install pandas
pip install tqdm
pip install tensorboardX
pip install Unidecode
```

### Protein Functions Prediction
A pro_sequence.csv file needs to be prepared with one column ['pro_id'] and one column ['Sequence']. DeepFRI Pretrained models can be downloaded from:
https://users.flatironinstitute.org/~renfrew/DeepFRI_data/trained_models.tar.gz (run DeepFRI on GPU). Uncompress tar.gz file into the DeepFRI directory (tar xvzf trained_models.tar.gz -C /path/to/DeepFRI). 
```
cd DeepFRI
python predict_protein_functions.py sample_protein.csv
```
### General model Prediction
```
cd KANO
mkdir general_model_prediciton
python model_prediction.py --gpu 0 --test_path examples/sample_data.csv --preds_path general_model_prediction  --checkpoint_dir ../saved_models/0610data_5FFN_3models/
```
### Reproduction of evaluation results
download and decompress LIT_PCBA_EF.zip from [https://doi.org/10.5281/zenodo.16825387](https://doi.org/10.5281/zenodo.16825387) and put it under evaluate_dataset
```
cd evaluate_dataset
python LIT-PCBA_data_process.py
cd ../DeepFRI
python predict_protein_functions.py lit-pcba-protein.csv
cd ../KANO
python model_prediction.py --gpu 0 --test_path ../evaluate_dataset/LIT_PCBA_EF/data_per_target/ADRB2_protein_ligands.csv --preds_path ../evaluate_dataset/LIT_PCBA_EF/data_per_target --checkpoint_dir ../saved_models/0610data_5FFN_3models/
###After done with other targets
python LIT_PCBA_EF.py
```

### Retrain General model
Proteins used to train the general ProMol_Func model are provided in general_model_protein_sequence.csv (available at the Zenodo repository: https://doi.org/10.5281/zenodo.16825387). Protein molecular functional scores are first predicted for these protein sequences, and the resulting function embeddings are then combined with compound data from ProMol_Func_general_model_data_06102024_add_decoys.csv (from the same repository) to train the general model. The general ProMol_Func model is initialized from a pretrained KANO checkpoint, and three ensemble models are trained using different random split seeds (1, 2, and 3).

Step-by-step setup
1. Protein functional score prediction. Download general_model_protein_sequence.csv, which contains the protein sequences used for training the model, from the Zenodo repository (https://doi.org/10.5281/zenodo.16825387).
Place this file in the folder: ProMol_Func/DeepFRI/
```
cd DeepFRI
python predict_protein_functions.py general_model_protein_sequence.csv
```
This will produce output files of the form {pro_id}_MF_pred_scores.json under: ProMol_Func/DeepFRI/DeepFRI_outputs. These JSON files are then automatically read during the ProMol_Func model training process.

2. Download training data file
From the same Zenodo repository, download ProMol_Func_general_model_data_06102024_add_decoys.csv.zip.
Place this file in the folder: ProMol_Func/KANO/
In the ProMol_Func/KANO directory, unzip the file:
```
gunzip ProMol_Func_general_model_data_06102024_add_decoys.csv.zip
```
After this step, the file ProMol_Func_general_model_data_06102024_add_decoys.csv should be present in ProMol_Func/KANO/.

3. Train the general ProMol_Func model
ProMol_Func_general_model_data_06102024_add_decoys.csv, is used as inputs to train the general model.
Initialize the model from the pretrained KANO checkpoint.
Train three ensemble models using different random split seeds (1, 2, and 3). Users may adjust the ensemble size according to their needs.
```
cd KANO
python model_train.py --gpu 0 --data_path './ProMol_Func_general_model_data_06102024_add_decoys.csv' --metric accuracy --dataset_type classification --epochs 10 --gpu 0 --batch_size 256 --ensemble_size 1 --num_runs 1 --seed 1 --init_lr '1e-4' --split_type scaffold_balanced --step functional_prompt --ffn_num_layers 5 --exp_name ProMol_func_general_0610data_5FFN_1 --exp_id ProMol_func_general_0610data_5FFN_1 --checkpoint_path './original_CMPN_0623_1350_14000th_epoch.pkl'
```

### Acknowledgements
Thanks for the following released code bases:
[KANO](https://github.com/HICAI-ZJU/KANO), [DeepFRI](https://github.com/flatironinstitute/DeepFRI), [RDKit](https://github.com/rdkit/rdkit), [pdb2fasta](https://github.com/alexholehouse/pdb2fasta)

