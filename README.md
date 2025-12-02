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
We provide a ready-to-use overlay environment via Zenodo (https://doi.org/10.5281/zenodo.16825387
).
Please note that DeepFRI dependencies are installed in a dedicated environment; users should activate deepfri env when performing protein-function prediction.

Alternatively, users may create their own environments and install all required dependencies (Python 3.7 and associated packages).
For best practice, we strongly recommend installing DeepFRI in a separate conda environment to avoid conflicts and ensure reproducibility.
```
# 1. Create and activate a fresh conda environment (recommended)
conda create -n ProMol_Func python=3.7 -y
conda activate ProMol_Func

# 2. Core dependencies
pip install torch torchvision          
pip install rdkit                   
pip install numpy==1.20.3
pip install gensim==4.2.0
pip install nltk==3.4.5
pip install Owlready2==0.37
pip install torch-scatter==2.0.9 -f https://data.pyg.org/whl/torch-1.13.1+cu117.html
pip install pandas
pip install tqdm
pip install tensorboardX
pip install Unidecode

# 3. Create and activate the DeepFRI environment
conda create -n deepfri python=3.7 -y
conda activate deepfri

# 4. DeepFRI requirements
pip install tensorflow-gpu==2.3.1
pip install biopython==1.76
pip install scikit-learn==0.23.1

```

### Protein Functions Prediction
A protein sequence file, pro_sequence.csv (e.g., sample_protein.csv), must be prepared with two columns: 'pro_id' and 'Sequence'.
DeepFRI pretrained models can be downloaded from:
https://users.flatironinstitute.org/~renfrew/DeepFRI_data/trained_models.tar.gz
(DeepFRI requires a GPU to run efficiently.)
After downloading, extract the pretrained models into the DeepFRI directory: 
```
tar xvzf trained_models.tar.gz -C /path/to/DeepFRI
```
```
conda activate deepfri
cd DeepFRI
python predict_protein_functions.py sample_protein.csv

```
This will produce output files of the form {pro_id}_MF_pred_scores.json under: ProMol_Func/DeepFRI/DeepFRI_outputs. These JSON files are then automatically read during the ProMol_Func inference process.

### General model Prediction
Prepare a protein_ligand.csv file (e.g., sample_data.csv) that includes three required columns: 'smiles', 'Class', and 'pro_id'.

  'smiles': lists the SMILES strings of the compounds.
  
  'Class': specifies the known binding class (e.g., active/inactive). This column may be left empty for real-world screening scenarios where binding labels are not available.
  
  'pro_id': specifies the protein identifier and is used to retrieve the corresponding protein functional scores.
ProMol_Func requires a GPU to run efficiently.
```
cd KANO
mkdir general_model_prediciton
python model_prediction.py --gpu 0 --test_path examples/sample_data.csv --preds_path general_model_prediction  --checkpoint_dir ../saved_models/0610data_5FFN_3models/
```
After model prediction, a new column named 'preds' will be generated, containing the predicted protein–small-molecule binding probabilities.
### Reproduction of evaluation results
download and decompress LIT_PCBA_EF.zip from [https://doi.org/10.5281/zenodo.16825387](https://doi.org/10.5281/zenodo.16825387) and put it under evaluate_dataset.
ProMol_Func requires a GPU to run efficiently.
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
ProMol_Func requires a GPU to run efficiently.

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
Train three ensemble models using different random split seeds (1, 2, and 3). Users may adjust the hyperparameters such as ensemble_size, batch_size and metric, according to their needs.
```
cd KANO
python model_train.py --gpu 0 --data_path './ProMol_Func_general_model_data_06102024_add_decoys.csv' --metric accuracy --dataset_type classification --epochs 10 --gpu 0 --batch_size 256 --ensemble_size 1 --num_runs 1 --seed 1 --init_lr '1e-4' --split_type scaffold_balanced --step functional_prompt --ffn_num_layers 5 --exp_name ProMol_func_general_0610data_5FFN_1 --exp_id ProMol_func_general_0610data_5FFN_1 --checkpoint_path './original_CMPN_0623_1350_14000th_epoch.pkl'
```

### Acknowledgements
Thanks for the following released code bases:
[KANO](https://github.com/HICAI-ZJU/KANO), [DeepFRI](https://github.com/flatironinstitute/DeepFRI), [RDKit](https://github.com/rdkit/rdkit), [pdb2fasta](https://github.com/alexholehouse/pdb2fasta)

