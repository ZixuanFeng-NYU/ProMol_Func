import warnings
warnings.filterwarnings('ignore')
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from argparse import Namespace
from logging import Logger
import os
from typing import Tuple
import numpy as np
import pandas as pd

from chemprop.train.make_predictions_add import make_predictions
from chemprop.train.run_training_with_add import run_training
from chemprop.data_with_pro.utils import get_task_names
from chemprop.utils import makedirs
from chemprop.parsing import parse_predict_args
from chemprop.torchlight import initialize_exp
from sklearn.metrics import precision_score, recall_score,f1_score,confusion_matrix
from sklearn.metrics import accuracy_score

def run_stat(args: Namespace) -> Tuple[list, list]:
    """Run predictions"""

    # Call the make_predictions function and collect the predictions and smiles
    model_preds, smiles,pro_ids,targets = make_predictions(args)
    test_path=args.test_path
    file_part=test_path.split('/')[-1].split('.')[0]

    return model_preds,smiles,pro_ids,targets,file_part


if __name__ == '__main__':
    args = parse_predict_args()

    # Call the run_stat function to get the predictions and smiles
    preds,smiles,pro_ids,targets,file_part= run_stat(args)

    df_out=pd.DataFrame()
    df_out['smiles']=smiles
    Targets=[]    
    for target in targets:
        if target:
            Targets.append(target[0])
        else:
            Targets.append(None)
    df_out['Targets']=Targets
    df_out['pro_id']=pro_ids
    df_out['preds'] = [pred[0] for pred in preds]
    df_out['preds']=df_out['preds'].astype(float) 
    save_path=args.preds_path
    df_out.to_csv(os.path.join(save_path,f'{file_part}_ProMol_Func_prediction.csv'))
