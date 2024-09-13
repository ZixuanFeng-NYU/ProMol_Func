import warnings
warnings.filterwarnings('ignore')
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from argparse import Namespace
from logging import Logger
from typing import Tuple
import pandas as pd

from chemprop.train.make_predictions_add import get_embs
from chemprop.data_with_pro.utils import get_task_names
from chemprop.parsing import parse_predict_args


def run_stat(args: Namespace) -> Tuple[list, list]:
    """Run predictions"""

    # Call the make_predictions function and collect the predictions and smiles
    model_preds, smiles,pro_ids,targets = get_embs(args)
    test_path=args.test_path
    file_part=test_path.split('/')[-1].split('.')[0]

    return model_preds,smiles,pro_ids,targets,file_part


if __name__ == '__main__':
    args = parse_predict_args()

    # Call the run_stat function to get the predictions and smiles
    preds,smiles,pro_ids,targets,file_part= run_stat(args)
    print(preds)
    df=pd.DataFrame(preds)
    df.to_csv(f"{file_part}_emb_after_general_model_training.csv")
    print("model_embs for small molecule shape:",len(preds),len(preds[0]))
