import pandas as pd
import subprocess
import sys

def run_DeepFRI(dataset):
    df=pd.read_csv(dataset)
    for i, row in df.iterrows():
        seq = row['Sequence'].upper()
        index = row['pro_id']
        print("seq:",seq)
        # Construct the command to run predict.py
        command = [
            'python', 'predict.py',
            '--seq', seq,
            '-ont', 'mf',
            '--verbose',
            '-o', f'DeepFRI_outputs/{index}'  # Output file for each row
        ]
        # Execute the command
        subprocess.run(command)

if __name__=='__main__':
    dataset=sys.argv[1]
    run_DeepFRI(dataset)
