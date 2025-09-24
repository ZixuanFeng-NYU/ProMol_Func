#!/usr/bin/env python3
import os
from multiprocessing import Pool, cpu_count
import pandas as pd

# --- Config ---
DATA_DIR = "LIT_PCBA_EF/data_per_target"
# If you want to match only on certain columns (e.g., SMILES), set:
# KEY_COLUMNS = ["smiles"]
KEY_COLUMNS = None  # use all shared columns by default
N_PROCESSES = int(os.environ.get("N_PROCESSES", cpu_count()))

def process_over_file(over_filename: str) -> str:
    """
    For an {target}_over{threshold}.csv file, subtract rows from {target}_protein_ligands.csv
    and write {target}_no_{threshold}.csv. Returns a status string.
    """
    try:
        base = over_filename[:-4]  # drop .csv
        target, threshold = base.split("_over", 1)
        protein_file = os.path.join(DATA_DIR, f"{target}_protein_ligands.csv")
        over_file = os.path.join(DATA_DIR, over_filename)
        out_file = os.path.join(DATA_DIR, f"{target}_no_{threshold}.csv")

        if not os.path.exists(protein_file):
            return f"⚠️ Skip {over_filename}: missing {os.path.basename(protein_file)}"

        # Load data
        df_protein = pd.read_csv(protein_file)
        df_over = pd.read_csv(over_file)

        # Decide merge keys
        if KEY_COLUMNS is None:
            # Use intersection of columns (safe even if column order differs)
            keys = [c for c in df_protein.columns if c in df_over.columns]
            if not keys:
                return f"⚠️ Skip {over_filename}: no shared columns between protein and over files."
            # Use keys to avoid massive full-row matching artifacts
            on = keys
        else:
            missing = [k for k in KEY_COLUMNS if k not in df_protein.columns or k not in df_over.columns]
            if missing:
                return f"⚠️ Skip {over_filename}: missing key columns {missing}"
            on = KEY_COLUMNS

        # Drop duplicates on keys to make anti-join cleaner
        df_over_keys = df_over[on].drop_duplicates()

        # Anti-join: keep rows from df_protein whose keys are NOT in df_over
        # Use a left merge with indicator to filter left_only
        df_remaining = (
            df_protein.merge(df_over_keys.assign(_marker=1), how="left", on=on)
                      .loc[lambda d: d["_marker"].isna()]
                      .drop(columns=["_marker"])
        )

        # Save
        df_remaining.to_csv(out_file, index=False)
        return f"✅ {os.path.basename(out_file)} saved ({len(df_remaining)} rows)"

    except Exception as e:
        return f"❌ Error processing {over_filename}: {e}"
def main():
    files = os.listdir(DATA_DIR)
    over_files = [f for f in files if f.endswith(".csv") and "_over" in f]

    if not over_files:
        print("No *_over*.csv files found.")
        return

    print(f"Using {N_PROCESSES} processes...")
    with Pool(processes=N_PROCESSES) as pool:
        for msg in pool.imap_unordered(process_over_file, over_files):
            print(msg)

if __name__ == "__main__":
    main()

