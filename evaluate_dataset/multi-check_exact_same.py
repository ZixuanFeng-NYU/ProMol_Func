#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pathlib import Path
import os
import pandas as pd
from multiprocessing import Pool, cpu_count

# --- Config ---
GENERAL_TRAIN_CSV = "../general_model_training_data/ProMol_Func_general_model_data_06102024_add_decoys.csv"
TARGET_DIR = "LIT_PCBA_EF/data_per_target"
N_PROCESSES = int(os.environ.get("N_PROCESSES", cpu_count()))

# --- Globals populated in workers via initializer ---
_G_GENERAL_EXACT_PAIRS = None
_G_TARGET_DIR = None


# --- Helpers ---
def canonicalize_smiles(smi: str):
    """Return RDKit canonical (isomeric) SMILES or None if invalid/empty."""
    from rdkit import Chem
    if not isinstance(smi, str) or not smi.strip():
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
  
def _init_worker(general_exact_pairs, target_dir, out_dir):
    # Store in module-level globals for fast access in each worker
    global _G_GENERAL_EXACT_PAIRS, _G_TARGET_DIR, _G_OUT_DIR
    _G_GENERAL_EXACT_PAIRS = general_exact_pairs
    _G_TARGET_DIR = Path(target_dir)
  
  

def _process_one_target(tpath_str: str):
    """
    Worker: process a single target file.
    Returns (target_name, n_dup, dup_csv_or_None, no_dup_csv_or_None, error_or_None).
    """
    try:
        tpath = Path(tpath_str)
        target_name = tpath.name.split("_protein")[0]

        df_t = pd.read_csv(tpath, usecols=["smiles", "Sequence"]).copy()
        if df_t.empty:
            # Still emit an empty "no_duplicates" file for consistency
            no_dup_path = _G_OUT_DIR / f"{target_name}_no_duplicates_protein_ligands.csv"
            df_t.to_csv(no_dup_path, index=False)
            return (target_name, 0, None, str(no_dup_path), None)

        # Canonicalize + validate
        df_t["canon_smiles"] = df_t["smiles"].apply(canonicalize_smiles)
        df_t = df_t[df_t["canon_smiles"].notna() & df_t["Sequence"].astype(bool)].reset_index(drop=True)

        # Exact duplicates against GENERAL
        pairs = list(zip(df_t["Sequence"].tolist(), df_t["canon_smiles"].tolist()))
        is_dup = [pair in _G_GENERAL_EXACT_PAIRS for pair in pairs]

        df_dup = df_t[is_dup].reset_index(drop=True)
        df_no_dup = df_t[~pd.Series(is_dup)].reset_index(drop=True)

        # Write duplicates (if any) beside the target file
        dup_path = None
        if not df_dup.empty:
            dup_path = _G_TARGET_DIR / f"{target_name}_duplicate.csv"
            df_dup[["smiles", "Sequence"]].to_csv(dup_path, index=False)

        
        return (target_name, len(df_dup), str(dup_path) if dup_path else None)

    except Exception as e:
        return (Path(tpath_str).stem, 0, None, None, str(e))

def main():
    # Load GENERAL and build exact-identity set: (Sequence, canon_smiles)
    General = pd.read_csv(GENERAL_TRAIN_CSV, usecols=["smiles", "Sequence"]).copy()
    General["canon_smiles"] = General["smiles"].apply(canonicalize_smiles)
    General = General[General["canon_smiles"].notna() & General["Sequence"].astype(bool)].reset_index(drop=True)

    general_exact_pairs = set(zip(General["Sequence"], General["canon_smiles"]))
    print(
        f"GENERAL: {len(General)} valid (Sequence, canonical SMILES) rows; "
        f"{len(general_exact_pairs)} unique exact pairs."
    )

    # Collect target files
    tdir = Path(TARGET_DIR)
    target_files = [str(p) for p in tdir.iterdir() if p.name.endswith("protein_ligands.csv")]
    if not target_files:
        print(f"No files ending with 'protein_ligands.csv' in {tdir}")
        return

    print(f"Found {len(target_files)} target files. Using {N_PROCESSES} processes.")

    # Parallel map
    with Pool(
        processes=N_PROCESSES,
        initializer=_init_worker,
        initargs=(general_exact_pairs, str(tdir), NO_DUP_DIR),
    ) as pool:
        for target_name, n_dup, dup_path, no_dup_path, err in pool.imap_unordered(_process_one_target, target_files):
            if err:
                print(f"[{target_name}] ERROR: {err}")
            else:
                if n_dup == 0:
                    print(f"[{target_name}] No exact duplicates with GENERAL found.")
                else:
                    print(f"[{target_name}] Found {n_dup} exact duplicates. Saved to {Path(dup_path).name}")
                print(f"[{target_name}] Wrote no-duplicate file: {Path(no_dup_path).name}")

if __name__ == "__main__":
    main()


         
