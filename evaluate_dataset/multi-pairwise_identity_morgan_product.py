#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
from pathlib import Path
import pandas as pd
from rdkit import DataStructs

# --- Config ---
GENERAL_TRAIN_CSV = "../general_model_training_data/ProMol_Func_general_model_data_06102024_add_decoys.csv"
TARGET_DIR = "LIT_PCBA_EF/data_per_target"
THRESHOLDS = [0.90, 0.70, 0.50]
PROT_IDENTITY_CUTOFF = 0.50
LIG_TANIMOTO_CUTOFF = 0.50

# RDKit Morgan FP params
MORGAN_RADIUS = 2
MORGAN_NBITS = 2048

# Alignment scoring
MATCH_SCORE = 1.0
MISMATCH_SCORE = 0.0
GAP_OPEN = -1.0
GAP_EXTEND = -0.5
# ------------------------------------------------------------

# Globals that workers will read (set via pool initializer)
_G_GENERAL_FPS = None
_G_GENERAL_SMILES = None
_G_GENERAL_SEQ = None
_G_SEQ_IDENTITY = None
_G_THRESHOLDS = THRESHOLDS
_G_LIG_TANIMOTO_CUTOFF = LIG_TANIMOTO_CUTOFF
def _pool_init(general_fps, general_smiles, general_seq, seq_identity, thresholds, lig_tani_cutoff):
    global _G_GENERAL_FPS, _G_GENERAL_SMILES, _G_GENERAL_SEQ, _G_SEQ_IDENTITY, _G_THRESHOLDS, _G_LIG_TANIMOTO_CUTOFF
    _G_GENERAL_FPS = general_fps
    _G_GENERAL_SMILES = general_smiles
    _G_GENERAL_SEQ = general_seq
    _G_SEQ_IDENTITY = seq_identity
    _G_THRESHOLDS = thresholds
    _G_LIG_TANIMOTO_CUTOFF = lig_tani_cutoff

# --- Functions (defined once, outside) ---
def compute_identity(seq1: str, seq2: str):
    """
    Global alignment identity with Biopython pairwise2:
      identity = matches / alignment_length (including gaps)
    """
    from Bio import pairwise2
    aln = pairwise2.align.globalms(
        seq1, seq2,
        MATCH_SCORE, MISMATCH_SCORE,
        GAP_OPEN, GAP_EXTEND,
        one_alignment_only=True
    )[0]
    a, b = aln.seqA, aln.seqB
    matches = sum(1 for x, y in zip(a, b) if x == y)
    ident = matches / len(a) if len(a) > 0 else 0.0
    return ident

def morgan_fp_from_smiles(smi):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import pandas as pd

    if pd.isna(smi) or not isinstance(smi, str) or not smi.strip():
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, MORGAN_RADIUS, nBits=MORGAN_NBITS)

def _parallel_identity_worker(args):
    target_sequence, seq = args
    if not seq:
        return seq, 0.0
    ident = compute_identity(target_sequence, seq)
    return seq, ident

def _parallel_target_worker(row):
    """
    row: (idx, target_smiles, target_fp)
    Returns: dict with keys 'over09','over07','over05' -> list of (smiles, sequence) or original row info
    """
    idx, t_smi, t_fp = row
    out = {"over09": [], "over07": [], "over05": []}
    if t_fp is None:
        return out

    # Bulk Tanimoto against all general FPs
    tanis = DataStructs.BulkTanimotoSimilarity(t_fp, _G_GENERAL_FPS)

    # Evaluate thresholds only where tanimoto passes lig cutoff and identity map passes cutoff
    for j, tani in enumerate(tanis):
        if tani <= _G_LIG_TANIMOTO_CUTOFF:
            continue
        g_seq = _G_GENERAL_SEQ[j]
        prot_ident = _G_SEQ_IDENTITY.get(g_seq, 0.0)
        if prot_ident < PROT_IDENTITY_CUTOFF:
            continue

        prod = prot_ident * tani
        # Place into highest matching bucket
        if prod > 0.9:
            out["over09"].append(idx)
        if prod > 0.7:
            out["over07"].append(idx)
        if prod > 0.5:
            out["over05"].append(idx)

    return out

def main():
    import multiprocessing as mp

    # ---- Load GENERAL once ----
    General = pd.read_csv(GENERAL_TRAIN_CSV, usecols=["smiles", "Sequence"]).copy()

    # Deduplicate protein sequences for identity screening
    general_unique_seqs = sorted({s for s in General["Sequence"] if isinstance(s, str) and s})

    # ---- Iterate target files ----
    tdir = Path(TARGET_DIR)
    target_files = [p for p in tdir.iterdir() if p.name.endswith("protein_ligands.csv")]
    if not target_files:
        print(f"No files ending with 'protein_ligands.csv' in {tdir}")
        return

    for tpath in target_files:
        target_name = tpath.name.split("_protein")[0]
        print(f"\n=== Processing {tpath.name} ===")
        df_target = pd.read_csv(tpath, usecols=["smiles", "Sequence"]).copy()

        if df_target.empty:
            print("  (skip) target file is empty")
            continue

        # Assume each target file contains only one protein sequence
        target_sequence = df_target["Sequence"].values[0]

        # 1) Compute identity(target_seq, general_seq) for all unique general sequences (PARALLEL)
        with mp.Pool(processes=max(1, mp.cpu_count() - 1)) as pool:
            ident_items = pool.map(
                _parallel_identity_worker,
                [(target_sequence, s) for s in general_unique_seqs]
            )
        # Map and filter by cutoff
        seq_identity_full = dict(ident_items)
        similar_sequences = {s for s, idv in ident_items if idv >= PROT_IDENTITY_CUTOFF}

        # Filter general dataset to matching similar sequences
        df_general_similar = General[General['Sequence'].isin(similar_sequences)].copy()
        if df_general_similar.empty:
            print("  (skip) no similar protein sequences above cutoff.")
            continue

        # 2) Precompute Morgan FPs (vectorized over both tables)
        df_general_similar["fp"] = df_general_similar["smiles"].apply(morgan_fp_from_smiles)
        df_general_similar = df_general_similar[~df_general_similar["fp"].isna()].reset_index(drop=True)

        df_target["fp"] = df_target["smiles"].apply(morgan_fp_from_smiles)
        df_target = df_target[~df_target["fp"].isna()].reset_index(drop=True)

        if df_general_similar.empty or df_target.empty:
            print("  (skip) no valid fingerprints after parsing.")
            continue

        # Pack arrays for workers
        general_fps = list(df_general_similar["fp"].values)
        general_smiles = list(df_general_similar["smiles"].values)
        general_seq = list(df_general_similar["Sequence"].values)

        # Identity map (for ONLY sequences in df_general_similar)
        seq_identity_subset = {s: seq_identity_full.get(s, 0.0) for s in set(general_seq)}

        # 3) Parallel over target rows
        rows = list(df_target[["smiles", "fp"]].itertuples(index=True, name=None))  # (idx, smiles, fp)
        rows = [(idx, smi, fp) for (idx, smi, fp) in rows]
        with mp.Pool(
            processes=max(1, mp.cpu_count() - 1),
            initializer=_pool_init,
            initargs=(general_fps, general_smiles, general_seq, seq_identity_subset, THRESHOLDS, LIG_TANIMOTO_CUTOFF)
        ) as pool:
            results = pool.map(_parallel_target_worker, rows)

        # 4) Collect per-threshold unique target rows
        over09_idx = set()
        over07_idx = set()
        over05_idx = set()
        for res in results:
            over09_idx.update(res["over09"])
            over07_idx.update(res["over07"])
            over05_idx.update(res["over05"])

        # Convert to DataFrames and save
        if over09_idx:
            pd.DataFrame(df_target.loc[sorted(over09_idx), ["smiles", "Sequence"]]).to_csv(f"LIT_PCBA_EF/data_per_target/{target_name}_over09.csv", index=False)
            print(f"  Wrote {len(over09_idx)} rows to {target_name}_over09.csv")
        if over07_idx:
            pd.DataFrame(df_target.loc[sorted(over07_idx), ["smiles", "Sequence"]]).to_csv(f"LIT_PCBA_EF/data_per_target/{target_name}_over07.csv", index=False)
            print(f"  Wrote {len(over07_idx)} rows to {target_name}_over07.csv")
        if over05_idx:
            pd.DataFrame(df_target.loc[sorted(over05_idx), ["smiles", "Sequence"]]).to_csv(f"LIT_PCBA_EF/data_per_target/{target_name}_over05.csv", index=False)
            print(f"  Wrote {len(over05_idx)} rows to {target_name}_over05.csv")

if __name__ == "__main__":
    main()
