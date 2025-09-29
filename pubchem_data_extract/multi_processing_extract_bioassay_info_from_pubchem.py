#!/usr/bin/env python3
import os
import time
import csv
import math
import json
import signal
import requests
import pandas as pd
from multiprocessing import Pool, current_process
from urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter

# -------------------
# Config
# -------------------
START_AID = 1
END_AID   = 1_964_000   # inclusive range end
PROCESSES = 8           # be nice to PubChem; avoid huge parallelism
CHUNK_SIZE = 5_000      # write to disk every N results
OUT_CSV = "pubchem_assay_type.csv"

# -------------------
# Per-process session
# -------------------
_session = None

def _init_worker():
    """Initializer for each worker process: set up a requests.Session with retries."""
    global _session
    s = requests.Session()
    retries = Retry(
        total=5,
        backoff_factor=1.0,    # 1s, 2s, 4s, 8s...
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=frozenset(["GET"])
    )
    adapter = HTTPAdapter(max_retries=retries, pool_connections=50, pool_maxsize=50)
    s.mount("http://", adapter)
    s.mount("https://", adapter)
    _session = s
def _fetch_one(aid: int):
    """Fetch one AID’s summary safely with timeouts and retries."""
    global _session
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/assay/{aid}/JSON"
    try:
        resp = _session.get(url, timeout=(5, 20))
    except Exception as e:
        # network error: return N/A row
        return {"AID": aid, "assay_type": "N/A", "Target": "N/A", "Target_URL": "N/A", "status": "net_err"}

    if resp.status_code == 404:
        # Non-existent AID
        return {"AID": aid, "assay_type": "N/A", "Target": "N/A", "Target_URL": "N/A", "status": "404"}

    if resp.status_code != 200:
        # e.g., 429 / 5xx after retries
        return {"AID": aid, "assay_type": "N/A", "Target": "N/A", "Target_URL": "N/A", "status": str(resp.status_code)}

    try:
        data = resp.json()
    except json.JSONDecodeError:
        return {"AID": aid, "assay_type": "N/A", "Target": "N/A", "Target_URL": "N/A", "status": "bad_json"}

    current_assay_type = "N/A"
    current_target = "N/A"
    current_target_url = "N/A"

    try:
        for item in data.get("Record", {}).get("Section", []):
            if item.get("TOCHeading") == "Target":
                # Guard against missing nested keys
                sec0 = item.get("Section", [{}])[0]
                info0 = sec0.get("Information", [{}])[0]
                swm = info0.get("Value", {}).get("StringWithMarkup", [{}])[0]
                current_target = swm.get("String", "N/A")
                current_target_url = info0.get("URL", "N/A")
            elif item.get("TOCHeading") == "BioAssay Annotations":
                info0 = item.get("Information", [{}])[0]
                swm = info0.get("Value", {}).get("StringWithMarkup", [{}])[0]
                current_assay_type = swm.get("String", "N/A")
    except Exception:
        # Schema irregularity: just return what we have
        pass

    return {
        "AID": aid,
        "assay_type": current_assay_type,
        "Target": current_target,
        "Target_URL": current_target_url,
        "status": "ok",
    }
def _write_header_if_needed(path):
    if not os.path.exists(path):
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["AID", "assay_type", "Target", "Target_URL", "status"])

def main():
    _write_header_if_needed(OUT_CSV)

    ids = range(START_AID, END_AID + 1)

    # Use imap_unordered to stream results as they complete.
    with Pool(processes=PROCESSES, initializer=_init_worker) as pool, \
         open(OUT_CSV, "a", newline="", encoding="utf-8") as f_out:

        writer = csv.writer(f_out)
        batch = []
        for idx, row in enumerate(pool.imap_unordered(_fetch_one, ids, chunksize=200), 1):
            batch.append(row)
            # Periodically flush to disk
            if len(batch) >= CHUNK_SIZE:
                for r in batch:
                    writer.writerow([r["AID"], r["assay_type"], r["Target"], r["Target_URL"], r["status"]])
                f_out.flush()
                batch.clear()

            # (Optional) light progress
            if idx % 10000 == 0:
                print(f"Processed {idx:,} AIDs")

        # write remaining
        if batch:
            for r in batch:
                writer.writerow([r["AID"], r["assay_type"], r["Target"], r["Target_URL"], r["status"]])

if __name__ == "__main__":
    main()
