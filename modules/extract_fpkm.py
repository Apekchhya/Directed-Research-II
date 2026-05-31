"""
extract_fpkm.py
---------------
Single function to extract FPKM values from all StringTie GTF files
and combine them into one CSV.
"""

import os
import re
import glob
import pandas as pd


def extract_fpkm_matrix(stringtie_dir, output_csv):
    """
    Parse all StringTie GTF files under stringtie_dir and write a combined
    FPKM matrix to output_csv.

    Parameters
    ----------
    stringtie_dir : str  – folder containing per-sample GTF files (searched recursively)
    output_csv    : str  – destination CSV path

    Returns
    -------
    pd.DataFrame  – rows = transcript_id, columns = sample names
    """

    # ── 1. Find all GTF files ──────────────────────────────────────────────
    gtf_files = sorted(glob.glob(os.path.join(stringtie_dir, "**/*.gtf"), recursive=True))
    if not gtf_files:
        # fallback: flat directory
        gtf_files = sorted(glob.glob(os.path.join(stringtie_dir, "*.gtf")))
    if not gtf_files:
        raise FileNotFoundError(f"No GTF files found under: {stringtie_dir}")

    # ── 2. Parse FPKM from each file ──────────────────────────────────────
    series_list = []
    for gtf_path in gtf_files:
        # derive sample name: prefer parent folder name, else filename stem
        parent = os.path.basename(os.path.dirname(gtf_path))
        sample = parent if parent != os.path.basename(stringtie_dir) \
                        else os.path.splitext(os.path.basename(gtf_path))[0]

        fpkm_dict = {}
        with open(gtf_path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) < 9 or fields[2] != "transcript":
                    continue
                tid   = re.search(r'transcript_id "([^"]+)"', fields[8])
                fpkm  = re.search(r'FPKM "([^"]+)"',          fields[8])
                if tid and fpkm:
                    fpkm_dict[tid.group(1)] = float(fpkm.group(1))

        s = pd.Series(fpkm_dict, name=sample, dtype=float)
        s.index.name = "transcript_id"
        series_list.append(s)
        print(f"  {sample}: {len(s):,} transcripts")

    # ── 3. Combine and save ────────────────────────────────────────────────
    df = pd.concat(series_list, axis=1).fillna(0.0).sort_index()
    os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
    df.to_csv(output_csv)
    print(f"\nSaved {df.shape[0]:,} transcripts × {df.shape[1]} samples → {output_csv}")
    return df