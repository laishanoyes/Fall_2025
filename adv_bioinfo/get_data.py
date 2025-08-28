#!/usr/bin/env python3
"""
get_us_sars2_accessions.py

Download SARS-CoV-2 metadata from NCBI (taxon 2697049),
filter to USA samples, and select random accession IDs.
"""

import subprocess
import pandas as pd
import numpy as np
import sys
from pathlib import Path

def main(n_samples=150, outdir="sars2_meta", seed=42):
    outdir = Path(outdir)
    zip_file = outdir.with_suffix(".zip")
    tsv_file = outdir.with_suffix(".tsv")



    print("[INFO] Converting metadata to TSV...")
    with open(tsv_file, "w") as f_out:
        subprocess.run([
            "dataformat", "tsv", "virus-genome",
            "--package", str(zip_file),
            #Need to update the below once we know the actual fields
			"--fields", "accession,virus-pangolin,isolate-collection-date,geo-location"
        ], check=True, stdout=f_out)

    print("[INFO] Filtering to USA samples...")
    df = pd.read_csv(tsv_file, sep="\t", dtype=str)
    mask = df["geo-location"].fillna("").str.contains("USA|United States", case=False)
    us = df[mask].copy()
    print(f"[INFO] Found {len(us)} US records")

    if us.empty:
        raise RuntimeError("No USA records found in metadata.")


    rng = np.random.default_rng(seed)
    n = min(n_samples, len(us))
    sampled = us.sample(n=n, random_state=int(rng.integers(0, 1_000_000)))


    sampled["accession"].to_csv("us_sars2_accessions.txt", index=False, header=False)
    sampled.to_csv("us_sars2_audit.tsv", sep="\t", index=False)

    print(f"[INFO] Selected {len(sampled)} US accessions")
    print("  - us_sars2_accessions.txt (list of IDs)")
    print("  - us_sars2_audit.tsv (with lineage, date, location)")

if __name__ == "__main__":

    n = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    main(n_samples=n)

