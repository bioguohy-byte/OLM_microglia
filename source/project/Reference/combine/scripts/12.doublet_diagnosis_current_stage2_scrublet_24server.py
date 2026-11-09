#!/usr/bin/env python3
"""Stage 3, part 2: Scrublet on current reanalysis raw-count partitions.

Run this program in the visible RStudio Terminal.  All cache directories are
inside the project result directory, never in the account home directory.
"""

import importlib.metadata
import json
import os
from pathlib import Path

PROJECT_ROOT = Path(__import__("os").path.join(__import__("os").environ["MG_PROJECT_ROOT"], "Reference"))
RESULT_ROOT = PROJECT_ROOT / "combine/results/full_count_reanalysis/03_current_MG0_doublet_diagnosis"
RUNTIME_DIR = RESULT_ROOT / "runtime"

os.environ["PYTHONPYCACHEPREFIX"] = str(RUNTIME_DIR / "python_pycache")
os.environ["MPLCONFIGDIR"] = str(RUNTIME_DIR / "matplotlib")
os.environ["NUMBA_CACHE_DIR"] = str(RUNTIME_DIR / "numba_cache")
for cache_path in (Path(os.environ["PYTHONPYCACHEPREFIX"]), Path(os.environ["MPLCONFIGDIR"]), Path(os.environ["NUMBA_CACHE_DIR"])):
    cache_path.mkdir(parents=True, exist_ok=True)

import numpy as np
import pandas as pd
import scipy.io
import scrublet as scr


def main() -> None:
    input_root = RESULT_ROOT / "scrublet_input"
    output_root = RESULT_ROOT / "scrublet_output"
    output_root.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(input_root / "scrublet_manifest.csv")
    summaries = []
    for row in manifest.itertuples(index=False):
        partition_dir = input_root / row.safe_partition
        counts_cells_by_genes = scipy.io.mmread(partition_dir / "counts_genes_by_cells.mtx").tocsc().T
        cells = pd.read_csv(partition_dir / "cells.tsv", sep="\t")
        if counts_cells_by_genes.shape[0] != len(cells):
            raise RuntimeError(f"{row.DF_partition}: matrix/cell list mismatch")
        n_components = min(30, counts_cells_by_genes.shape[0] - 1, counts_cells_by_genes.shape[1] - 1)
        if n_components < 2:
            raise RuntimeError(f"{row.DF_partition}: too few cells/features for Scrublet")
        scrub = scr.Scrublet(
            counts_cells_by_genes,
            expected_doublet_rate=float(row.expected_doublet_rate),
            random_state=20260907,
        )
        scores, calls = scrub.scrub_doublets(
            min_counts=2, min_cells=3, min_gene_variability_pctl=85,
            n_prin_comps=n_components,
        )
        result = pd.DataFrame({
            "Cell": cells["cell_id"].astype(str),
            "DF_partition": row.DF_partition,
            "Scrublet_score": scores,
            "Scrublet_call": np.where(calls, "Doublet", "Singlet"),
            "Scrublet_threshold": float(scrub.threshold_),
            "expected_doublet_rate": float(row.expected_doublet_rate),
        })
        result.to_csv(output_root / f"Scrublet_{row.safe_partition}.csv", index=False)
        summaries.append({
            "DF_partition": row.DF_partition,
            "n_cells": int(counts_cells_by_genes.shape[0]),
            "n_genes": int(counts_cells_by_genes.shape[1]),
            "n_prin_comps": int(n_components),
            "expected_doublet_rate": float(row.expected_doublet_rate),
            "threshold": float(scrub.threshold_),
            "n_called_doublet": int(calls.sum()),
        })
        print(f"{row.DF_partition}: Scrublet complete ({int(calls.sum())}/{counts_cells_by_genes.shape[0]} doublets)", flush=True)
    pd.DataFrame(summaries).to_csv(output_root / "Scrublet_run_summary.csv", index=False)
    with (RUNTIME_DIR / "stage2_scrublet_versions.json").open("w") as handle:
        json.dump({
            "python": os.sys.version,
            "scrublet": importlib.metadata.version("scrublet"),
            "numpy": importlib.metadata.version("numpy"),
            "scipy": importlib.metadata.version("scipy"),
            "scikit_learn": importlib.metadata.version("scikit-learn"),
        }, handle, indent=2)


if __name__ == "__main__":
    main()
