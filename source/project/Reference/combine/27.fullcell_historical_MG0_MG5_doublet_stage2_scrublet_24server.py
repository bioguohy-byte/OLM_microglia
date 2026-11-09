#!/usr/bin/env python3
"""Corrected full-cell diagnosis, stage 2: Scrublet per study/batch partition.

Every input row comes from the all-QC-passed real-cell universe saved by stage 1.
The historical MG0--MG5 labels are not provided to Scrublet and do not affect its
simulation, score, threshold or call.
"""

import importlib.metadata
import json
import os
from pathlib import Path

PROJECT_ROOT = Path(__import__("os").path.join(__import__("os").environ["MG_PROJECT_ROOT"], "Reference"))
RESULT_ROOT = PROJECT_ROOT / "combine/results/fullcell_historical_MG0_MG5_doublet_diagnosis_v2"
RUNTIME_DIR = RESULT_ROOT / "runtime"
os.environ["PYTHONPYCACHEPREFIX"] = str(RUNTIME_DIR / "python_pycache")
os.environ["MPLCONFIGDIR"] = str(RUNTIME_DIR / "matplotlib")
os.environ["NUMBA_CACHE_DIR"] = str(RUNTIME_DIR / "numba_cache")
for folder in (Path(os.environ["PYTHONPYCACHEPREFIX"]), Path(os.environ["MPLCONFIGDIR"]), Path(os.environ["NUMBA_CACHE_DIR"])):
    folder.mkdir(parents=True, exist_ok=True)

import numpy as np
import pandas as pd
import scipy.io
import scrublet as scr


def main() -> None:
    input_root = RESULT_ROOT / "scrublet_input"
    output_root = RESULT_ROOT / "scrublet_output"
    output_root.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(input_root / "scrublet_manifest.csv")
    summary = []
    for row in manifest.itertuples(index=False):
        folder = input_root / row.safe_partition
        cells_by_genes = scipy.io.mmread(folder / "counts_genes_by_cells.mtx").tocsc().T
        cells = pd.read_csv(folder / "cells.tsv", sep="\t")
        if cells_by_genes.shape[0] != len(cells) or cells_by_genes.shape[0] != int(row.n_cells):
            raise RuntimeError(f"{row.Study} / {row.DF_partition}: all-QC-passed count matrix and manifest disagree")
        n_components = min(30, cells_by_genes.shape[0] - 1, cells_by_genes.shape[1] - 1)
        if n_components < 2:
            raise RuntimeError(f"{row.Study} / {row.DF_partition}: too few cells/genes for Scrublet")
        expected_rate = float(row.expected_doublet_rate)
        caller = scr.Scrublet(cells_by_genes, expected_doublet_rate=expected_rate, random_state=20260913)
        scores, calls = caller.scrub_doublets(
            min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=n_components
        )
        # Scrublet sometimes cannot identify a valley between the simulated
        # doublet and observed-cell score modes.  It returns `threshold_ =
        # None` in that documented situation.  We do not silently accept its
        # all-negative provisional result: use the pre-specified library-size
        # prior for this all-QC-passed technical partition as a reproducible
        # upper-score quantile fallback, and record that choice explicitly.
        # Versions of Scrublet differ: some set threshold_ to None and others
        # omit the attribute entirely when automatic valley detection fails.
        threshold = getattr(caller, "threshold_", None)
        if threshold is None or not np.isfinite(threshold):
            threshold = float(np.quantile(scores, max(0.0, min(1.0, 1.0 - expected_rate))))
            calls = np.asarray(scores) >= threshold
            threshold_source = "expected_doublet_rate_quantile_fallback"
        else:
            threshold = float(threshold)
            calls = np.asarray(calls, dtype=bool)
            threshold_source = "Scrublet_automatic"
        output = pd.DataFrame({
            "Cell": cells["cell_id"].astype(str), "Study": row.Study,
            "DF_partition": row.DF_partition, "Scrublet_score": scores,
            "Scrublet_call": np.where(calls, "Doublet", "Singlet"),
            "Scrublet_threshold": threshold,
            "Scrublet_threshold_source": threshold_source,
            "expected_doublet_rate": expected_rate,
            "n_all_QC_passed_cells_in_partition": int(row.n_cells),
        })
        output.to_csv(output_root / f"Scrublet_{row.safe_partition}.csv", index=False)
        summary.append({"Study": row.Study, "DF_partition": row.DF_partition,
                        "n_all_QC_passed_cells": int(row.n_cells), "n_genes": int(row.n_genes),
                        "expected_doublet_rate": expected_rate, "threshold": threshold,
                        "threshold_source": threshold_source, "n_called_doublet": int(calls.sum())})
        print(f"{row.Study} / {row.DF_partition}: Scrublet complete ({int(calls.sum())}/{int(row.n_cells)} all-QC-passed cells called doublet; {threshold_source})", flush=True)
    pd.DataFrame(summary).to_csv(output_root / "Scrublet_run_summary.csv", index=False)
    with (RUNTIME_DIR / "stage2_scrublet_versions.json").open("w") as handle:
        json.dump({"python": os.sys.version, "scrublet": importlib.metadata.version("scrublet"),
                   "numpy": importlib.metadata.version("numpy"), "scipy": importlib.metadata.version("scipy"),
                   "scikit_learn": importlib.metadata.version("scikit-learn")}, handle, indent=2)


if __name__ == "__main__":
    main()
