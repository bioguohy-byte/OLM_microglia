#!/usr/bin/env python3
"""Artificial-doublet visualisation, part 2 of 3: Scrublet simulation only.

This program calls Scrublet.simulate_doublets() on the project-exported raw
integer count matrix.  It does not recalculate Scrublet scores or alter the
already completed real-cell Scrublet calls.
"""

import importlib.metadata
import json
import os
from pathlib import Path

PROJECT_ROOT = Path(__import__("os").path.join(__import__("os").environ["MG_PROJECT_ROOT"], "Reference"))
RESULT_ROOT = PROJECT_ROOT / "combine/results/artificial_doublet_visualization"
INPUT_ROOT = RESULT_ROOT / "input"
OUTPUT_ROOT = RESULT_ROOT / "scrublet_output"
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
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(INPUT_ROOT / "partition_manifest.csv")
    summary_rows: list[dict[str, object]] = []
    for row_index, row in manifest.reset_index(drop=True).iterrows():
        partition_dir = INPUT_ROOT / str(row["safe_partition"])
        input_dir = partition_dir / "scrublet_input"
        counts_cells_by_genes = scipy.io.mmread(input_dir / "counts_genes_by_cells.mtx").tocsc().T
        cell_ids = pd.read_csv(input_dir / "cells.tsv", sep="\t")["cell_id"].astype(str).to_numpy()
        gene_ids = pd.read_csv(input_dir / "genes.tsv", sep="\t")["gene_id"].astype(str).to_numpy()
        if counts_cells_by_genes.shape != (len(cell_ids), len(gene_ids)):
            raise RuntimeError(f"{row.DF_partition}: input matrix dimensions do not match cell/gene manifests")

        # The ratio equals the artificial/real ratio implied by DoubletFinder's
        # pN=0.25 call and keeps every artificial class at an interpretable,
        # manageable scale. Scrublet's own method supplies the parent indices.
        sim_ratio = float(row["Scrublet_sim_doublet_ratio"])
        scrub = scr.Scrublet(
            counts_cells_by_genes,
            sim_doublet_ratio=sim_ratio,
            random_state=2026092800 + int(row_index),
        )
        scrub.simulate_doublets(sim_doublet_ratio=sim_ratio, synthetic_doublet_umi_subsampling=1.0)
        simulated_cells_by_genes = scrub._E_sim.tocsc()
        parent_indices = np.asarray(scrub.doublet_parents_, dtype=int)
        if simulated_cells_by_genes.shape[0] != parent_indices.shape[0]:
            raise RuntimeError(f"{row.DF_partition}: Scrublet artificial matrix and parent table disagree")

        artificial_ids = np.asarray([
            f"SCRUBLETart__{row.safe_partition}__{i + 1}" for i in range(simulated_cells_by_genes.shape[0])
        ])
        out_dir = OUTPUT_ROOT / str(row["safe_partition"])
        out_dir.mkdir(parents=True, exist_ok=True)
        scipy.io.mmwrite(out_dir / "Scrublet_artificial_counts_genes_by_cells.mtx", simulated_cells_by_genes.T)
        pd.DataFrame({"cell_id": artificial_ids}).to_csv(out_dir / "Scrublet_artificial_cells.tsv", sep="\t", index=False)
        pd.DataFrame({"gene_id": gene_ids}).to_csv(out_dir / "Scrublet_artificial_genes.tsv", sep="\t", index=False)
        parents = pd.DataFrame({
            "artificial_cell_id": artificial_ids,
            "Study": str(row["Study"]),
            "DF_partition": str(row["DF_partition"]),
            "parent_1": cell_ids[parent_indices[:, 0]],
            "parent_2": cell_ids[parent_indices[:, 1]],
            "construction": "Scrublet.simulate_doublets: raw parent-1 + raw parent-2",
        })
        parents.to_csv(out_dir / "Scrublet_artificial_parent_pairs.csv", index=False)
        summary_rows.append({
            "Study": str(row["Study"]),
            "DF_partition": str(row["DF_partition"]),
            "safe_partition": str(row["safe_partition"]),
            "n_real_cells": int(counts_cells_by_genes.shape[0]),
            "n_detected_features": int(counts_cells_by_genes.shape[1]),
            "sim_doublet_ratio": sim_ratio,
            "n_Scrublet_artificial": int(simulated_cells_by_genes.shape[0]),
        })
        print(
            f"{row.Study} / {row.DF_partition}: simulated {simulated_cells_by_genes.shape[0]} Scrublet artificial profiles "
            f"from {counts_cells_by_genes.shape[0]} real cells",
            flush=True,
        )

    pd.DataFrame(summary_rows).to_csv(OUTPUT_ROOT / "Scrublet_artificial_profile_summary.csv", index=False)
    with (RUNTIME_DIR / "stage2_scrublet_versions.json").open("w") as handle:
        json.dump({
            "python": os.sys.version,
            "scrublet": importlib.metadata.version("scrublet"),
            "numpy": importlib.metadata.version("numpy"),
            "scipy": importlib.metadata.version("scipy"),
        }, handle, indent=2)


if __name__ == "__main__":
    main()
