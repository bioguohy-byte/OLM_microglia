#!/usr/bin/env python3
"""Corrected full-cell artificial validation, stage 2: Scrublet simulation.

The input to every Scrublet simulation is the complete QC-passed real-cell
partition, not an MG0/MG00 or historical-Microglia subset.
"""

import importlib.metadata
import json
import os
from pathlib import Path

PROJECT_ROOT = Path(__import__("os").path.join(__import__("os").environ["MG_PROJECT_ROOT"], "Reference"))
RESULT_ROOT = PROJECT_ROOT / "combine/results/fullcell_artificial_doublet_validation_v2"
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
    input_root = RESULT_ROOT / "input"
    output_root = RESULT_ROOT / "scrublet_output"
    output_root.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(input_root / "partition_manifest.csv")
    summary = []
    for row in manifest.itertuples(index=False):
        source_dir = input_root / row.safe_partition / "scrublet_input"
        cells_by_genes = scipy.io.mmread(source_dir / "counts_genes_by_cells.mtx").tocsc().T
        cells = pd.read_csv(source_dir / "cells.tsv", sep="\t")
        genes = pd.read_csv(source_dir / "genes.tsv", sep="\t")
        if cells_by_genes.shape != (len(cells), len(genes)) or len(cells) != int(row.n_all_QC_passed_real_cells):
            raise RuntimeError(f"{row.Study} / {row.DF_partition}: full-cell Scrublet input does not match the manifest")
        caller = scr.Scrublet(cells_by_genes, random_state=20260913)
        caller.simulate_doublets(sim_doublet_ratio=float(row.Scrublet_sim_doublet_ratio), synthetic_doublet_umi_subsampling=1.0)
        simulated = caller._E_sim
        # Scrublet's public simulate_doublets() stores sampled observed-cell
        # indices as `doublet_parents_` in the installed release.  Do not
        # substitute a reconstructed random pairing: the exported audit must
        # exactly match the method-generated artificial profile.
        parents = np.asarray(caller.doublet_parents_)
        if simulated.shape[0] != parents.shape[0]:
            raise RuntimeError(f"{row.Study} / {row.DF_partition}: Scrublet artificial profiles and parent index disagree")
        identifiers = np.asarray([f"Scrubart__{row.safe_partition}__{i + 1}" for i in range(simulated.shape[0])])
        destination = output_root / row.safe_partition
        destination.mkdir(parents=True, exist_ok=True)
        scipy.io.mmwrite(destination / "Scrublet_artificial_counts_genes_by_cells.mtx", simulated.T)
        pd.DataFrame({"cell_id": identifiers}).to_csv(destination / "Scrublet_artificial_cells.tsv", sep="\t", index=False)
        pd.DataFrame({"gene_id": genes["gene_id"].astype(str)}).to_csv(destination / "Scrublet_artificial_genes.tsv", sep="\t", index=False)
        parent_table = pd.DataFrame({
            "artificial_cell_id": identifiers, "Study": row.Study, "DF_partition": row.DF_partition,
            "parent_1": cells.iloc[parents[:, 0], 0].astype(str).values,
            "parent_2": cells.iloc[parents[:, 1], 0].astype(str).values,
            "parent_universe": "all QC-passed real cells in this study/technical partition",
            "construction": "Scrublet.simulate_doublets: raw parent-1 + raw parent-2",
        })
        parent_table.to_csv(destination / "Scrublet_artificial_parent_pairs.csv", index=False)
        summary.append({"Study": row.Study, "DF_partition": row.DF_partition,
                        "n_all_QC_passed_real_cells": int(row.n_all_QC_passed_real_cells),
                        "n_Scrublet_artificial": int(simulated.shape[0]),
                        "sim_doublet_ratio": float(row.Scrublet_sim_doublet_ratio)})
        print(f"{row.Study} / {row.DF_partition}: simulated {simulated.shape[0]} Scrublet profiles from all QC-passed real cells", flush=True)
    pd.DataFrame(summary).to_csv(output_root / "Scrublet_artificial_profile_summary.csv", index=False)
    with (RUNTIME_DIR / "stage2_scrublet_versions.json").open("w") as handle:
        json.dump({"python": os.sys.version, "scrublet": importlib.metadata.version("scrublet"),
                   "numpy": importlib.metadata.version("numpy"), "scipy": importlib.metadata.version("scipy"),
                   "scikit_learn": importlib.metadata.version("scikit-learn")}, handle, indent=2)


if __name__ == "__main__":
    main()
