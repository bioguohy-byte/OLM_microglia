#!/usr/bin/env python3
"""Recalibrate Scrublet for the historical MG0--MG5 diagnostic analysis.

This is a new, non-destructive analysis branch.  It keeps the original
DoubletFinder/scDblFinder calls and the historical MG0--MG5 labels unchanged.

Why this branch differs from the historical Scrublet run:
  * calls use the narrowest available technical proxy (``Sample``), rather
    than treating an entire study as one capture library;
  * samples with fewer than 300 QC-passed cells are explicitly pooled only for
    Scrublet scoring, with the original sample identity retained;
  * the automatic Scrublet threshold is retained for audit, while the formal
    recalibrated call uses a per-partition threshold at the empirical score
    quantile corresponding to that partition's expected doublet prior.

The sample column is a technical-partition proxy because the input object has
no more granular capture-channel/library field.  A later rerun should replace
it when authentic library/channel metadata become available.
"""

import csv
import importlib.metadata
import json
import math
import os
from pathlib import Path

PROJECT_ROOT = Path(__import__("os").path.join(__import__("os").environ["MG_PROJECT_ROOT"], "Reference"))
SOURCE_ROOT = PROJECT_ROOT / "combine/results/MG0_MG5_doublet_diagnosis"
RESULT_ROOT = PROJECT_ROOT / "combine/results/MG0_MG5_doublet_diagnosis_recalibrated_scrublet"
RUNTIME_DIR = RESULT_ROOT / "runtime"
MIN_PARTITION_CELLS = 300
SEED = 20260908

# Keep every cache beneath the project directory, never in a home directory.
os.environ["PYTHONPYCACHEPREFIX"] = str(RUNTIME_DIR / "python_pycache")
os.environ["MPLCONFIGDIR"] = str(RUNTIME_DIR / "matplotlib")
os.environ["NUMBA_CACHE_DIR"] = str(RUNTIME_DIR / "numba_cache")
for cache in (
    Path(os.environ["PYTHONPYCACHEPREFIX"]),
    Path(os.environ["MPLCONFIGDIR"]),
    Path(os.environ["NUMBA_CACHE_DIR"]),
):
    cache.mkdir(parents=True, exist_ok=True)

import numpy as np
import pandas as pd
import scipy.io
import scrublet as scr


def safe_name(value: str) -> str:
    return "".join(char if char.isalnum() else "_" for char in str(value)).strip("_")


def expected_rate(n_cells: int) -> float:
    """Empirical 10x-style prior retained from the historical diagnostic run."""
    return min(0.20, 0.008 * n_cells / 1000.0)


def nearest_expected_count(rate: float, n_cells: int) -> int:
    """Return the nearest whole-cell expectation without forcing a call for tiny inputs."""
    return int(math.floor(rate * n_cells + 0.5))


def make_partitions(study_meta: pd.DataFrame) -> tuple[dict[str, np.ndarray], pd.DataFrame]:
    """Use Sample as the available technical proxy and pool only undersized samples."""
    if study_meta["Sample"].isna().any() or (study_meta["Sample"].astype(str).str.strip() == "").any():
        raise RuntimeError("Missing Sample values prevent sample-proxy partitioning")
    sample_sizes = study_meta.groupby("Sample", sort=True).size()
    large_samples = sample_sizes[sample_sizes >= MIN_PARTITION_CELLS].index.tolist()
    small_samples = sample_sizes[sample_sizes < MIN_PARTITION_CELLS].index.tolist()
    partitions: dict[str, np.ndarray] = {}
    mapping_rows: list[dict[str, object]] = []

    for sample in large_samples:
        label = f"sample__{safe_name(sample)}"
        indices = np.flatnonzero(study_meta["Sample"].astype(str).to_numpy() == str(sample))
        partitions[label] = indices
        for cell in study_meta.iloc[indices]["Cell"]:
            mapping_rows.append({
                "Cell": cell,
                "Study": study_meta["Study"].iat[0],
                "Sample": sample,
                "Scrublet_partition": label,
                "partition_strategy": "Sample as technical proxy (n >= 300)",
            })

    if small_samples:
        label = "pooled_small_samples"
        indices = np.flatnonzero(study_meta["Sample"].isin(small_samples).to_numpy())
        partitions[label] = indices
        for cell, sample in study_meta.iloc[indices][["Cell", "Sample"]].itertuples(index=False):
            mapping_rows.append({
                "Cell": cell,
                "Study": study_meta["Study"].iat[0],
                "Sample": sample,
                "Scrublet_partition": label,
                "partition_strategy": "Pooled samples with n < 300; original Sample retained",
            })

    if not partitions:
        raise RuntimeError(f"No partitions created for {study_meta['Study'].iat[0]}")
    return partitions, pd.DataFrame(mapping_rows)


def run_partition(counts, cell_rows: pd.DataFrame, partition: str) -> tuple[pd.DataFrame, dict[str, object]]:
    n_cells, n_genes = counts.shape
    if n_cells < 10 or n_genes < 10:
        raise RuntimeError(f"{partition}: too few cells or genes for Scrublet")
    n_components = min(30, n_cells - 1, n_genes - 1)
    rate = expected_rate(n_cells)
    caller = scr.Scrublet(counts, expected_doublet_rate=rate, random_state=SEED)
    scores, auto_calls = caller.scrub_doublets(
        min_counts=2,
        min_cells=3,
        min_gene_variability_pctl=85,
        n_prin_comps=n_components,
    )
    auto_threshold = float(caller.threshold_)
    auto_detected_rate = float(caller.detected_doublet_rate_)
    auto_detectable_fraction = float(caller.detectable_doublet_fraction_)
    auto_overall_rate = float(caller.overall_doublet_rate_)

    target_n = nearest_expected_count(rate, n_cells)
    if target_n == 0:
        recal_threshold = np.nan
        recal_calls = np.zeros(n_cells, dtype=bool)
    else:
        ordered_scores = np.sort(np.asarray(scores))[::-1]
        recal_threshold = float(ordered_scores[target_n - 1])
        recal_calls = caller.call_doublets(threshold=recal_threshold, verbose=True)
        recal_calls = np.asarray(recal_calls, dtype=bool)

    result = cell_rows[["Cell", "Study", "Sample"]].copy()
    result["Scrublet_partition"] = partition
    result["Scrublet_score"] = np.asarray(scores, dtype=float)
    result["Scrublet_auto_threshold"] = auto_threshold
    result["Scrublet_auto_call"] = np.where(np.asarray(auto_calls, dtype=bool), "Doublet", "Singlet")
    result["Scrublet_recalibrated_threshold"] = recal_threshold
    result["Scrublet_recalibrated_call"] = np.where(recal_calls, "Doublet", "Singlet")
    result["expected_doublet_rate"] = rate
    result["expected_doublet_n"] = target_n

    summary = {
        "Study": cell_rows["Study"].iat[0],
        "Scrublet_partition": partition,
        "n_cells": n_cells,
        "n_genes": n_genes,
        "n_prin_comps": n_components,
        "expected_doublet_rate": rate,
        "expected_doublet_n": target_n,
        "auto_threshold": auto_threshold,
        "auto_n_called_doublet": int(np.asarray(auto_calls, dtype=bool).sum()),
        "auto_detected_doublet_rate": auto_detected_rate,
        "auto_detectable_doublet_fraction": auto_detectable_fraction,
        "auto_estimated_overall_doublet_rate": auto_overall_rate,
        "recalibrated_threshold": recal_threshold,
        "recalibrated_n_called_doublet": int(recal_calls.sum()),
        "recalibrated_detected_doublet_rate": float(recal_calls.mean()),
    }
    return result, summary


def main() -> None:
    input_root = SOURCE_ROOT / "scrublet_input"
    source_table = SOURCE_ROOT / "table/stage1_cell_metadata_before_scrublet.csv"
    table_dir = RESULT_ROOT / "table"
    table_dir.mkdir(parents=True, exist_ok=True)
    RUNTIME_DIR.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(source_table, dtype={"Cell": str, "Study": str, "Sample": str})
    required = {"Cell", "Study", "Sample"}
    if missing := required.difference(metadata.columns):
        raise RuntimeError(f"Historical metadata lacks: {', '.join(sorted(missing))}")
    manifest = pd.read_csv(input_root / "scrublet_manifest.csv")

    call_rows: list[pd.DataFrame] = []
    mapping_rows: list[pd.DataFrame] = []
    parameter_rows: list[dict[str, object]] = []

    for manifest_row in manifest.itertuples(index=False):
        study = str(manifest_row.Study)
        study_dir = input_root / str(manifest_row.safe_study)
        counts = scipy.io.mmread(study_dir / "counts_genes_by_cells.mtx").tocsc().T
        cells = pd.read_csv(study_dir / "cells.tsv", sep="\t", dtype={"cell_id": str})
        if counts.shape[0] != len(cells):
            raise RuntimeError(f"{study}: count matrix and cells.tsv disagree")
        study_meta = metadata.loc[metadata["Study"] == study].copy()
        study_meta = cells.rename(columns={"cell_id": "Cell"}).merge(
            study_meta[["Cell", "Study", "Sample"]], on="Cell", how="left", validate="one_to_one"
        )
        if study_meta["Study"].isna().any() or study_meta["Sample"].isna().any():
            raise RuntimeError(f"{study}: could not match all count-matrix cells to historical metadata")
        if not (study_meta["Study"] == study).all():
            raise RuntimeError(f"{study}: cell-to-study mapping mismatch")

        partitions, mapping = make_partitions(study_meta)
        mapping_rows.append(mapping)
        print(f"\n===== {study}: {len(partitions)} sample-proxy partition(s) =====", flush=True)
        for partition, indices in partitions.items():
            partition_meta = study_meta.iloc[indices].reset_index(drop=True)
            partition_counts = counts[indices, :]
            calls, parameters = run_partition(partition_counts, partition_meta, partition)
            strategy = mapping.loc[mapping["Scrublet_partition"] == partition, "partition_strategy"].iat[0]
            parameters["partition_strategy"] = strategy
            parameters["source_samples"] = ";".join(sorted(partition_meta["Sample"].astype(str).unique()))
            call_rows.append(calls)
            parameter_rows.append(parameters)
            print(
                f"{study} / {partition}: recalibrated "
                f"{parameters['recalibrated_n_called_doublet']}/{parameters['n_cells']} doublets "
                f"at score {parameters['recalibrated_threshold']}",
                flush=True,
            )

    calls_all = pd.concat(call_rows, ignore_index=True)
    mapping_all = pd.concat(mapping_rows, ignore_index=True)
    parameters_all = pd.DataFrame(parameter_rows)
    if calls_all["Cell"].duplicated().any() or mapping_all["Cell"].duplicated().any():
        raise RuntimeError("A cell received more than one recalibrated Scrublet partition")
    if set(calls_all["Cell"]) != set(metadata["Cell"]):
        raise RuntimeError("Recalibrated calls do not cover the historical diagnostic cell set exactly")

    calls_all.to_csv(table_dir / "recalibrated_Scrublet_per_cell.csv", index=False)
    mapping_all.to_csv(table_dir / "recalibrated_Scrublet_partition_mapping.csv", index=False)
    parameters_all.to_csv(table_dir / "recalibrated_Scrublet_parameters_by_partition.csv", index=False)
    study_summary = parameters_all.groupby("Study", as_index=False).agg(
        n_partitions=("Scrublet_partition", "count"),
        n_cells=("n_cells", "sum"),
        expected_doublet_n=("expected_doublet_n", "sum"),
        auto_n_called_doublet=("auto_n_called_doublet", "sum"),
        recalibrated_n_called_doublet=("recalibrated_n_called_doublet", "sum"),
    )
    study_summary["recalibrated_doublet_rate"] = (
        study_summary["recalibrated_n_called_doublet"] / study_summary["n_cells"]
    )
    study_summary.to_csv(table_dir / "recalibrated_Scrublet_summary_by_study.csv", index=False)

    manifest_out = {
        "analysis": "Historical MG0--MG5 Scrublet recalibration",
        "source_result_root": str(SOURCE_ROOT),
        "input": "The exact QC-passed historical MG0--MG5 integer-count diagnostic cells",
        "partition_policy": (
            "Sample is the narrowest available technical proxy. Samples with fewer than "
            f"{MIN_PARTITION_CELLS} cells are pooled only for Scrublet scoring."
        ),
        "threshold_policy": (
            "Retain Scrublet automatic threshold for audit; define recalibrated call at the "
            "observed-score quantile that selects the partition's expected doublet count."
        ),
        "expected_rate_formula": "min(0.20, 0.008 * partition_cells / 1000)",
        "seed": SEED,
        "n_cells": int(len(calls_all)),
        "n_recalibrated_doublets": int((calls_all["Scrublet_recalibrated_call"] == "Doublet").sum()),
        "packages": {
            "python": os.sys.version,
            "scrublet": importlib.metadata.version("scrublet"),
            "numpy": importlib.metadata.version("numpy"),
            "scipy": importlib.metadata.version("scipy"),
            "pandas": importlib.metadata.version("pandas"),
            "scikit_learn": importlib.metadata.version("scikit-learn"),
        },
    }
    with (RUNTIME_DIR / "recalibrated_Scrublet_manifest.json").open("w") as handle:
        json.dump(manifest_out, handle, indent=2)
    print("Recalibrated Scrublet calls complete. Run the stage-3 R plotting script next.", flush=True)


if __name__ == "__main__":
    main()
