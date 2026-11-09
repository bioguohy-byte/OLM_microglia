#!/usr/bin/env bash
set -euo pipefail
task_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
case "$task_dir" in ${MG_WORKSPACE_ROOT}/*) ;; *) echo "Run in the server24 task directory under ${MG_WORKSPACE_ROOT}" >&2; exit 1;; esac
mkdir -p "$task_dir/runtime/tmp" "$task_dir/runtime/cache" "$task_dir/runtime/config" "$task_dir/runtime/data"
export TMPDIR="$task_dir/runtime/tmp"
export XDG_CACHE_HOME="$task_dir/runtime/cache"
export XDG_CONFIG_HOME="$task_dir/runtime/config"
export XDG_DATA_HOME="$task_dir/runtime/data"
export R_LIBS_USER=${MG_WORKSPACE_ROOT}/R_libs
export PIP_CACHE_DIR="$task_dir/runtime/cache/pip"
export UV_CACHE_DIR="$task_dir/runtime/cache/uv"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$task_dir"
/usr/bin/Rscript --vanilla "$task_dir/plot_celltype_labels_on_umap.R" "$task_dir" > "$task_dir/runtime/plotting.log" 2>&1
