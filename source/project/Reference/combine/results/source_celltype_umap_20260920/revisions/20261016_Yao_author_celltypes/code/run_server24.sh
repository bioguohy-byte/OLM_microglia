#!/usr/bin/env bash
set -euo pipefail
task_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
case "$task_dir" in ${MG_WORKSPACE_ROOT}/*) ;; *) echo "Run under ${MG_WORKSPACE_ROOT} on server24" >&2; exit 1;; esac
export R_LIBS_USER=${MG_WORKSPACE_ROOT}/R_libs
export TMPDIR="$task_dir/runtime/tmp"
export XDG_CACHE_HOME="$task_dir/runtime/cache"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
mkdir -p "$TMPDIR" "$XDG_CACHE_HOME"
if [[ ! -f "$task_dir/inputs/author_metadata_umaps.rds" ]]; then
  /usr/bin/Rscript --vanilla "$task_dir/code/extract_author_metadata.R" "$task_dir" > "$task_dir/runtime/extract_author_metadata.log" 2>&1
fi
/usr/bin/Rscript --vanilla "$task_dir/code/prepare_author_labels.R" "$task_dir" > "$task_dir/runtime/prepare_author_labels.log" 2>&1
/usr/bin/Rscript --vanilla "$task_dir/code/plot_author_labels.R" "$task_dir" > "$task_dir/runtime/plot_author_labels.log" 2>&1
baseline=${MG_PROJECT_ROOT}/Reference/combine/results/source_celltype_umap_20260920/revisions/20261016_celltype_labels_on_umap_no_boxes
pdfseparate -f 1 -l 5 "$baseline/Six_study_celltypes_labeled_on_UMAPs.pdf" "$task_dir/runtime/base-main-%d.pdf"
pdfunite "$task_dir/runtime/base-main-1.pdf" "$task_dir/runtime/base-main-2.pdf" "$task_dir/runtime/base-main-3.pdf" "$task_dir/runtime/base-main-4.pdf" "$task_dir/runtime/base-main-5.pdf" "$task_dir/output/pdf/Yao_current_UMAP_author_labels.pdf" "$task_dir/output/pdf/Six_study_celltypes_with_Yao_author_labels.pdf"
echo AUTHOR_LABEL_FIGURES_COMPLETE
