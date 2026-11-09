#!/usr/bin/env bash
set -euo pipefail
task=$(cd "$(dirname "$0")/.." && pwd)
combine=${MG_PROJECT_ROOT}/Reference/combine
export R_LIBS_USER=${MG_WORKSPACE_ROOT}/R_libs
export TMPDIR="$task/runtime/tmp"
export XDG_CACHE_HOME="$task/runtime/cache"
export FONTCONFIG_FILE="$task/runtime/fontconfig.xml"
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
mkdir -p "$TMPDIR" "$XDG_CACHE_HOME" "$task/runtime/fonts"
cat > "$FONTCONFIG_FILE" <<EOF
<?xml version="1.0"?>
<!DOCTYPE fontconfig SYSTEM "fonts.dtd">
<fontconfig>
  <dir>$task/runtime/fonts</dir>
  <cachedir>$task/runtime/cache/fontconfig</cachedir>
  <alias><family>sans-serif</family><prefer><family>Arial</family></prefer></alias>
</fontconfig>
EOF
fc-cache -f
fc-match Arial > "$task/audit/fontconfig_match.txt"
if [ "${1:-}" != "--reuse-inputs" ]; then
  /usr/bin/python3 "$task/code/prepare_inputs.py" "$task" "$combine"
fi
/usr/bin/Rscript --vanilla "$task/code/plot_uniform_umaps.R" "$task"
/usr/bin/python3 "$task/code/verify_and_render.py" "$task"
