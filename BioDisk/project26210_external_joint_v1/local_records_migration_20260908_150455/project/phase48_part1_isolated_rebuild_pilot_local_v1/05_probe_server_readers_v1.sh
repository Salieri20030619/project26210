#!/usr/bin/env bash
set -u

echo "probe_started=$(date --iso-8601=seconds)"
echo "system_R=$(command -v Rscript || true)"
if command -v Rscript >/dev/null 2>&1; then
  Rscript --vanilla - <<'RS'
packages <- c("Matrix", "rhdf5", "hdf5r", "data.table", "Seurat", "SeuratObject")
for (pkg in packages) {
  present <- requireNamespace(pkg, quietly = TRUE)
  version <- if (present) as.character(packageVersion(pkg)) else "NA"
  cat(pkg, present, version, sep = "\t")
  cat("\n")
}
RS
fi

env_root=${1:?Usage: 05_probe_server_readers_v1.sh CONDA_ENVS_ROOT}
if [[ -d "$env_root" ]]; then
  while IFS= read -r env_name; do
    py="$env_root/$env_name/bin/python"
    [[ -x "$py" ]] || continue
    if timeout 10 "$py" -c 'import h5py,numpy,scipy; print("PASS\t"+h5py.__version__+"\t"+numpy.__version__+"\t"+scipy.__version__)' 2>/dev/null; then
      echo "python_env=$env_name"
    fi
  done < <(find "$env_root" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)
fi
echo "probe_finished=$(date --iso-8601=seconds)"
