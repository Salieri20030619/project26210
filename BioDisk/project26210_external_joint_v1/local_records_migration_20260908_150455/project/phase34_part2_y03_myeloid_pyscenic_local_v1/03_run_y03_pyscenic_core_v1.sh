#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 3 ]]; then
  echo "usage: $0 WORK_ROOT PROJECT_ROOT THREADS" >&2
  exit 2
fi
WORK_ROOT="$1"
PROJECT_ROOT="$2"
THREADS="$3"
PY_SCENIC="$PROJECT_ROOT/.venv-pyscenic/bin/python"
WRAPPER="$PROJECT_ROOT/run_pyscenic_compat.py"
REF="$PROJECT_ROOT/references/pyscenic/hg38"
DB500="$REF/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
DB10K="$REF/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
ANNOT="$REF/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
TFS="$REF/allTFs_hg38.txt"
LOOMS="$WORK_ROOT/looms_v1"
OUT="$WORK_ROOT/pyscenic_core_v1"

if [[ -e "$OUT" ]]; then
  echo "Refusing to overwrite $OUT" >&2
  exit 1
fi
mkdir -p "$OUT"
"$PY_SCENIC" "$WRAPPER" grn "$LOOMS/balanced_discovery_counts_v1.loom" "$TFS" \
  --method grnboost2 --num_workers "$THREADS" --seed 20260814 --cell_id_attribute CellID --gene_attribute Gene --sparse \
  -o "$OUT/adjacencies_grnboost2_v1.tsv"
"$PY_SCENIC" "$WRAPPER" ctx "$OUT/adjacencies_grnboost2_v1.tsv" "$DB500" "$DB10K" \
  --annotations_fname "$ANNOT" --expression_mtx_fname "$LOOMS/balanced_discovery_counts_v1.loom" \
  --num_workers "$THREADS" --cell_id_attribute CellID --gene_attribute Gene \
  -o "$OUT/regulons_grnboost2_v1.dat"
"$PY_SCENIC" "$WRAPPER" aucell "$LOOMS/all_myeloid_filtered_counts_v1.loom" "$OUT/regulons_grnboost2_v1.dat" \
  --num_workers "$THREADS" --seed 20260814 --cell_id_attribute CellID --gene_attribute Gene \
  -o "$OUT/auc_grnboost2_all_myeloid_v1.csv"
date --iso-8601=seconds > "$OUT/RUN_COMPLETE"
echo Y03_PYSCENIC_CORE_PASS
