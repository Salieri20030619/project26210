#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase64_large_immune_validation_v1
INPUT="$PROJECT/inputs_v1/GSE207493/tf_regulon_de_novo_v1/scenic_inputs_v1"
OUT="$PROJECT/results_v1/GSE207493/tf_regulon_de_novo_v1/scenic_de_novo_v1"
LOG="$PROJECT/logs_v1/gse207493_de_novo_scenic_v1.log"
PY=/mnt/BioDisk/ludingyu/project26210/.venv-pyscenic/bin/python
WRAPPER=/mnt/BioDisk/ludingyu/project26210/run_pyscenic_compat.py
REF=/mnt/BioDisk/ludingyu/project26210/references/pyscenic/hg38
TFS="$REF/allTFs_hg38.txt"
DB500="$REF/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
DB10K="$REF/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
ANNOT="$REF/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

mkdir -p "$OUT"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMEXPR_NUM_THREADS=8
printf 'time\tlineage\tstage\tstatus\n' > "$OUT/run_status_v1.tsv"
record() { printf '%s\t%s\t%s\t%s\n' "$(date -Is)" "$1" "$2" "$3" >> "$OUT/run_status_v1.tsv"; }

for lineage in tumor myeloid_dc treg; do
  in_dir="$INPUT/$lineage"
  out_dir="$OUT/$lineage"
  mkdir -p "$out_dir"
  adjacency="$out_dir/adjacencies_grnboost2_de_novo_v1.tsv"
  regulons="$out_dir/regulons_grnboost2_de_novo_v1.dat"
  auc="$out_dir/auc_grnboost2_all_cells_v1.csv"

  if [[ ! -s "$adjacency" ]]; then
    record "$lineage" grn STARTED
    "$PY" "$WRAPPER" grn "$in_dir/balanced_grn_counts_v1.loom" "$TFS" \
      --method grnboost2 --num_workers 8 --seed 26210 \
      --cell_id_attribute CellID --gene_attribute Gene --sparse -o "$adjacency"
    record "$lineage" grn PASS
  else
    record "$lineage" grn REUSED
  fi

  if [[ ! -s "$regulons" ]]; then
    record "$lineage" ctx STARTED
    "$PY" "$WRAPPER" ctx "$adjacency" "$DB500" "$DB10K" \
      --annotations_fname "$ANNOT" --expression_mtx_fname "$in_dir/balanced_grn_counts_v1.loom" \
      --num_workers 8 --cell_id_attribute CellID --gene_attribute Gene -o "$regulons"
    record "$lineage" ctx PASS
  else
    record "$lineage" ctx REUSED
  fi

  if [[ ! -s "$auc" ]]; then
    record "$lineage" aucell STARTED
    "$PY" "$WRAPPER" aucell "$in_dir/all_cells_filtered_counts_v1.loom" "$regulons" \
      --num_workers 8 --seed 26210 --cell_id_attribute CellID --gene_attribute Gene -o "$auc"
    record "$lineage" aucell PASS
  else
    record "$lineage" aucell REUSED
  fi
done

date -Is > "$OUT/RUN_COMPLETE"
echo DE_NOVO_SCENIC_V1_PASS
