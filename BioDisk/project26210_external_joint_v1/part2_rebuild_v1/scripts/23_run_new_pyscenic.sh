#!/usr/bin/env bash
set -euo pipefail
R=/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1
Q=/mnt/BioDisk/ludingyu/project26210
O="$R/08_regulons/network_v1"
test ! -e "$O"
mkdir "$O"
PY="$Q/.venv-pyscenic/bin/python"
W="$Q/run_pyscenic_compat.py"
REF="$Q/references/pyscenic/hg38"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$PY" "$W" grn "$R/08_regulons/looms_v1/discovery.loom" "$REF/allTFs_hg38.txt" --method grnboost2 --num_workers 8 --seed 262148 --cell_id_attribute CellID --gene_attribute Gene --sparse -o "$O/adjacencies.tsv"
"$PY" "$W" ctx "$O/adjacencies.tsv" "$REF/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" "$REF/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" --annotations_fname "$REF/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl" --expression_mtx_fname "$R/08_regulons/looms_v1/discovery.loom" --num_workers 8 --cell_id_attribute CellID --gene_attribute Gene -o "$O/regulons.dat"
"$PY" "$W" aucell "$R/08_regulons/looms_v1/all_cells.loom" "$O/regulons.dat" --num_workers 8 --seed 262148 --cell_id_attribute CellID --gene_attribute Gene -o "$O/AUC_ALL_CELLS.csv"
date --iso-8601=seconds > "$O/RUN_COMPLETE"
