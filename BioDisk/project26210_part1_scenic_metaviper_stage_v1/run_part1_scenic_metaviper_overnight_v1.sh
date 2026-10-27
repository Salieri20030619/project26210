#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
ROOT="$PROJECT/results/part1_rebuild_v1"
SCENIC="$ROOT/07_pyscenic_v1"
VIPER="$ROOT/08_metaviper_v1"
STATUS="$ROOT/07_pyscenic_v1/overnight_status_v1.tsv"
PY_BATCH=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
PY_SCENIC="$PROJECT/.venv-pyscenic/bin/python"
R_SCRIPT=/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript
R_BIN=/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/R
R_LIB="$ROOT/08_metaviper_v1/r_library"
REF="$PROJECT/references/pyscenic/hg38"
WRAPPER="$PROJECT/run_pyscenic_compat.py"
DB500="$REF/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
DB10K="$REF/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
ANNOT="$REF/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
TFS="$REF/allTFs_hg38.txt"
LOGDIR="$SCENIC/logs"

mkdir -p "$LOGDIR" "$SCENIC/intermediate" "$SCENIC/tables" "$VIPER/logs" "$R_LIB"
printf 'time\tstage\tstatus\tdetail\n' > "$STATUS"
record() { printf '%s\t%s\t%s\t%s\n' "$(date -Is)" "$1" "$2" "$3" >> "$STATUS"; }
gate() { "$PROJECT/check_memory_gate.sh" "$1"; }
run_stage() {
  local stage="$1"; shift
  record "$stage" STARTED "$*"
  gate "$stage"
  "$@"
  record "$stage" PASS completed
}

record overnight STARTED "recommended balanced pySCENIC plus metaVIPER"

if [[ ! -f "$SCENIC/input/RUN_COMPLETE" ]]; then
  run_stage prepare_input "$PY_BATCH" "$SCENIC/prepare_part1_scenic_metaviper_input_v1.py"
else
  record prepare_input REUSED RUN_COMPLETE
fi

if [[ ! -f "$SCENIC/input/LOOMS_COMPLETE" ]]; then
  run_stage build_looms "$PY_SCENIC" "$SCENIC/build_part1_pyscenic_looms_v1.py"
else
  record build_looms REUSED LOOMS_COMPLETE
fi

pushd "$REF" >/dev/null
sha1sum -c "$(basename "$DB500.sha1sum.txt")" | tee "$LOGDIR/cistarget_db500_sha1_v1.log"
sha1sum -c "$(basename "$DB10K.sha1sum.txt")" | tee "$LOGDIR/cistarget_db10k_sha1_v1.log"
popd >/dev/null

if [[ ! -s "$SCENIC/intermediate/adjacencies_genie3_v1.tsv" ]]; then
  run_stage pyscenic_genie3 "$PY_SCENIC" "$WRAPPER" grn "$SCENIC/input/focused_genie3_counts_v1.loom" "$TFS" --method genie3 --num_workers 12 --seed 26210 --cell_id_attribute CellID --gene_attribute Gene --sparse -o "$SCENIC/intermediate/adjacencies_genie3_v1.tsv"
else record pyscenic_genie3 REUSED adjacency; fi

if [[ ! -s "$SCENIC/intermediate/regulons_genie3_v1.dat" ]]; then
  run_stage pyscenic_ctx_genie3 "$PY_SCENIC" "$WRAPPER" ctx "$SCENIC/intermediate/adjacencies_genie3_v1.tsv" "$DB500" "$DB10K" --annotations_fname "$ANNOT" --expression_mtx_fname "$SCENIC/input/focused_genie3_counts_v1.loom" --num_workers 12 --cell_id_attribute CellID --gene_attribute Gene -o "$SCENIC/intermediate/regulons_genie3_v1.dat"
else record pyscenic_ctx_genie3 REUSED regulons; fi

if [[ ! -s "$SCENIC/intermediate/auc_genie3_all_cells_v1.csv" ]]; then
  run_stage pyscenic_aucell_genie3 "$PY_SCENIC" "$WRAPPER" aucell "$SCENIC/input/all_cells_filtered_counts_v1.loom" "$SCENIC/intermediate/regulons_genie3_v1.dat" --num_workers 12 --seed 26210 --cell_id_attribute CellID --gene_attribute Gene -o "$SCENIC/intermediate/auc_genie3_all_cells_v1.csv"
else record pyscenic_aucell_genie3 REUSED auc; fi

if [[ ! -s "$SCENIC/intermediate/adjacencies_grnboost2_v1.tsv" ]]; then
  run_stage pyscenic_grnboost2 "$PY_SCENIC" "$WRAPPER" grn "$SCENIC/input/balanced_grn_counts_v1.loom" "$TFS" --method grnboost2 --num_workers 12 --seed 26210 --cell_id_attribute CellID --gene_attribute Gene --sparse -o "$SCENIC/intermediate/adjacencies_grnboost2_v1.tsv"
else record pyscenic_grnboost2 REUSED adjacency; fi

if [[ ! -s "$SCENIC/intermediate/regulons_grnboost2_v1.dat" ]]; then
  run_stage pyscenic_ctx_grnboost2 "$PY_SCENIC" "$WRAPPER" ctx "$SCENIC/intermediate/adjacencies_grnboost2_v1.tsv" "$DB500" "$DB10K" --annotations_fname "$ANNOT" --expression_mtx_fname "$SCENIC/input/balanced_grn_counts_v1.loom" --num_workers 12 --cell_id_attribute CellID --gene_attribute Gene -o "$SCENIC/intermediate/regulons_grnboost2_v1.dat"
else record pyscenic_ctx_grnboost2 REUSED regulons; fi

if [[ ! -s "$SCENIC/intermediate/auc_grnboost2_all_cells_v1.csv" ]]; then
  run_stage pyscenic_aucell_grnboost2 "$PY_SCENIC" "$WRAPPER" aucell "$SCENIC/input/all_cells_filtered_counts_v1.loom" "$SCENIC/intermediate/regulons_grnboost2_v1.dat" --num_workers 12 --seed 26210 --cell_id_attribute CellID --gene_attribute Gene -o "$SCENIC/intermediate/auc_grnboost2_all_cells_v1.csv"
else record pyscenic_aucell_grnboost2 REUSED auc; fi

touch "$SCENIC/PYSCENIC_CORE_COMPLETE"

if [[ ! -d "$R_LIB/aracne.networks" ]]; then
  run_stage install_aracne_local "$R_BIN" CMD INSTALL --library="$R_LIB" --no-docs --no-html "$PROJECT/references/viper/aracne.networks_1.36.0.tar.gz"
else record install_aracne_local REUSED project_library; fi

if [[ ! -f "$VIPER/RUN_COMPLETE" ]]; then
  record metaviper STARTED patient_stratified_KIRC_KIRP
  gate metaviper_v1
  R_LIBS_USER="$R_LIB" "$R_SCRIPT" "$VIPER/run_part1_metaviper_v1.R" 2>&1 | tee "$VIPER/logs/run_part1_metaviper_v1.log"
  record metaviper PASS completed
else record metaviper REUSED RUN_COMPLETE; fi

gate overnight_terminal_v1 | tee "$SCENIC/logs/terminal_memory_gate_v1.log"
record overnight PASS all_core_stages_complete
touch "$SCENIC/OVERNIGHT_COMPLETE"
echo PART1_SCENIC_METAVIPER_OVERNIGHT_V1_PASS
