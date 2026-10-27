#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 11 ]]; then
  echo "usage: $0 PYSCENIC_PY WRAPPER PREP_DIR GENIE3_LOOM GRNBOOST2_LOOM ALL_CELLS_LOOM DB500 DB10K ANNOT INSPECT_SCRIPT OUTPUT_DIR" >&2
  exit 2
fi

PY="$1"; WRAPPER="$2"; PREP="$3"; GLOOM="$4"; BLOOM="$5"; ALLLOOM="$6"
DB500="$7"; DB10K="$8"; ANNOT="$9"; INSPECT="${10}"; OUTPUT="${11}"
PART="${OUTPUT}.part"
[[ ! -e "$OUTPUT" && ! -e "$PART" ]] || { echo "refusing existing output" >&2; exit 3; }
mkdir -p "$PART/logs" "$PART/regulons" "$PART/auc"
printf 'method\tvariant\tnes_threshold\tmin_genes\tctx_status\thif1a_positive_present\tauc_status\n' > "$PART/run_summary_v1.tsv"

run_one() {
  local method="$1" adjacency="$2" loom="$3" variant="$4" nes="$5" min_genes="$6"
  local stem="${method}_${variant}"
  set +e
  "$PY" "$WRAPPER" ctx "$adjacency" "$DB500" "$DB10K" \
    --annotations_fname "$ANNOT" --expression_mtx_fname "$loom" --num_workers 12 \
    --cell_id_attribute CellID --gene_attribute Gene \
    --rank_threshold 5000 --auc_threshold 0.05 --nes_threshold "$nes" \
    --thresholds 0.75 0.90 --top_n_targets 50 --top_n_regulators 5 10 50 \
    --min_genes "$min_genes" -o "$PART/regulons/${stem}.dat" \
    > "$PART/logs/${stem}_ctx.log" 2>&1
  local ctx_rc=$?
  set -e
  if [[ "$ctx_rc" -ne 0 ]]; then
    if [[ ! -s "$PART/regulons/${stem}.dat" ]] && grep -q 'Signatures dataframe is empty' "$PART/logs/${stem}_ctx.log"; then
      printf 'regulon\tn_targets\tcontext\n' > "$PART/regulons/${stem}_summary.tsv"
      printf 'regulon\ttarget\tweight\n' > "$PART/regulons/${stem}_targets.tsv"
      printf '%s\t%s\t%s\t%s\tPASS_NO_REGULON\tfalse\tNOT_RUN_NO_HIF1A_REGULON\n' "$method" "$variant" "$nes" "$min_genes" >> "$PART/run_summary_v1.tsv"
      return 0
    fi
    echo "unexpected ctx failure for ${stem}; rc=${ctx_rc}" >&2
    return "$ctx_rc"
  fi
  "$PY" "$INSPECT" --dat "$PART/regulons/${stem}.dat" \
    --summary "$PART/regulons/${stem}_summary.tsv" \
    --targets "$PART/regulons/${stem}_targets.tsv" \
    --present-marker "$PART/regulons/${stem}_HIF1A_PRESENT.json" \
    > "$PART/logs/${stem}_inspect.log" 2>&1
  local present="false" auc_status="NOT_RUN_NO_HIF1A_REGULON"
  if [[ -f "$PART/regulons/${stem}_HIF1A_PRESENT.json" ]]; then
    present="true"
    "$PY" "$WRAPPER" aucell "$ALLLOOM" "$PART/regulons/${stem}.dat" \
      --num_workers 12 --seed 26210 --cell_id_attribute CellID --gene_attribute Gene \
      -o "$PART/auc/${stem}_auc.csv" > "$PART/logs/${stem}_aucell.log" 2>&1
    auc_status="PASS"
  fi
  printf '%s\t%s\t%s\t%s\tPASS_REGULON\t%s\t%s\n' "$method" "$variant" "$nes" "$min_genes" "$present" "$auc_status" >> "$PART/run_summary_v1.tsv"
}

for spec in 'default 3.0 20' 'min10_only 3.0 10' 'nes2_5_only 2.5 20' 'combined_light 2.5 10'; do
  read -r variant nes min_genes <<< "$spec"
  run_one GENIE3 "$PREP/genie3_hif1a_adjacencies_v1.tsv" "$GLOOM" "$variant" "$nes" "$min_genes"
  run_one GRNBoost2 "$PREP/grnboost2_hif1a_adjacencies_v1.tsv" "$BLOOM" "$variant" "$nes" "$min_genes"
done

touch "$PART/RUN_COMPLETE"
mv "$PART" "$OUTPUT"
echo PHASE28_HIF1A_CTX_PASS
