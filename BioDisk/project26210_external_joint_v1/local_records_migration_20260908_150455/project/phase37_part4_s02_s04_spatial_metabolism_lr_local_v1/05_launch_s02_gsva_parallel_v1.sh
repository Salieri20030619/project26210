#!/usr/bin/env bash
set -euo pipefail
root="$1"
rscript="$2"
manifest="$3"
custom="$4"
kegg="$5"
mkdir -p "$root/s02_gsva_sections_v1" "$root/logs"
samples=(GLMF2 H46_T J38_T R114_T R29_T R51_T S15_T X49_T X98_T Y12_T Y27_T Y7_T Z43_T)
pids=()
for sample in "${samples[@]}"; do
  "$rscript" --vanilla "$root/provenance/04_run_s02_gsva_section_v1.R" \
    --manifest "$manifest" --sample "$sample" --custom-gmt "$custom" --kegg-gmt "$kegg" \
    --output "$root/s02_gsva_sections_v1/$sample" > "$root/logs/s02_gsva_${sample}.log" 2>&1 &
  pids+=("$!")
done
failed=0
for i in "${!pids[@]}"; do
  if ! wait "${pids[$i]}"; then
    echo "FAIL ${samples[$i]}" >&2; failed=1
  else
    echo "PASS ${samples[$i]}"
  fi
done
exit "$failed"
