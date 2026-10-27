#!/usr/bin/env bash
set -uo pipefail

ROOT=/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1
STAGE="$ROOT/phase11_external_cnv_consensus_v1"
MAMBA=/mnt/BioDisk/wangmingyuan/miniconda3/bin/mamba
ENV_PREFIX=/mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2
VALIDATOR="$STAGE/scripts/02_validate_terminal_external_cnv_rds_v1.R"
RUN="$STAGE/results/rds_validation_runs/20260812_final_v5"
STATUS="$STAGE/results/FINAL_BACKGROUND_CHAIN_STATUS_v3.tsv"
PROVENANCE="$STAGE/results/infercnv_conda_v2_repair_v3_provenance"

if [[ ! -x "$ENV_PREFIX/bin/Rscript" ]]; then
  echo "Existing isolated R environment is unavailable: $ENV_PREFIX" >&2
  exit 2
fi
if [[ -e "$RUN" ]]; then
  echo "Refusing existing validation run: $RUN" >&2
  exit 2
fi
if [[ -e "$STATUS" ]]; then
  echo "Refusing existing chain status: $STATUS" >&2
  exit 2
fi

status_part="${STATUS}.part"
printf 'stage\tresult\ttimestamp\tdetail\n' > "$status_part"
record() {
  printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$(date -Is)" "$3" >> "$status_part"
}
finish_fail() {
  local stage="$1" code="$2" detail="$3"
  record "$stage" FAIL "$detail"
  record OVERALL FAIL "failed_stage=${stage};exit=${code}"
  mv "$status_part" "$STATUS"
  exit "$code"
}

record DEPENDENCY_REPAIR RUNNING r-r.utils
"$MAMBA" install -y --strict-channel-priority --override-channels \
  -c conda-forge -c bioconda -p "$ENV_PREFIX" r-r.utils \
  || finish_fail DEPENDENCY_REPAIR "$?" mamba_install_rutils_failed
record DEPENDENCY_REPAIR PASS r-r.utils

mkdir -p "$PROVENANCE"
record PREFLIGHT RUNNING package_load_and_gzip_fread
"$ENV_PREFIX/bin/Rscript" --vanilla -e \
  'suppressPackageStartupMessages({library(infercnv);library(rjags);library(data.table);library(jsonlite);library(R.utils)}); p <- "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase10_external_cnv_v1/full_inputs_v2/infercnv"; u <- sort(list.dirs(p,recursive=FALSE,full.names=TRUE))[1]; f <- file.path(u,"cell_id_mapping_v1.tsv.gz"); x <- fread(f,nrows=3,colClasses="character"); stopifnot(nrow(x)>0); cat("R=",R.version.string,"\n",sep=""); cat("infercnv=",as.character(packageVersion("infercnv")),"\n",sep=""); cat("R.utils=",as.character(packageVersion("R.utils")),"\n",sep=""); cat("gzip_fread=PASS\n")' \
  > "$PROVENANCE/preflight.txt.part" 2>&1 \
  || finish_fail PREFLIGHT "$?" package_or_gzip_fread_failed
mv "$PROVENANCE/preflight.txt.part" "$PROVENANCE/preflight.txt"
"$MAMBA" list -p "$ENV_PREFIX" --explicit > "$PROVENANCE/conda_explicit_spec.txt.part" \
  || finish_fail PROVENANCE "$?" explicit_spec_failed
mv "$PROVENANCE/conda_explicit_spec.txt.part" "$PROVENANCE/conda_explicit_spec.txt"
record PREFLIGHT PASS "$PROVENANCE/preflight.txt"

record RDS_VALIDATION RUNNING expected_infercnv=49_expected_copykat=62
"$ENV_PREFIX/bin/Rscript" --vanilla "$VALIDATOR" \
  --phase10-root "$ROOT/phase10_external_cnv_v1" \
  --out-dir "$RUN" \
  || finish_fail RDS_VALIDATION "$?" validator_failed

if [[ ! -s "$RUN/TERMINAL_RDS_VALIDATION_STATUS.json" ]] || \
   [[ ! -s "$RUN/terminal_rds_validation_v1.tsv" ]]; then
  finish_fail RDS_VALIDATION 4 required_outputs_missing
fi
record RDS_VALIDATION PASS "$RUN"
record OVERALL PASS repaired_R.utils_and_111_RDS_validation_complete
mv "$status_part" "$STATUS"

