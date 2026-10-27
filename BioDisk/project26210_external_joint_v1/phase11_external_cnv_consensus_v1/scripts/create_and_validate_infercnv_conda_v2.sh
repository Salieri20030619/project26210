#!/usr/bin/env bash
set -uo pipefail

ROOT=/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1
STAGE="$ROOT/phase11_external_cnv_consensus_v1"
MAMBA=/mnt/BioDisk/wangmingyuan/miniconda3/bin/mamba
ENV_PREFIX=/mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2
VALIDATOR="$STAGE/scripts/02_validate_terminal_external_cnv_rds_v1.R"
RUN="$STAGE/results/rds_validation_runs/20260812_final_v4"
STATUS="$STAGE/results/FINAL_BACKGROUND_CHAIN_STATUS_v2.tsv"
PROVENANCE="$STAGE/results/infercnv_conda_v2_provenance"

if [[ -e "$ENV_PREFIX" ]]; then
  echo "Refusing existing environment: $ENV_PREFIX" >&2
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
if [[ ! -x "$MAMBA" ]]; then
  echo "Mamba is not executable: $MAMBA" >&2
  exit 2
fi

status_part="${STATUS}.part"
printf 'stage\tresult\ttimestamp\tdetail\n' > "$status_part"
record() {
  printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$(date -Is)" "$3" >> "$status_part"
}
finish_fail() {
  local stage="$1"
  local code="$2"
  local detail="$3"
  record "$stage" FAIL "$detail"
  record OVERALL FAIL "failed_stage=${stage};exit=${code}"
  mv "$status_part" "$STATUS"
  exit "$code"
}

record ENV_CREATE RUNNING "prefix=$ENV_PREFIX"
"$MAMBA" create -y --strict-channel-priority --override-channels \
  -c conda-forge -c bioconda \
  -p "$ENV_PREFIX" \
  bioconductor-infercnv r-data.table r-jsonlite r-rjags jags \
  || finish_fail ENV_CREATE "$?" mamba_create_failed
record ENV_CREATE PASS "prefix=$ENV_PREFIX"

mkdir -p "$PROVENANCE"
"$MAMBA" list -p "$ENV_PREFIX" --explicit > "$PROVENANCE/conda_explicit_spec.txt.part" \
  || finish_fail PROVENANCE "$?" explicit_spec_failed
mv "$PROVENANCE/conda_explicit_spec.txt.part" "$PROVENANCE/conda_explicit_spec.txt"
"$MAMBA" env export -p "$ENV_PREFIX" > "$PROVENANCE/environment.yml.part" \
  || finish_fail PROVENANCE "$?" environment_export_failed
mv "$PROVENANCE/environment.yml.part" "$PROVENANCE/environment.yml"

record SELF_TEST RUNNING infercnv_rjags_data.table_jsonlite
"$ENV_PREFIX/bin/Rscript" --vanilla -e \
  'suppressPackageStartupMessages({library(infercnv);library(rjags);library(data.table);library(jsonlite)}); cat("R=",R.version.string,"\n",sep=""); cat("infercnv=",as.character(packageVersion("infercnv")),"\n",sep=""); cat("rjags=",as.character(packageVersion("rjags")),"\n",sep=""); cat("data.table=",as.character(packageVersion("data.table")),"\n",sep=""); cat("jsonlite=",as.character(packageVersion("jsonlite")),"\n",sep=""); cat("libPaths=",paste(.libPaths(),collapse="|"),"\n",sep=""); sessionInfo()' \
  > "$PROVENANCE/r_self_test.txt.part" 2>&1 \
  || finish_fail SELF_TEST "$?" package_load_or_abi_test_failed
mv "$PROVENANCE/r_self_test.txt.part" "$PROVENANCE/r_self_test.txt"
record SELF_TEST PASS "$PROVENANCE/r_self_test.txt"

record RDS_VALIDATION RUNNING "expected_infercnv=49;expected_copykat=62"
"$ENV_PREFIX/bin/Rscript" --vanilla "$VALIDATOR" \
  --phase10-root "$ROOT/phase10_external_cnv_v1" \
  --out-dir "$RUN" \
  || finish_fail RDS_VALIDATION "$?" validator_failed

if [[ ! -s "$RUN/TERMINAL_RDS_VALIDATION_STATUS.json" ]] || \
   [[ ! -s "$RUN/terminal_rds_validation_v1.tsv" ]]; then
  finish_fail RDS_VALIDATION 4 required_outputs_missing
fi
record RDS_VALIDATION PASS "$RUN"
record OVERALL PASS install_self_test_and_111_RDS_validation_complete
mv "$status_part" "$STATUS"

