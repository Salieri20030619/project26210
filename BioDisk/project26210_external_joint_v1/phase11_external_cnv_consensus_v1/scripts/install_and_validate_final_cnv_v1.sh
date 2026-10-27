#!/usr/bin/env bash
set -uo pipefail

ROOT=/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1
STAGE="$ROOT/phase11_external_cnv_consensus_v1"
ENV_ROOT=/mnt/BioDisk/wangmingyuan/project26210_envs
LIB="$ENV_ROOT/phase11_infercnv_R_lib_v1"
INSTALLER="$STAGE/scripts/install_infercnv_project_lib_v1.R"
VALIDATOR="$STAGE/scripts/02_validate_terminal_external_cnv_rds_v1.R"
RUN="$STAGE/results/rds_validation_runs/20260812_120500_final_v3"
STATUS="$STAGE/results/FINAL_BACKGROUND_CHAIN_STATUS_v1.tsv"

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

if Rscript "$INSTALLER" "$LIB"; then
  printf 'infercnv_install\tPASS\t%s\t%s\n' "$(date -Is)" "$LIB" >> "$status_part"
else
  code=$?
  printf 'infercnv_install\tFAIL\t%s\texit_%s\n' "$(date -Is)" "$code" >> "$status_part"
  mv "$status_part" "$STATUS"
  exit "$code"
fi

export R_LIBS_USER="$LIB"
if Rscript "$VALIDATOR" --phase10-root "$ROOT/phase10_external_cnv_v1" --out-dir "$RUN"; then
  printf 'terminal_rds_validation\tPASS\t%s\t%s\n' "$(date -Is)" "$RUN" >> "$status_part"
  printf 'OVERALL\tPASS\t%s\tinstall_and_111_RDS_validation_complete\n' "$(date -Is)" >> "$status_part"
  mv "$status_part" "$STATUS"
else
  code=$?
  printf 'terminal_rds_validation\tFAIL\t%s\texit_%s\n' "$(date -Is)" "$code" >> "$status_part"
  printf 'OVERALL\tFAIL\t%s\tsee_background_log\n' "$(date -Is)" >> "$status_part"
  mv "$status_part" "$STATUS"
  exit "$code"
fi
