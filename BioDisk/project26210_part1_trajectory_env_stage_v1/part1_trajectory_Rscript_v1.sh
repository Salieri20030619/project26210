#!/usr/bin/env bash
set -euo pipefail
export R_LIBS_USER=/mnt/BioDisk/ludingyu/project26210/envs/part2_6_2_7_20260726_v1/R_library
R453_LIB=/mnt/BioDisk/ludingyu/miniconda3/envs/r453/lib
TRAJECTORY_SYSTEM_LIB=/mnt/BioDisk/ludingyu/project26210/envs/part2_6_2_7_20260726_v1/system_libs
export LD_LIBRARY_PATH="$R453_LIB:$TRAJECTORY_SYSTEM_LIB/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export UDUNITS2_XML_PATH="$TRAJECTORY_SYSTEM_LIB/share/udunits/udunits2.xml"
exec /mnt/BioDisk/ludingyu/miniconda3/envs/r453/bin/Rscript --no-init-file "$@"
