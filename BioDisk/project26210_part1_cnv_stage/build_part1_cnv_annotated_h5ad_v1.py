#!/usr/bin/env python3
"""Create a versioned H5AD carrying primary and isolated mixed CNV annotations."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1"
SOURCE = ROOT / "02_integration_annotation/cell_annotation_v1/part1_scvi_integrated_cell_annotated_v1.h5ad"
CNV_ROOT = ROOT / "03_malignancy_cnv"
PRIMARY = CNV_ROOT / "malignancy_cnv_summary_v1/tumor_epithelial_malignancy_calls_v1.tsv.gz"
MIXED = CNV_ROOT / "mixed_cnv_sensitivity_v1/summary_v2_with_control/tumor_mixed_cnv_calls_with_control_v2.tsv"
OUT_DIR = CNV_ROOT / "cnv_annotated_h5ad_v1"
OUTPUT = OUT_DIR / "part1_scvi_integrated_cell_annotated_cnv_v1.h5ad"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def memory_gate(stage: str, records: list[dict[str, str]]) -> None:
    result = subprocess.run([str(MEMORY_GATE)], text=True, capture_output=True)
    output = " ".join((result.stdout + " " + result.stderr).split())
    records.append({"stage": stage, "exit_status": str(result.returncode), "gate_output": output})
    print(output, flush=True)
    if result.returncode != 0 or "PASS" not in output:
        raise RuntimeError(f"Memory gate blocked {stage}: {output}")


def main() -> None:
    if not (CNV_ROOT / "malignancy_cnv_summary_v1/RUN_COMPLETE").is_file():
        raise RuntimeError("Primary malignancy summary is incomplete")
    if not (CNV_ROOT / "mixed_cnv_sensitivity_v1/summary_v1/RUN_COMPLETE").is_file():
        raise RuntimeError("Mixed CNV sensitivity summary is incomplete")
    if not (CNV_ROOT / "mixed_cnv_sensitivity_v1/summary_v2_with_control/RUN_COMPLETE").is_file():
        raise RuntimeError("Mixed CNV matched-control summary is incomplete")
    if not (CNV_ROOT / "mixed_cnv_sensitivity_verification_v1/RUN_COMPLETE").is_file():
        raise RuntimeError("Mixed CNV sensitivity verification is incomplete")
    if OUT_DIR.exists() and any(OUT_DIR.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty H5AD output directory: {OUT_DIR}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    gates: list[dict[str, str]] = []
    memory_gate("read_primary_and_mixed_calls", gates)
    primary = pd.read_csv(PRIMARY, sep="\t", dtype={"exact_cell_id": str})
    mixed = pd.read_csv(MIXED, sep="\t", dtype={"exact_cell_id": str})
    if len(primary) != 38670 or primary["exact_cell_id"].duplicated().any():
        raise RuntimeError(f"Unexpected primary call registry: {len(primary)}")
    if len(mixed) != 2 or mixed["exact_cell_id"].duplicated().any():
        raise RuntimeError(f"Unexpected mixed call registry: {len(mixed)}")
    if set(primary["exact_cell_id"]).intersection(mixed["exact_cell_id"]):
        raise RuntimeError("Primary and mixed call registries overlap")

    memory_gate("load_source_h5ad", gates)
    adata = ad.read_h5ad(SOURCE)
    if adata.shape != (258120, 37795) or adata.obs_names.has_duplicates or adata.var_names.has_duplicates:
        raise RuntimeError(f"Unexpected source H5AD registry: {adata.shape}")
    original_obs_columns = list(adata.obs.columns)
    primary_index = adata.obs_names.get_indexer(primary["exact_cell_id"])
    mixed_index = adata.obs_names.get_indexer(mixed["exact_cell_id"])
    if np.any(primary_index < 0) or np.any(mixed_index < 0):
        raise RuntimeError("CNV calls contain cell IDs absent from source H5AD")
    primary_obs = adata.obs.iloc[primary_index]
    mixed_obs = adata.obs.iloc[mixed_index]
    if not (
        primary_obs["source_tissue"].astype(str).eq("Tumor").all()
        and primary_obs["cell_type_broad_v1"].astype(str).eq("Epithelial cell").all()
    ):
        raise RuntimeError("Primary calls do not map exclusively to clear tumor epithelial cells")
    if not (
        mixed_obs["source_tissue"].astype(str).eq("Tumor").all()
        and mixed_obs["cell_type_broad_v1"].astype(str).eq("Epithelial/Myeloid mixed").all()
    ):
        raise RuntimeError("Mixed calls do not map exclusively to tumor mixed-lineage cells")

    obs = adata.obs
    obs["cnv_assessed_v1"] = False
    obs["cnv_analysis_scope_v1"] = "not_assessed"
    obs["infercnv_call_v1"] = "not_assessed"
    obs["copykat_call_v1"] = "not_assessed"
    obs["malignancy_status_v1"] = "not_assessed"
    obs["malignant_conservative_v1"] = False
    obs["malignancy_evidence_v1"] = "not_assessed"
    obs["infercnv_unit_v1"] = "not_assessed"
    obs["copykat_unit_v1"] = "not_assessed"
    obs["included_in_primary_malignancy_summary_v1"] = False
    obs["eligible_for_subclone_selection_v1"] = False
    obs["mixed_cnv_raw_status_v1"] = "not_applicable"
    obs["mixed_cnv_infercnv_reliability_v2"] = "not_applicable"
    obs["mixed_cnv_copykat_reliability_v2"] = "not_applicable"
    obs["mixed_cnv_interpretation_v2"] = "not_applicable"

    pids = primary["exact_cell_id"].to_numpy(str)
    obs.loc[pids, "cnv_assessed_v1"] = True
    obs.loc[pids, "cnv_analysis_scope_v1"] = "primary_clear_tumor_epithelial"
    obs.loc[pids, "infercnv_call_v1"] = primary["infercnv_call"].fillna("not_assessed").to_numpy(str)
    obs.loc[pids, "copykat_call_v1"] = primary["copykat_call"].fillna("not_assessed").to_numpy(str)
    obs.loc[pids, "malignancy_status_v1"] = primary["malignancy_status_v1"].fillna("not_assessed").to_numpy(str)
    obs.loc[pids, "malignant_conservative_v1"] = primary["malignant_conservative_v1"].astype(bool).to_numpy()
    obs.loc[pids, "malignancy_evidence_v1"] = primary["malignancy_evidence_v1"].fillna("not_assessed").to_numpy(str)
    obs.loc[pids, "infercnv_unit_v1"] = primary["infercnv_unit"].fillna("not_eligible").to_numpy(str)
    obs.loc[pids, "copykat_unit_v1"] = primary["copykat_unit"].fillna("not_eligible").to_numpy(str)
    obs.loc[pids, "included_in_primary_malignancy_summary_v1"] = True
    obs.loc[pids, "eligible_for_subclone_selection_v1"] = primary["malignant_conservative_v1"].astype(bool).to_numpy()

    mids = mixed["exact_cell_id"].to_numpy(str)
    obs.loc[mids, "cnv_assessed_v1"] = True
    obs.loc[mids, "cnv_analysis_scope_v1"] = "mixed_tumor_sensitivity_only"
    obs.loc[mids, "infercnv_call_v1"] = mixed["infercnv_call"].to_numpy(str)
    obs.loc[mids, "copykat_call_v1"] = mixed["copykat_call"].to_numpy(str)
    obs.loc[mids, "malignancy_status_v1"] = "mixed_sensitivity_excluded_unstable_infercnv"
    obs.loc[mids, "malignant_conservative_v1"] = False
    obs.loc[mids, "malignancy_evidence_v1"] = "copykat_stable_support_but_infercnv_plus_one_background_unstable"
    obs.loc[mids, "infercnv_unit_v1"] = mixed["infercnv_unit"].to_numpy(str)
    obs.loc[mids, "copykat_unit_v1"] = mixed["copykat_unit"].to_numpy(str)
    obs.loc[mids, "mixed_cnv_raw_status_v1"] = mixed["malignancy_status_mixed_sensitivity_v1"].to_numpy(str)
    obs.loc[mids, "mixed_cnv_infercnv_reliability_v2"] = mixed["infercnv_sensitivity_reliability_v2"].to_numpy(str)
    obs.loc[mids, "mixed_cnv_copykat_reliability_v2"] = mixed["copykat_sensitivity_reliability_v2"].to_numpy(str)
    obs.loc[mids, "mixed_cnv_interpretation_v2"] = mixed["combined_sensitivity_interpretation_v2"].to_numpy(str)

    new_columns = list(obs.columns[len(original_obs_columns):])
    if int(obs["cnv_assessed_v1"].sum()) != 38672:
        raise RuntimeError("Unexpected number of assessed cells")
    if int(obs["included_in_primary_malignancy_summary_v1"].sum()) != 38670:
        raise RuntimeError("Unexpected number of primary-summary cells")
    if int(obs["malignant_conservative_v1"].sum()) != 11950:
        raise RuntimeError("Conservative malignant count does not reconcile")
    if int(obs["eligible_for_subclone_selection_v1"].sum()) != 11950:
        raise RuntimeError("Unexpected number of cells eligible for subclone selection")

    adata.uns["part1_malignancy_cnv_v1"] = {
        "schema_version": "part1_malignancy_cnv_annotation_v1",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_h5ad": str(SOURCE.resolve()),
        "primary_call_table": str(PRIMARY.resolve()),
        "mixed_sensitivity_control_call_table": str(MIXED.resolve()),
        "n_primary_clear_tumor_epithelial": 38670,
        "n_mixed_tumor_sensitivity": 2,
        "n_mixed_normal_not_tumor_queries": 916,
        "n_primary_malignant_conservative": 11950,
        "mixed_included_in_primary_summary": False,
        "mixed_eligible_for_subclone_selection": False,
        "selection_rule": "eligible_for_subclone_selection_v1 == True; this is equivalent to included_in_primary_malignancy_summary_v1 and malignant_conservative_v1",
        "mixed_interpretation": "CopyKAT support was stable, but inferCNV clear-cell calls were unstable after adding one mixed query; mixed cells remain excluded",
        "new_obs_columns": np.asarray(new_columns, dtype=str),
    }
    memory_gate("write_cnv_annotated_h5ad", gates)
    temporary = OUTPUT.with_suffix(".h5ad.tmp")
    adata.write_h5ad(temporary, compression="gzip", compression_opts=4)
    temporary.replace(OUTPUT)
    pd.DataFrame(gates).to_csv(OUT_DIR / "memory_gate_log_v1.tsv", sep="\t", index=False, lineterminator="\n")
    metadata = {
        "schema_version": "part1_cnv_annotated_h5ad_build_v1",
        "completed_at": datetime.now().astimezone().isoformat(),
        "source_h5ad": str(SOURCE.resolve()),
        "source_h5ad_size_bytes": SOURCE.stat().st_size,
        "source_h5ad_sha256": sha256(SOURCE),
        "output_h5ad": str(OUTPUT.resolve()),
        "output_h5ad_size_bytes": OUTPUT.stat().st_size,
        "output_h5ad_sha256": sha256(OUTPUT),
        "shape": list(adata.shape),
        "n_original_obs_columns": len(original_obs_columns),
        "n_new_obs_columns": len(new_columns),
        "new_obs_columns": new_columns,
        "n_cnv_assessed": 38672,
        "n_primary_summary_cells": 38670,
        "n_mixed_sensitivity_cells": 2,
        "n_normal_mixed_not_assessed": 916,
        "n_malignant_conservative": 11950,
        "n_eligible_for_subclone_selection": 11950,
    }
    (OUT_DIR / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    (OUT_DIR / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))
    print("PART1_CNV_ANNOTATED_H5AD_V1_BUILD_PASS")


if __name__ == "__main__":
    main()
