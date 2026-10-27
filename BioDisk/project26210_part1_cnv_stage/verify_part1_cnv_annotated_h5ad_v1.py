#!/usr/bin/env python3
"""Independently verify the versioned CNV-annotated H5AD and source-data preservation."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1"
SOURCE = ROOT / "02_integration_annotation/cell_annotation_v1/part1_scvi_integrated_cell_annotated_v1.h5ad"
CNV_ROOT = ROOT / "03_malignancy_cnv"
PRIMARY = CNV_ROOT / "malignancy_cnv_summary_v1/tumor_epithelial_malignancy_calls_v1.tsv.gz"
MIXED = CNV_ROOT / "mixed_cnv_sensitivity_v1/summary_v2_with_control/tumor_mixed_cnv_calls_with_control_v2.tsv"
H5AD_DIR = CNV_ROOT / "cnv_annotated_h5ad_v1"
OUTPUT = H5AD_DIR / "part1_scvi_integrated_cell_annotated_cnv_v1.h5ad"
VERIFY = CNV_ROOT / "cnv_annotated_h5ad_verification_v1"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if VERIFY.exists() and any(VERIFY.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty verification directory: {VERIFY}")
    VERIFY.mkdir(parents=True, exist_ok=True)
    checks: list[dict[str, str]] = []
    gates: list[dict[str, str]] = []

    def check(name: str, passed: bool, evidence: str) -> None:
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "evidence": evidence})
        if not passed:
            raise RuntimeError(f"Verification failed: {name}: {evidence}")

    def gate(stage: str) -> None:
        result = subprocess.run([str(MEMORY_GATE)], text=True, capture_output=True)
        output = " ".join((result.stdout + " " + result.stderr).split())
        gates.append({"stage": stage, "exit_status": str(result.returncode), "gate_output": output})
        print(output, flush=True)
        check(f"memory_gate_{stage}", result.returncode == 0 and "PASS" in output, output)

    def arrays_equal(left: object, right: object) -> bool:
        if sp.issparse(left) or sp.issparse(right):
            a = sp.csr_matrix(left)
            b = sp.csr_matrix(right)
            if a.shape != b.shape:
                return False
            diff = a - b
            diff.eliminate_zeros()
            return diff.nnz == 0
        a = np.asarray(left)
        b = np.asarray(right)
        return a.shape == b.shape and np.array_equal(a, b, equal_nan=True)

    def compare_axis0(name: str, left: object, right: object, n_rows: int, block: int = 2048) -> None:
        gate(f"compare_{name}")
        for start in range(0, n_rows, block):
            stop = min(start + block, n_rows)
            if not arrays_equal(left[start:stop], right[start:stop]):
                check(f"preserve_{name}", False, f"difference in rows {start}:{stop}")
        check(f"preserve_{name}", True, f"all {n_rows} rows exact")

    gate("open_h5ad")
    check("source_exists", SOURCE.is_file() and SOURCE.stat().st_size > 0, str(SOURCE))
    check("annotated_exists", OUTPUT.is_file() and OUTPUT.stat().st_size > 0, str(OUTPUT))
    source = ad.read_h5ad(SOURCE, backed="r")
    annotated = ad.read_h5ad(OUTPUT, backed="r")
    check("shape", source.shape == annotated.shape == (258120, 37795), str(annotated.shape))
    check("obs_names", source.obs_names.equals(annotated.obs_names), f"n={annotated.n_obs}")
    check("var_names", source.var_names.equals(annotated.var_names), f"n={annotated.n_vars}")
    pd.testing.assert_frame_equal(source.obs, annotated.obs[source.obs.columns], check_dtype=True)
    check("original_obs_columns_exact", True, f"n={len(source.obs.columns)}")
    pd.testing.assert_frame_equal(source.var, annotated.var, check_dtype=True)
    check("var_table_exact", True, f"shape={source.var.shape}")
    check("layers_keys", set(source.layers.keys()) == set(annotated.layers.keys()), str(list(annotated.layers.keys())))
    check("obsm_keys", set(source.obsm.keys()) == set(annotated.obsm.keys()), str(list(annotated.obsm.keys())))
    check("varm_keys", set(source.varm.keys()) == set(annotated.varm.keys()), str(list(annotated.varm.keys())))
    check("obsp_keys", set(source.obsp.keys()) == set(annotated.obsp.keys()), str(list(annotated.obsp.keys())))
    check("varp_keys", set(source.varp.keys()) == set(annotated.varp.keys()), str(list(annotated.varp.keys())))
    check("uns_keys", set(source.uns.keys()).issubset(annotated.uns.keys()), str(list(annotated.uns.keys())))
    check("cnv_uns_present", "part1_malignancy_cnv_v1" in annotated.uns, "new CNV provenance block")

    compare_axis0("X", source.X, annotated.X, source.n_obs)
    for key in source.layers.keys():
        compare_axis0(f"layer_{key}", source.layers[key], annotated.layers[key], source.n_obs)
    for key in source.obsm.keys():
        compare_axis0(f"obsm_{key}", source.obsm[key], annotated.obsm[key], source.n_obs, block=8192)
    for key in source.obsp.keys():
        compare_axis0(f"obsp_{key}", source.obsp[key], annotated.obsp[key], source.n_obs, block=2048)
    for key in source.varm.keys():
        compare_axis0(f"varm_{key}", source.varm[key], annotated.varm[key], source.n_vars, block=8192)
    for key in source.varp.keys():
        compare_axis0(f"varp_{key}", source.varp[key], annotated.varp[key], source.n_vars, block=2048)
    check("raw_presence", (source.raw is None) == (annotated.raw is None), f"raw={source.raw is not None}")
    if source.raw is not None:
        pd.testing.assert_frame_equal(source.raw.var, annotated.raw.var, check_dtype=True)
        check("raw_var_exact", True, f"shape={source.raw.var.shape}")
        compare_axis0("raw_X", source.raw.X, annotated.raw.X, source.n_obs)

    obs = annotated.obs
    required = [
        "cnv_assessed_v1", "cnv_analysis_scope_v1", "infercnv_call_v1", "copykat_call_v1",
        "malignancy_status_v1", "malignant_conservative_v1", "malignancy_evidence_v1",
        "infercnv_unit_v1", "copykat_unit_v1", "included_in_primary_malignancy_summary_v1",
        "eligible_for_subclone_selection_v1", "mixed_cnv_raw_status_v1",
        "mixed_cnv_infercnv_reliability_v2", "mixed_cnv_copykat_reliability_v2",
        "mixed_cnv_interpretation_v2",
    ]
    check("new_obs_columns", all(column in obs.columns for column in required), ",".join(required))
    check("n_cnv_assessed", int(obs["cnv_assessed_v1"].sum()) == 38672, str(int(obs["cnv_assessed_v1"].sum())))
    check(
        "n_primary_summary",
        int(obs["included_in_primary_malignancy_summary_v1"].sum()) == 38670,
        str(int(obs["included_in_primary_malignancy_summary_v1"].sum())),
    )
    primary = pd.read_csv(PRIMARY, sep="\t")
    mixed = pd.read_csv(MIXED, sep="\t")
    p = obs.loc[primary["exact_cell_id"]]
    m = obs.loc[mixed["exact_cell_id"]]
    check("primary_scope", p["cnv_analysis_scope_v1"].astype(str).eq("primary_clear_tumor_epithelial").all(), f"n={len(p)}")
    check("mixed_scope", m["cnv_analysis_scope_v1"].astype(str).eq("mixed_tumor_sensitivity_only").all(), f"n={len(m)}")
    check("primary_infer_calls", np.array_equal(p["infercnv_call_v1"].astype(str), primary["infercnv_call"].astype(str)), "exact")
    check("primary_copykat_calls", np.array_equal(p["copykat_call_v1"].astype(str), primary["copykat_call"].astype(str)), "exact")
    check("primary_status", np.array_equal(p["malignancy_status_v1"].astype(str), primary["malignancy_status_v1"].astype(str)), "exact")
    check("primary_malignant_bool", np.array_equal(p["malignant_conservative_v1"].astype(bool), primary["malignant_conservative_v1"].astype(bool)), "exact")
    check("mixed_infer_calls", np.array_equal(m["infercnv_call_v1"].astype(str), mixed["infercnv_call"].astype(str)), "exact")
    check("mixed_copykat_calls", np.array_equal(m["copykat_call_v1"].astype(str), mixed["copykat_call"].astype(str)), "exact")
    check("mixed_not_primary", not m["included_in_primary_malignancy_summary_v1"].any(), "2/2 false")
    check("mixed_not_conservative_malignant", not m["malignant_conservative_v1"].any(), "2/2 false")
    check("mixed_not_subclone_eligible", not m["eligible_for_subclone_selection_v1"].any(), "2/2 false")
    check(
        "mixed_infercnv_unstable_recorded",
        m["mixed_cnv_infercnv_reliability_v2"].astype(str).eq("unstable_clear_background").all(),
        "2/2 unstable_clear_background",
    )
    check(
        "mixed_copykat_stable_recorded",
        m["mixed_cnv_copykat_reliability_v2"].astype(str).eq("stable_clear_background").all(),
        "2/2 stable_clear_background",
    )
    normal_mixed = obs["cell_type_broad_v1"].astype(str).eq("Epithelial/Myeloid mixed") & obs["source_tissue"].astype(str).eq("Normal")
    check("normal_mixed_count", int(normal_mixed.sum()) == 916, str(int(normal_mixed.sum())))
    check("normal_mixed_not_assessed", not obs.loc[normal_mixed, "cnv_assessed_v1"].any(), "916/916 false")
    check("primary_conservative_count", int(p["malignant_conservative_v1"].sum()) == 11950, str(int(p["malignant_conservative_v1"].sum())))
    check("global_conservative_count", int(obs["malignant_conservative_v1"].sum()) == 11950, str(int(obs["malignant_conservative_v1"].sum())))
    check("subclone_eligible_count", int(obs["eligible_for_subclone_selection_v1"].sum()) == 11950, str(int(obs["eligible_for_subclone_selection_v1"].sum())))
    source.file.close()
    annotated.file.close()

    pd.DataFrame(checks).to_csv(VERIFY / "verification_checks_v1.tsv", sep="\t", index=False, lineterminator="\n")
    pd.DataFrame(gates).to_csv(VERIFY / "memory_gate_log_v1.tsv", sep="\t", index=False, lineterminator="\n")
    manifest = []
    for path in [SOURCE, PRIMARY, MIXED, OUTPUT, H5AD_DIR / "run_metadata_v1.json"]:
        manifest.append({"path": str(path.resolve()), "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest).to_csv(VERIFY / "critical_file_manifest_v1.tsv", sep="\t", index=False, lineterminator="\n")
    metadata = {
        "schema_version": "part1_cnv_annotated_h5ad_verification_v1",
        "completed_at": datetime.now().astimezone().isoformat(),
        "n_checks": len(checks),
        "n_pass": sum(row["status"] == "PASS" for row in checks),
        "n_fail": sum(row["status"] == "FAIL" for row in checks),
    }
    (VERIFY / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    (VERIFY / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))
    print("PART1_CNV_ANNOTATED_H5AD_V1_VERIFICATION_PASS")


if __name__ == "__main__":
    main()
