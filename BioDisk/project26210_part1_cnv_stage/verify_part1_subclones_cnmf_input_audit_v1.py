#!/usr/bin/env python3
"""Independent verification of the strict subclone/cNMF input audit."""

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
INPUT = ROOT / "03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad"
AUDIT = ROOT / "04_subclones_tymp/input_audit_v1"
OUT = ROOT / "04_subclones_tymp/input_audit_verification_v1"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if OUT.exists() and any(OUT.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty verification directory: {OUT}")
    OUT.mkdir(parents=True, exist_ok=True)
    checks: list[dict[str, str]] = []

    def check(name: str, passed: bool, evidence: object) -> None:
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "evidence": str(evidence)})
        if not passed:
            raise RuntimeError(f"Verification failed: {name}: {evidence}")

    gate = subprocess.run([str(MEMORY_GATE)], text=True, capture_output=True)
    gate_text = " ".join((gate.stdout + " " + gate.stderr).split())
    print(gate_text, flush=True)
    check("memory_gate", gate.returncode == 0 and "PASS" in gate_text, gate_text)
    check("audit_complete", (AUDIT / "RUN_COMPLETE").is_file(), AUDIT / "RUN_COMPLETE")

    adata = ad.read_h5ad(INPUT, backed="r")
    eligible = adata.obs["eligible_for_subclone_selection_v1"].astype(bool)
    expected = adata.obs["included_in_primary_malignancy_summary_v1"].astype(bool) & adata.obs["malignant_conservative_v1"].astype(bool)
    check("input_shape", adata.shape == (258120, 37795), adata.shape)
    check("selection_rule_exact", eligible.equals(expected), int((eligible != expected).sum()))
    check("eligible_count", int(eligible.sum()) == 11950, int(eligible.sum()))
    selected_ids = adata.obs_names[eligible].astype(str)
    selected = adata.obs.loc[eligible]
    check("all_tumor", selected["source_tissue"].astype(str).eq("Tumor").all(), len(selected))
    check("all_clear_epithelial", selected["cell_type_broad_v1"].astype(str).eq("Epithelial cell").all(), len(selected))
    check("mixed_excluded", not adata.obs.loc[adata.obs["cell_type_broad_v1"].astype(str).eq("Epithelial/Myeloid mixed"), "eligible_for_subclone_selection_v1"].astype(bool).any(), "0 selected")
    adata.file.close()

    registry = pd.read_csv(AUDIT / "subclone_cnmf_cell_registry_v1.tsv.gz", sep="\t", dtype={"exact_cell_id": str})
    check("registry_rows", len(registry) == 11950, len(registry))
    check("registry_unique", not registry["exact_cell_id"].duplicated().any(), registry["exact_cell_id"].nunique())
    check("registry_order", np.array_equal(registry["exact_cell_id"].to_numpy(str), selected_ids.to_numpy(str)), "exact")
    check("registry_integer_counts", registry["raw_library_size_v1"].ge(0).all() and np.allclose(registry["raw_library_size_v1"], np.rint(registry["raw_library_size_v1"])), "all")
    check("registry_detected_genes", registry["raw_detected_genes_v1"].between(1, 37795).all(), [registry["raw_detected_genes_v1"].min(), registry["raw_detected_genes_v1"].max()])
    check("registry_tymp_nonnegative_integer", registry["TYMP_raw_count_v1"].ge(0).all() and np.allclose(registry["TYMP_raw_count_v1"], np.rint(registry["TYMP_raw_count_v1"])), "all")
    check("registry_tymp_detected", np.array_equal(registry["TYMP_detected_v1"].astype(bool), registry["TYMP_raw_count_v1"].gt(0)), int(registry["TYMP_detected_v1"].sum()))

    patient = pd.read_csv(AUDIT / "subclone_cnmf_patient_inventory_v1.tsv", sep="\t")
    sample = pd.read_csv(AUDIT / "subclone_cnmf_sample_inventory_v1.tsv", sep="\t")
    study = pd.read_csv(AUDIT / "subclone_cnmf_study_inventory_v1.tsv", sep="\t")
    check("patient_rows", len(patient) == 26 and int(patient["n_cells"].sum()) == 11950, [len(patient), patient["n_cells"].sum()])
    check("sample_rows", len(sample) == 29 and int(sample["n_cells"].sum()) == 11950, [len(sample), sample["n_cells"].sum()])
    check("study_rows", len(study) == 5 and int(study["n_cells"].sum()) == 11950, [len(study), study["n_cells"].sum()])
    check("patient_min20", int(patient["eligible_for_patient_level_test_min20"].astype(bool).sum()) == 20, int(patient["eligible_for_patient_level_test_min20"].astype(bool).sum()))
    check("patient_below20", int((~patient["eligible_for_patient_level_test_min20"].astype(bool)).sum()) == 6, int((~patient["eligible_for_patient_level_test_min20"].astype(bool)).sum()))

    qc = pd.read_csv(AUDIT / "subclone_cnmf_input_qc_v1.tsv", sep="\t").set_index("metric")["value"]
    check("qc_count_nnz", int(float(qc["n_count_nnz"])) == int(registry["raw_detected_genes_v1"].sum()), [qc["n_count_nnz"], registry["raw_detected_genes_v1"].sum()])
    check("qc_tymp_positive", int(float(qc["TYMP_positive_cells"])) == int(registry["TYMP_detected_v1"].sum()) == 8178, qc["TYMP_positive_cells"])
    check("qc_invalid_counts_zero", all(int(float(qc[key])) == 0 for key in ["count_noninteger", "count_nonfinite", "count_negative"]), "0/0/0")
    check("qc_mixed_zero", int(float(qc["mixed_cells_selected"])) == 0, qc["mixed_cells_selected"])

    capability = pd.read_csv(AUDIT / "subclone_cnmf_environment_capabilities_v1.tsv", sep="\t", dtype=str)
    rcppml = capability.loc[capability["component"].eq("RcppML")].iloc[0]
    check("rcppml_available", rcppml["available"] == "TRUE" and rcppml["version"] == "0.3.7.1", rcppml.to_dict())
    metadata = json.loads((AUDIT / "run_metadata_v1.json").read_text(encoding="utf-8"))
    check("input_hash_frozen", metadata["input_h5ad_sha256"] == "e4779c46af13a8643e8dec584eeff4322258181f5edab3e55dcd8c20bc28a9e7", metadata["input_h5ad_sha256"])
    check("no_dependency_changes", metadata["dependency_changes"] is False, metadata["dependency_changes"])

    manifest = pd.read_csv(AUDIT / "file_manifest_v1.tsv", sep="\t", dtype=str)
    for row in manifest.itertuples(index=False):
        path = AUDIT / row.file
        check(f"hash_{row.file}", path.is_file() and sha256(path) == row.sha256, row.sha256)

    pd.DataFrame(checks).to_csv(OUT / "verification_checks_v1.tsv", sep="\t", index=False)
    run_metadata = {
        "schema_version": "part1_subclone_cnmf_input_audit_verification_v1",
        "completed_at": datetime.now().astimezone().isoformat(),
        "n_checks": len(checks),
        "n_pass": sum(row["status"] == "PASS" for row in checks),
        "n_fail": sum(row["status"] == "FAIL" for row in checks),
    }
    (OUT / "run_metadata_v1.json").write_text(json.dumps(run_metadata, indent=2) + "\n", encoding="utf-8")
    (OUT / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n", encoding="utf-8")
    print(json.dumps(run_metadata, indent=2))
    print("PART1_SUBCLONE_CNMF_INPUT_AUDIT_VERIFICATION_V1_PASS")


if __name__ == "__main__":
    main()
