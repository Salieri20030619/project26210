#!/usr/bin/env python3
"""Audit the strict Part 1 CNV-selected input before subclone and cNMF work."""

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
INPUT = ROOT / "03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad"
INPUT_METADATA = ROOT / "03_malignancy_cnv/cnv_annotated_h5ad_v1/run_metadata_v1.json"
OUT = ROOT / "04_subclones_tymp/input_audit_v1"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"
PYTHON = Path("/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python")
RSCRIPT = Path("/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript")
EXPECTED_SHA256 = "e4779c46af13a8643e8dec584eeff4322258181f5edab3e55dcd8c20bc28a9e7"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def gate(stage: str, records: list[dict[str, str]]) -> None:
    result = subprocess.run([str(MEMORY_GATE)], text=True, capture_output=True)
    output = " ".join((result.stdout + " " + result.stderr).split())
    records.append({"stage": stage, "exit_status": str(result.returncode), "gate_output": output})
    print(output, flush=True)
    if result.returncode != 0 or "PASS" not in output:
        raise RuntimeError(f"Memory gate blocked {stage}: {output}")


def r_capabilities() -> pd.DataFrame:
    expression = (
        'cat(R.version.string,"\\n"); '
        'for(p in c("Matrix","RcppML","data.table","cluster","edgeR","lme4")) '
        'cat(p,"\\t",requireNamespace(p,quietly=TRUE),"\\t",'
        'if(requireNamespace(p,quietly=TRUE)) as.character(packageVersion(p)) else "NA","\\n",sep="")'
    )
    result = subprocess.run([str(RSCRIPT), "-e", expression], text=True, capture_output=True, check=True)
    lines = [line for line in result.stdout.splitlines() if "\t" in line]
    rows = []
    for line in lines:
        package, available, version = line.split("\t")
        rows.append({"component": package, "available": available, "version": version})
    rows.append({"component": "R", "available": "TRUE", "version": result.stdout.splitlines()[0]})
    rows.append({"component": "Rscript_path", "available": "TRUE", "version": str(RSCRIPT)})
    rows.append({"component": "Python_path", "available": "TRUE", "version": str(PYTHON)})
    return pd.DataFrame(rows)


def main() -> None:
    if OUT.exists() and any(OUT.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty audit directory: {OUT}")
    OUT.mkdir(parents=True, exist_ok=True)
    gates: list[dict[str, str]] = []
    gate("open_input", gates)
    metadata = json.loads(INPUT_METADATA.read_text(encoding="utf-8"))
    if metadata["output_h5ad_sha256"] != EXPECTED_SHA256:
        raise RuntimeError("Input metadata SHA-256 does not match the frozen expected value")
    adata = ad.read_h5ad(INPUT, backed="r")
    if adata.shape != (258120, 37795):
        raise RuntimeError(f"Unexpected H5AD shape: {adata.shape}")
    if list(adata.layers.keys()) != ["counts"]:
        raise RuntimeError(f"Expected only raw counts layer, found {list(adata.layers.keys())}")

    obs = adata.obs
    eligible = obs["eligible_for_subclone_selection_v1"].astype(bool).to_numpy()
    primary = obs["included_in_primary_malignancy_summary_v1"].astype(bool).to_numpy()
    malignant = obs["malignant_conservative_v1"].astype(bool).to_numpy()
    expected = primary & malignant
    if not np.array_equal(eligible, expected) or int(eligible.sum()) != 11950:
        raise RuntimeError("Subclone eligibility is not exactly primary AND conservative malignant")
    selected = obs.loc[eligible].copy()
    if not selected["source_tissue"].astype(str).eq("Tumor").all():
        raise RuntimeError("Eligible registry includes non-tumor cells")
    if not selected["cell_type_broad_v1"].astype(str).eq("Epithelial cell").all():
        raise RuntimeError("Eligible registry includes non-clear-epithelial cells")
    mixed = obs["cell_type_broad_v1"].astype(str).eq("Epithelial/Myeloid mixed")
    if obs.loc[mixed, "eligible_for_subclone_selection_v1"].astype(bool).any():
        raise RuntimeError("Mixed-lineage cells entered the subclone selection")

    idx = np.flatnonzero(eligible)
    library_size = np.zeros(len(idx), dtype=np.float64)
    detected_genes = np.zeros(len(idx), dtype=np.int32)
    count_nnz = 0
    count_min = np.inf
    count_max = -np.inf
    noninteger = 0
    nonfinite = 0
    negative = 0
    block_size = 512
    for block_id, start in enumerate(range(0, len(idx), block_size), start=1):
        if block_id == 1 or block_id % 8 == 0:
            gate(f"scan_counts_block_{block_id}", gates)
        stop = min(start + block_size, len(idx))
        matrix = adata.layers["counts"][idx[start:stop], :]
        matrix = sp.csr_matrix(matrix)
        values = matrix.data
        count_nnz += int(matrix.nnz)
        if values.size:
            count_min = min(count_min, float(values.min()))
            count_max = max(count_max, float(values.max()))
            noninteger += int(np.count_nonzero(np.abs(values - np.rint(values)) > 1e-8))
            nonfinite += int(np.count_nonzero(~np.isfinite(values)))
            negative += int(np.count_nonzero(values < 0))
        library_size[start:stop] = np.asarray(matrix.sum(axis=1)).ravel()
        detected_genes[start:stop] = np.diff(matrix.indptr)
    if noninteger or nonfinite or negative or np.any(library_size <= 0):
        raise RuntimeError(
            f"Invalid selected raw counts: noninteger={noninteger}, nonfinite={nonfinite}, negative={negative}"
        )

    tymp_index = adata.var_names.get_indexer(["TYMP"])[0]
    if tymp_index < 0:
        raise RuntimeError("TYMP is absent from the input gene registry")
    tymp_counts = np.asarray(sp.csr_matrix(adata.layers["counts"][idx, tymp_index]).toarray()).ravel()
    adata.file.close()

    selected.insert(0, "exact_cell_id", selected.index.astype(str))
    selected["raw_library_size_v1"] = library_size
    selected["raw_detected_genes_v1"] = detected_genes
    selected["TYMP_raw_count_v1"] = tymp_counts
    selected["TYMP_detected_v1"] = tymp_counts > 0
    selected.to_csv(OUT / "subclone_cnmf_cell_registry_v1.tsv.gz", sep="\t", index=False, compression="gzip")

    def summarize(group: str) -> pd.DataFrame:
        table = selected.groupby(group, observed=True).agg(
            n_cells=("exact_cell_id", "size"),
            median_library_size=("raw_library_size_v1", "median"),
            median_detected_genes=("raw_detected_genes_v1", "median"),
            TYMP_positive_cells=("TYMP_detected_v1", "sum"),
            TYMP_positive_fraction=("TYMP_detected_v1", "mean"),
            TYMP_mean_raw_count=("TYMP_raw_count_v1", "mean"),
        ).reset_index()
        table["eligible_for_patient_level_test_min20"] = table["n_cells"] >= 20
        return table

    patient = summarize("curated_patient_id").sort_values(["n_cells", "curated_patient_id"])
    sample = summarize("sample_id").sort_values(["n_cells", "sample_id"])
    study = summarize("study").sort_values("study")
    patient.to_csv(OUT / "subclone_cnmf_patient_inventory_v1.tsv", sep="\t", index=False)
    sample.to_csv(OUT / "subclone_cnmf_sample_inventory_v1.tsv", sep="\t", index=False)
    study.to_csv(OUT / "subclone_cnmf_study_inventory_v1.tsv", sep="\t", index=False)

    capabilities = r_capabilities()
    capabilities.to_csv(OUT / "subclone_cnmf_environment_capabilities_v1.tsv", sep="\t", index=False)
    input_qc = pd.DataFrame(
        {
            "metric": [
                "n_cells", "n_genes", "n_patients", "n_samples", "n_studies",
                "n_patients_min20", "n_patients_below20", "n_count_nnz",
                "count_min_nonzero", "count_max_nonzero", "count_noninteger",
                "count_nonfinite", "count_negative", "library_size_min",
                "library_size_median", "library_size_max", "detected_genes_min",
                "detected_genes_median", "detected_genes_max", "TYMP_positive_cells",
                "TYMP_positive_fraction", "mixed_cells_selected", "selection_rule_exact",
            ],
            "value": [
                len(selected), 37795, selected["curated_patient_id"].nunique(),
                selected["sample_id"].nunique(), selected["study"].nunique(),
                int((patient["n_cells"] >= 20).sum()), int((patient["n_cells"] < 20).sum()),
                count_nnz, count_min, count_max, noninteger, nonfinite, negative,
                float(library_size.min()), float(np.median(library_size)), float(library_size.max()),
                int(detected_genes.min()), float(np.median(detected_genes)), int(detected_genes.max()),
                int((tymp_counts > 0).sum()), float(np.mean(tymp_counts > 0)), 0, True,
            ],
        }
    )
    input_qc.to_csv(OUT / "subclone_cnmf_input_qc_v1.tsv", sep="\t", index=False)
    pd.DataFrame(gates).to_csv(OUT / "memory_gate_log_v1.tsv", sep="\t", index=False)
    run_metadata = {
        "schema_version": "part1_subclone_cnmf_input_audit_v1",
        "completed_at": datetime.now().astimezone().isoformat(),
        "input_h5ad": str(INPUT),
        "input_h5ad_size_bytes": INPUT.stat().st_size,
        "input_h5ad_sha256": EXPECTED_SHA256,
        "selection_rule": "eligible_for_subclone_selection_v1 == True",
        "n_cells": len(selected),
        "n_patients": int(selected["curated_patient_id"].nunique()),
        "n_patients_min20": int((patient["n_cells"] >= 20).sum()),
        "n_patients_below20": int((patient["n_cells"] < 20).sum()),
        "cnmf_engine_available": "RcppML 0.3.7.1 in R_4.5.2",
        "python_cnmf_installed": False,
        "dependency_changes": False,
        "interpretation": "all 11950 cells may enter program discovery; patient-level summaries require an explicit minimum-cell rule",
    }
    (OUT / "run_metadata_v1.json").write_text(json.dumps(run_metadata, indent=2) + "\n", encoding="utf-8")
    manifest = []
    for path in sorted(OUT.iterdir()):
        if path.is_file() and path.name not in {"file_manifest_v1.tsv", "RUN_COMPLETE"}:
            manifest.append({"file": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest).to_csv(OUT / "file_manifest_v1.tsv", sep="\t", index=False)
    (OUT / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n", encoding="utf-8")
    print(json.dumps(run_metadata, indent=2))
    print("PART1_SUBCLONE_CNMF_INPUT_AUDIT_V1_PASS")


if __name__ == "__main__":
    main()
