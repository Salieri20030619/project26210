#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import pickle
import subprocess
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1/07_pyscenic_v1"
OUT = ROOT / "terminal_verification_v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_one_column(path: Path) -> list[str]:
    return pd.read_csv(path, sep="\t", dtype=str).iloc[:, 0].astype(str).tolist()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    gate = subprocess.run(
        [str(PROJECT / "check_memory_gate.sh"), "part1_pyscenic_terminal_verification_v1"],
        check=True,
        text=True,
        capture_output=True,
    )
    (OUT / "memory_gate_v1.txt").write_text(gate.stdout + gate.stderr, encoding="utf-8")
    checks: list[dict[str, str]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})

    cells = read_one_column(ROOT / "input/cells_v1.tsv")
    metadata = pd.read_csv(ROOT / "input/cell_metadata_v1.tsv.gz", sep="\t", dtype=str)
    check("frozen_cell_count", len(cells) == 11950 and len(set(cells)) == 11950, f"n={len(cells)};unique={len(set(cells))}")
    check("metadata_cell_order", metadata["exact_cell_id"].astype(str).tolist() == cells, f"rows={len(metadata)}")
    check("cohort_scope", metadata["curated_patient_id"].nunique() == 26 and metadata["study"].nunique() == 5, f"patients={metadata['curated_patient_id'].nunique()};studies={metadata['study'].nunique()}")
    check("completion_markers", (ROOT / "PYSCENIC_CORE_COMPLETE").exists() and (ROOT / "OVERNIGHT_COMPLETE").exists(), "PYSCENIC_CORE_COMPLETE;OVERNIGHT_COMPLETE")

    auc_by_method: dict[str, pd.DataFrame] = {}
    regulon_rows: list[dict[str, object]] = []
    key_files: list[Path] = []
    for method in ("genie3", "grnboost2"):
        adjacency_path = ROOT / f"intermediate/adjacencies_{method}_v1.tsv"
        regulon_path = ROOT / f"intermediate/regulons_{method}_v1.dat"
        auc_path = ROOT / f"intermediate/auc_{method}_all_cells_v1.csv"
        key_files.extend([adjacency_path, regulon_path, auc_path])

        adjacency = pd.read_csv(adjacency_path, sep="\t")
        required = ["TF", "target", "importance"]
        check(f"{method}_adjacency_schema", adjacency.columns.tolist() == required, f"columns={adjacency.columns.tolist()}")
        importance = pd.to_numeric(adjacency["importance"], errors="coerce").to_numpy(float)
        duplicates = int(adjacency.duplicated(["TF", "target"]).sum())
        valid_edges = np.isfinite(importance).all() and (importance >= 0).all() and adjacency["TF"].notna().all() and adjacency["target"].notna().all()
        check(f"{method}_adjacency_values", bool(valid_edges), f"rows={len(adjacency)};TFs={adjacency['TF'].nunique()};targets={adjacency['target'].nunique()}")
        check(f"{method}_adjacency_unique_edges", duplicates == 0, f"duplicate_TF_target={duplicates}")
        del adjacency, importance

        with regulon_path.open("rb") as handle:
            regulons = pickle.load(handle)
        names = [str(regulon.name) for regulon in regulons]
        target_counts = np.asarray([len(regulon.gene2weight) for regulon in regulons], dtype=int)
        check(f"{method}_regulon_objects", len(regulons) > 100 and len(set(names)) == len(names), f"regulons={len(regulons)};unique={len(set(names))}")
        check(f"{method}_regulon_target_counts", bool((target_counts > 0).all()), f"min={target_counts.min()};median={np.median(target_counts):.1f};max={target_counts.max()}")
        for name, count in zip(names, target_counts, strict=True):
            regulon_rows.append({"method": method, "regulon": name, "n_targets": int(count)})

        auc = pd.read_csv(auc_path, index_col=0)
        auc.index = auc.index.astype(str)
        values = auc.to_numpy(float)
        check(f"{method}_aucell_shape", auc.shape[0] == 11950 and auc.shape[1] == len(regulons), f"shape={auc.shape};regulons={len(regulons)}")
        check(f"{method}_aucell_cell_order", auc.index.tolist() == cells, f"rows={len(auc)}")
        check(f"{method}_aucell_regulon_names", auc.columns.astype(str).tolist() == names, f"columns={len(auc.columns)}")
        check(f"{method}_aucell_finite_nonnegative", bool(np.isfinite(values).all() and (values >= 0).all()), f"finite={np.isfinite(values).sum()};total={values.size};min={np.nanmin(values):.6g};max={np.nanmax(values):.6g}")
        nonconstant = int((np.nanstd(values, axis=0) > 0).sum())
        check(f"{method}_aucell_nonconstant", nonconstant == auc.shape[1], f"nonconstant={nonconstant}/{auc.shape[1]}")
        auc_by_method[method] = auc

    pd.DataFrame(regulon_rows).to_csv(OUT / "pyscenic_regulon_size_qc_v1.tsv", sep="\t", index=False)
    shared = sorted(set(auc_by_method["genie3"].columns) & set(auc_by_method["grnboost2"].columns))
    concordance = []
    for regulon in shared:
        x = auc_by_method["genie3"][regulon].to_numpy(float)
        y = auc_by_method["grnboost2"][regulon].to_numpy(float)
        rho = float(spearmanr(x, y).statistic) if np.std(x) > 0 and np.std(y) > 0 else np.nan
        concordance.append({"regulon": regulon, "n_cells": len(x), "spearman_rho": rho})
    concordance_frame = pd.DataFrame(concordance)
    concordance_frame["BH_FDR_not_applicable"] = "descriptive_cross_method_QC"
    concordance_frame.to_csv(OUT / "pyscenic_cross_method_regulon_concordance_v1.tsv", sep="\t", index=False)
    finite_rho = concordance_frame["spearman_rho"].dropna().to_numpy(float)
    check("cross_method_shared_regulons", len(shared) > 100, f"shared={len(shared)}")
    check("cross_method_concordance_finite", len(finite_rho) == len(shared), f"finite={len(finite_rho)}/{len(shared)}")

    sampling = pd.read_csv(ROOT / "tables/grn_inference_cells_v1.tsv.gz", sep="\t", dtype=str)
    check("grn_sampling_unique", len(sampling) == 3115 and sampling.iloc[:, 0].nunique() == 3115, f"rows={len(sampling)};unique={sampling.iloc[:,0].nunique()}")
    check("grn_sampling_subset", set(sampling.iloc[:, 0].astype(str)).issubset(set(cells)), f"rows={len(sampling)}")

    manifest = []
    for path in key_files:
        manifest.append({"relative_path": str(path.relative_to(PROJECT)), "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest).to_csv(OUT / "key_file_manifest_v1.tsv", sep="\t", index=False)
    pd.DataFrame(checks).to_csv(OUT / "verification_checks_v1.tsv", sep="\t", index=False)
    failures = [row["check"] for row in checks if row["status"] != "PASS"]
    metadata_out = {
        "schema_version": "part1_pyscenic_terminal_verification_v1",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "n_checks": len(checks),
        "n_failed": len(failures),
        "failed_checks": failures,
        "n_cells": 11950,
        "genie3_regulons": int(auc_by_method["genie3"].shape[1]),
        "grnboost2_regulons": int(auc_by_method["grnboost2"].shape[1]),
        "shared_regulons": len(shared),
        "cross_method_median_spearman_rho": float(np.median(finite_rho)),
        "cross_method_fraction_positive": float(np.mean(finite_rho > 0)),
        "interpretation": "descriptive algorithmic concordance; cells are not biological replicates",
    }
    (OUT / "verification_metadata_v1.json").write_text(json.dumps(metadata_out, indent=2) + "\n", encoding="utf-8")
    if failures:
        raise RuntimeError(f"pySCENIC terminal verification failed: {failures}")
    (OUT / "RUN_COMPLETE").write_text(metadata_out["completed_at"] + "\n", encoding="utf-8")
    print(f"PART1_PYSCENIC_TERMINAL_VERIFICATION_V1_PASS checks={len(checks)} shared={len(shared)}")


if __name__ == "__main__":
    main()
