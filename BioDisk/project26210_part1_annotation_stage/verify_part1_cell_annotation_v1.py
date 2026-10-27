#!/usr/bin/env python3
"""Independent verification for Part 1 broad cell annotation v1."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from PIL import Image
from scipy import sparse


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
SOURCE = PROJECT / "results/part1_rebuild_v1/02_integration_annotation/scvi_integration_v1/part1_scvi_integrated_v1.h5ad"
ANNOTATION_DIR = PROJECT / "results/part1_rebuild_v1/02_integration_annotation/cell_annotation_v1"
TARGET = ANNOTATION_DIR / "part1_scvi_integrated_cell_annotated_v1.h5ad"
VERIFY_DIR = PROJECT / "results/part1_rebuild_v1/02_integration_annotation/cell_annotation_verification_v1"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"


def gate(stage: str, rows: list[dict[str, object]]) -> None:
    proc = subprocess.run([str(MEMORY_GATE)], text=True, capture_output=True)
    output = (proc.stdout or proc.stderr).strip()
    print(f"memory_gate[{stage}]: {output}", flush=True)
    rows.append({"stage": stage, "returncode": proc.returncode, "output": output})
    if proc.returncode:
        raise RuntimeError(output)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def matrix_equal(left, right) -> bool:
    if sparse.issparse(left) or sparse.issparse(right):
        left = sparse.csr_matrix(left)
        right = sparse.csr_matrix(right)
        return left.shape == right.shape and (left != right).nnz == 0
    return np.array_equal(np.asarray(left), np.asarray(right), equal_nan=True)


def check(name: str, condition: bool, detail: str, rows: list[dict[str, object]]) -> None:
    rows.append({"check": name, "passed": bool(condition), "detail": detail})
    if not condition:
        raise AssertionError(f"{name}: {detail}")


def main() -> None:
    if (VERIFY_DIR / "RUN_COMPLETE").exists():
        raise SystemExit(f"Refusing to overwrite completed verification: {VERIFY_DIR}")
    if VERIFY_DIR.exists() and any(VERIFY_DIR.iterdir()):
        raise SystemExit(f"Verification output is non-empty: {VERIFY_DIR}")
    VERIFY_DIR.mkdir(parents=True, exist_ok=True)
    started = datetime.now().astimezone().isoformat()
    gates: list[dict[str, object]] = []
    checks: list[dict[str, object]] = []

    for path in [SOURCE, TARGET, ANNOTATION_DIR / "RUN_COMPLETE", ANNOTATION_DIR / "cluster_annotation_mapping_v1.tsv", ANNOTATION_DIR / "cell_annotation_v1.tsv.gz", ANNOTATION_DIR / "celltype_counts_v1.tsv"]:
        check(f"exists::{path.name}", path.exists(), str(path), checks)

    gate("before_source_open", gates)
    source = ad.read_h5ad(SOURCE, backed="r")
    gate("before_target_open", gates)
    target = ad.read_h5ad(TARGET, backed="r")

    check("shape_exact", source.shape == target.shape == (258120, 37795), f"source={source.shape}; target={target.shape}", checks)
    check("obs_names_exact", source.obs_names.equals(target.obs_names), f"n={source.n_obs}", checks)
    check("var_names_exact", source.var_names.equals(target.var_names), f"n={source.n_vars}", checks)
    check("source_obs_columns_preserved", list(source.obs.columns) == list(target.obs.columns[: len(source.obs.columns)]), f"source_n={len(source.obs.columns)}; target_n={len(target.obs.columns)}", checks)
    try:
        pd.testing.assert_frame_equal(source.obs, target.obs[source.obs.columns], check_dtype=True, check_categorical=True)
        obs_equal = True
        obs_detail = f"all {len(source.obs.columns)} source obs columns exact"
    except AssertionError as exc:
        obs_equal = False
        obs_detail = str(exc)
    check("source_obs_values_exact", obs_equal, obs_detail, checks)

    required_fields = ["cell_type_broad_v1", "annotation_confidence_v1", "annotation_ambiguous_v1", "annotation_evidence_v1", "annotation_cluster_v1"]
    check("annotation_fields_present", all(x in target.obs for x in required_fields), " | ".join(required_fields), checks)
    check("annotation_no_missing", not target.obs[required_fields].isna().any().any(), str(target.obs[required_fields].isna().sum().to_dict()), checks)
    check("ambiguous_exact_cluster23", np.array_equal(target.obs["annotation_ambiguous_v1"].to_numpy(), target.obs["leiden04"].astype(str).eq("23").to_numpy()), f"n={int(target.obs['annotation_ambiguous_v1'].sum())}", checks)
    check("low_confidence_exact_cluster23", target.obs.loc[target.obs["annotation_confidence_v1"].astype(str).eq("Low"), "leiden04"].astype(str).eq("23").all(), f"low_n={int(target.obs['annotation_confidence_v1'].astype(str).eq('Low').sum())}", checks)
    check("legacy_cell_type_not_copied", "cell_type" not in target.obs, "cell_type column absent", checks)

    mapping = pd.read_csv(ANNOTATION_DIR / "cluster_annotation_mapping_v1.tsv", sep="\t", dtype={"leiden04": str})
    check("mapping_27_clusters", len(mapping) == 27 and set(mapping["leiden04"]) == set(target.obs["leiden04"].astype(str).unique()), f"rows={len(mapping)}", checks)
    mapped_labels = target.obs["leiden04"].astype(str).map(mapping.set_index("leiden04")["cell_type_broad_v1"])
    check("cell_labels_match_mapping", mapped_labels.astype(str).equals(target.obs["cell_type_broad_v1"].astype(str)), f"n={target.n_obs}", checks)

    counts_file = pd.read_csv(ANNOTATION_DIR / "celltype_counts_v1.tsv", sep="\t")
    observed_counts = target.obs["cell_type_broad_v1"].astype(str).value_counts().rename_axis("cell_type_broad_v1").reset_index(name="n_cells")
    merged_counts = counts_file.merge(observed_counts, on="cell_type_broad_v1", suffixes=("_file", "_h5ad"), how="outer").fillna(-1)
    check("celltype_counts_match", (merged_counts["n_cells_file"].astype(int) == merged_counts["n_cells_h5ad"].astype(int)).all(), merged_counts.to_json(orient="records"), checks)
    check("celltype_count_sum", int(counts_file["n_cells"].sum()) == target.n_obs, f"sum={int(counts_file['n_cells'].sum())}", checks)
    check("fourteen_nonempty_labels", int((counts_file["n_cells"] > 0).sum()) == 14, f"n={(counts_file['n_cells'] > 0).sum()}", checks)

    cell_table = pd.read_csv(ANNOTATION_DIR / "cell_annotation_v1.tsv.gz", sep="\t", index_col=0)
    check("cell_table_ids_exact", pd.Index(cell_table.index.astype(str)).equals(target.obs_names.astype(str)), f"rows={len(cell_table)}", checks)
    check("cell_table_labels_exact", np.array_equal(cell_table["cell_type_broad_v1"].astype(str).to_numpy(), target.obs["cell_type_broad_v1"].astype(str).to_numpy()), f"rows={len(cell_table)}", checks)

    for start in [0, source.n_obs // 2, source.n_obs - 1024]:
        end = min(start + 1024, source.n_obs)
        check(f"X_exact_{start}_{end}", matrix_equal(source.X[start:end], target.X[start:end]), f"slice={start}:{end}", checks)
        check(f"counts_exact_{start}_{end}", matrix_equal(source.layers["counts"][start:end], target.layers["counts"][start:end]), f"slice={start}:{end}", checks)
    for key in ["X_scVI", "X_pca", "X_umap"]:
        left = np.asarray(source.obsm[key])
        right = np.asarray(target.obsm[key])
        check(f"obsm_exact::{key}", np.array_equal(left, right, equal_nan=True), f"shape={left.shape}", checks)
    for key in ["connectivities", "distances"]:
        left = source.obsp[key]
        right = target.obsp[key]
        check(f"obsp_shape_nnz::{key}", left.shape == right.shape and left.nnz == right.nnz, f"shape={left.shape}; nnz={left.nnz}", checks)

    annotation_uns = target.uns.get("part1_cell_annotation_v1", {})
    check("annotation_uns_present", bool(annotation_uns), str(annotation_uns.keys()), checks)
    check("uns_legacy_copy_false", annotation_uns.get("legacy_labels_copied") in {False, np.bool_(False)}, str(annotation_uns.get("legacy_labels_copied")), checks)
    check("uns_evidence_v2", "annotation_evidence_v2" in str(annotation_uns.get("evidence_directory", "")), str(annotation_uns.get("evidence_directory", "")), checks)

    for stem in ["umap_cell_type_broad_v1", "umap_annotation_confidence_v1", "celltype_counts_v1"]:
        png = ANNOTATION_DIR / f"{stem}.png"
        pdf = ANNOTATION_DIR / f"{stem}.pdf"
        check(f"figure_files::{stem}", png.exists() and pdf.exists() and pdf.stat().st_size > 1000, f"png={png.stat().st_size if png.exists() else 0}; pdf={pdf.stat().st_size if pdf.exists() else 0}", checks)
        with Image.open(png) as image:
            check(f"figure_dimensions::{stem}", image.width >= 1000 and image.height >= 800, f"{image.width}x{image.height}", checks)

    sidecar_hash = (ANNOTATION_DIR / "part1_scvi_integrated_cell_annotated_v1.h5ad.sha256").read_text().split()[0]
    actual_hash = sha256(TARGET)
    check("h5ad_sha256_sidecar", actual_hash == sidecar_hash, actual_hash, checks)
    check("all_checks_passed", all(row["passed"] for row in checks), f"checks={len(checks)}", checks)

    source.file.close()
    target.file.close()
    frame = pd.DataFrame(checks)
    frame.to_csv(VERIFY_DIR / "verification_checks_v1.tsv", sep="\t", index=False)
    metadata = {
        "run_id": "part1_cell_annotation_verification_v1",
        "started_at": started,
        "completed_at": datetime.now().astimezone().isoformat(),
        "python_executable": sys.executable,
        "source_h5ad": str(SOURCE),
        "target_h5ad": str(TARGET),
        "target_sha256": actual_hash,
        "n_checks": len(checks),
        "n_passed": int(frame["passed"].sum()),
        "memory_gate": gates,
    }
    (VERIFY_DIR / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (VERIFY_DIR / "RUN_COMPLETE").write_text(metadata["completed_at"] + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
