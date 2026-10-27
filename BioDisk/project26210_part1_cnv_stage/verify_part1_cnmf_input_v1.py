#!/usr/bin/env python3
"""Independent verification of the fresh strict cNMF input matrix."""

from __future__ import annotations

import gzip
import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1"
H5AD = ROOT / "03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad"
AUDIT = ROOT / "04_subclones_tymp/input_audit_v1"
INPUT = ROOT / "05_cnmf/input_v1"
OUT = ROOT / "05_cnmf/input_verification_v1"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"
EXPECTED_H5AD_SHA256 = "e4779c46af13a8643e8dec584eeff4322258181f5edab3e55dcd8c20bc28a9e7"


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
            pd.DataFrame(checks).to_csv(OUT / "verification_checks_v1.tsv", sep="\t", index=False)
            raise RuntimeError(f"Verification failed: {name}: {evidence}")

    gate = subprocess.run([str(MEMORY_GATE)], text=True, capture_output=True)
    gate_text = " ".join((gate.stdout + " " + gate.stderr).split())
    print(gate_text, flush=True)
    check("memory_gate", gate.returncode == 0 and "PASS" in gate_text, gate_text)
    check("input_complete", (INPUT / "RUN_COMPLETE").is_file(), INPUT / "RUN_COMPLETE")

    cells = pd.read_csv(INPUT / "cnmf_cells_v1.tsv.gz", sep="\t", dtype={"exact_cell_id": str})
    genes = pd.read_csv(INPUT / "cnmf_hvg_genes_v1.tsv", sep="\t", dtype={"gene": str})
    cell_meta = pd.read_csv(INPUT / "cnmf_cell_metadata_v1.tsv.gz", sep="\t", dtype={"exact_cell_id": str})
    gene_qc = pd.read_csv(INPUT / "cnmf_gene_qc_v1.tsv.gz", sep="\t", dtype={"gene": str})
    audit_registry = pd.read_csv(AUDIT / "subclone_cnmf_cell_registry_v1.tsv.gz", sep="\t", dtype={"exact_cell_id": str})

    check("cell_rows", len(cells) == 11950, len(cells))
    check("cell_ids_unique", cells["exact_cell_id"].is_unique, cells["exact_cell_id"].nunique())
    check("cell_indices", np.array_equal(cells["cell_index_1based"].to_numpy(), np.arange(1, 11951)), "1..11950")
    check("hvg_rows", len(genes) == 2500, len(genes))
    check("hvg_unique", genes["gene"].is_unique, genes["gene"].nunique())
    check("hvg_indices", np.array_equal(genes["gene_index_1based"].to_numpy(), np.arange(1, 2501)), "1..2500")
    check("gene_qc_rows", len(gene_qc) == 37795 and gene_qc["gene"].is_unique, [len(gene_qc), gene_qc["gene"].nunique()])

    adata = ad.read_h5ad(H5AD, backed="r")
    eligible = adata.obs["eligible_for_subclone_selection_v1"].astype(bool)
    expected = adata.obs["included_in_primary_malignancy_summary_v1"].astype(bool) & adata.obs["malignant_conservative_v1"].astype(bool)
    expected_ids = adata.obs_names[eligible].astype(str).to_numpy()
    selected_obs = adata.obs.loc[eligible]
    check("h5ad_shape", adata.shape == (258120, 37795), adata.shape)
    check("h5ad_selection_rule_exact", eligible.equals(expected), int((eligible != expected).sum()))
    check("h5ad_eligible_count", int(eligible.sum()) == 11950, int(eligible.sum()))
    check("cell_order_matches_h5ad", np.array_equal(cells["exact_cell_id"].to_numpy(str), expected_ids), "exact")
    check("selected_all_tumor", selected_obs["source_tissue"].astype(str).eq("Tumor").all(), len(selected_obs))
    check("selected_all_clear_epithelial", selected_obs["cell_type_broad_v1"].astype(str).eq("Epithelial cell").all(), len(selected_obs))
    mixed_selected = adata.obs.loc[adata.obs["cell_type_broad_v1"].astype(str).eq("Epithelial/Myeloid mixed"), "eligible_for_subclone_selection_v1"].astype(bool)
    check("mixed_cells_excluded", not mixed_selected.any(), int(mixed_selected.sum()))
    check("hvg_genes_exist_in_h5ad", set(genes["gene"]).issubset(set(adata.var_names.astype(str))), len(genes))
    adata.file.close()

    check("cell_order_matches_audit", np.array_equal(cells["exact_cell_id"].to_numpy(str), audit_registry["exact_cell_id"].to_numpy(str)), "exact")
    check("metadata_rows", len(cell_meta) == 11950 and cell_meta["exact_cell_id"].is_unique, [len(cell_meta), cell_meta["exact_cell_id"].nunique()])
    check("metadata_order", np.array_equal(cell_meta["exact_cell_id"].to_numpy(str), cells["exact_cell_id"].to_numpy(str)), "exact")
    status_counts = cell_meta["malignancy_status_v1"].value_counts().to_dict()
    allowed_status = {"malignant_consensus", "malignant_infercnv_only"}
    check("metadata_allowed_malignancy_status", set(status_counts).issubset(allowed_status), status_counts)
    check("metadata_malignancy_status_counts", status_counts == {"malignant_consensus": 10844, "malignant_infercnv_only": 1106}, status_counts)
    check("metadata_all_conservative_malignant", cell_meta["malignant_conservative_v1"].astype(bool).all(), int(cell_meta["malignant_conservative_v1"].astype(bool).sum()))
    check("metadata_all_primary_included", cell_meta["included_in_primary_malignancy_summary_v1"].astype(bool).all(), int(cell_meta["included_in_primary_malignancy_summary_v1"].astype(bool).sum()))
    check("metadata_all_subclone_eligible", cell_meta["eligible_for_subclone_selection_v1"].astype(bool).all(), int(cell_meta["eligible_for_subclone_selection_v1"].astype(bool).sum()))
    burden = pd.to_numeric(cell_meta["infercnv_cnv_gene_burden"], errors="coerce")
    threshold = pd.to_numeric(cell_meta["infercnv_high_cnv_threshold"], errors="coerce")
    ratio = pd.to_numeric(cell_meta["infercnv_burden_over_threshold_v1"], errors="coerce")
    check("metadata_cnv_complete", burden.notna().all() and threshold.notna().all() and ratio.notna().all(), [burden.isna().sum(), threshold.isna().sum(), ratio.isna().sum()])
    check("metadata_cnv_ratio", np.allclose(ratio, burden / threshold, rtol=1e-12, atol=1e-12), float(np.max(np.abs(ratio - burden / threshold))))
    tymp = pd.to_numeric(cell_meta["TYMP_log1p_CP10K_v1"], errors="coerce")
    check("metadata_tymp_finite_nonnegative", np.isfinite(tymp).all() and tymp.ge(0).all(), [tymp.min(), tymp.max()])
    check("metadata_tymp_positive", int(tymp.gt(0).sum()) == 8178, int(tymp.gt(0).sum()))

    selected_qc = gene_qc.loc[gene_qc["selected_hvg"].astype(bool)].copy()
    selected_qc = selected_qc.sort_values("dispersion_rank_among_eligible", kind="stable")
    check("gene_qc_selected_count", len(selected_qc) == 2500, len(selected_qc))
    check("hvg_order_matches_qc_rank", np.array_equal(genes["gene"].to_numpy(str), selected_qc["gene"].to_numpy(str)), "exact")
    check("selected_hvg_eligible", selected_qc["eligible"].astype(bool).all(), int(selected_qc["eligible"].astype(bool).sum()))
    check("selected_hvg_detection_floor", selected_qc["detected_cells"].ge(119).all(), int(selected_qc["detected_cells"].min()))
    upper = genes["gene"].str.upper()
    mito_unprefixed = {"COX1", "COX2", "COX3", "CYTB", "ATP6", "ATP8", "ND1", "ND2", "ND3", "ND4", "ND5", "ND6", "ND4L"}
    technical = (
        upper.str.startswith("MT-") | upper.str.startswith("MTRNR") |
        upper.str.startswith("RPL") | upper.str.startswith("RPS") |
        upper.str.startswith("HBA") | upper.str.startswith("HBB") |
        upper.str.startswith("HBM") | upper.str.startswith("HBD") |
        upper.str.startswith("HBG") | upper.isin(mito_unprefixed)
    )
    check("selected_hvg_no_technical_genes", not technical.any(), genes.loc[technical, "gene"].tolist())
    check("tymp_not_selected_hvg", "TYMP" not in set(genes["gene"]), "TYMP excluded from factorization matrix; retained in metadata")

    with gzip.open(INPUT / "cnmf_hvg_scaled_genes_by_cells_v1.mtx.gz", "rb") as handle:
        matrix = mmread(handle)
    check("matrix_sparse", sp.issparse(matrix), type(matrix).__name__)
    matrix = matrix.tocsr()
    check("matrix_shape", matrix.shape == (2500, 11950), matrix.shape)
    check("matrix_nnz", matrix.nnz == 2416110, matrix.nnz)
    check("matrix_finite", np.isfinite(matrix.data).all(), int(np.count_nonzero(~np.isfinite(matrix.data))))
    check("matrix_nonnegative", np.all(matrix.data >= 0), float(matrix.data.min()))
    check("matrix_scale_cap", float(matrix.data.max()) <= 10.0 + 1e-12, float(matrix.data.max()))
    check("matrix_nonempty_rows", np.all(np.diff(matrix.indptr) > 0), int(np.count_nonzero(np.diff(matrix.indptr) == 0)))

    qc = pd.read_csv(INPUT / "cnmf_input_qc_v1.tsv", sep="\t").set_index("metric")["value"]
    check("qc_matrix_shape", int(float(qc["scaled_rows"])) == matrix.shape[0] and int(float(qc["scaled_columns"])) == matrix.shape[1], [qc["scaled_rows"], qc["scaled_columns"]])
    check("qc_matrix_nnz", int(float(qc["scaled_nnz"])) == matrix.nnz, qc["scaled_nnz"])
    check("qc_hvg_count", int(float(qc["n_hvg"])) == len(genes), qc["n_hvg"])
    check("qc_tymp_positive", int(float(qc["TYMP_positive_cells"])) == int(tymp.gt(0).sum()), qc["TYMP_positive_cells"])
    check("qc_invalid_matrix_zero", int(float(qc["scaled_nonfinite"])) == 0 and int(float(qc["scaled_negative"])) == 0, [qc["scaled_nonfinite"], qc["scaled_negative"]])
    check("qc_missing_cnv_zero", int(float(qc["metadata_missing_cnv_burden"])) == 0, qc["metadata_missing_cnv_burden"])

    metadata = json.loads((INPUT / "run_metadata_v1.json").read_text(encoding="utf-8"))
    check("input_hash_frozen", metadata["input_h5ad_sha256"] == EXPECTED_H5AD_SHA256, metadata["input_h5ad_sha256"])
    check("selection_metadata", metadata["selection_rule"] == "eligible_for_subclone_selection_v1 == True", metadata["selection_rule"])
    check("orientation_metadata", metadata["matrix_orientation"] == "genes_by_cells", metadata["matrix_orientation"])
    check("no_dependency_changes", metadata["dependency_changes"] is False, metadata["dependency_changes"])

    manifest = pd.read_csv(INPUT / "file_manifest_v1.tsv", sep="\t", dtype=str)
    check("manifest_unique", manifest["file"].is_unique, manifest["file"].nunique())
    for row in manifest.itertuples(index=False):
        path = INPUT / row.file
        check(f"hash_{row.file}", path.is_file() and path.stat().st_size == int(row.size_bytes) and sha256(path) == row.sha256, row.sha256)

    pd.DataFrame(checks).to_csv(OUT / "verification_checks_v1.tsv", sep="\t", index=False)
    run_metadata = {
        "schema_version": "part1_cnmf_input_verification_v1",
        "completed_at": datetime.now().astimezone().isoformat(),
        "n_checks": len(checks),
        "n_pass": sum(row["status"] == "PASS" for row in checks),
        "n_fail": sum(row["status"] == "FAIL" for row in checks),
        "dependency_changes": False,
    }
    (OUT / "run_metadata_v1.json").write_text(json.dumps(run_metadata, indent=2) + "\n", encoding="utf-8")
    (OUT / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n", encoding="utf-8")
    print(json.dumps(run_metadata, indent=2))
    print("PART1_CNMF_INPUT_VERIFICATION_V1_PASS")


if __name__ == "__main__":
    main()
