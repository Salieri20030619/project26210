#!/usr/bin/env python3
from __future__ import annotations

import gzip
import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
STAGE_V = Path("/mnt/BioDisk/wangmingyuan/project26210_tymp_vascular_niche_stage_v1")
STAGE_T = Path("/mnt/BioDisk/wangmingyuan/project26210_tymp_overnight_stage_v1")
VASCULAR = STAGE_V / "output_v1"
TRAJECTORY = STAGE_T / "trajectory_v1"
QC = STAGE_T / "joint_qc_v2"
MEMORY_GATE = ROOT / "check_memory_gate.sh"


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    result = np.empty_like(ranked)
    result[order] = np.minimum(ranked, 1.0)
    return result


def main() -> None:
    if QC.exists():
        raise SystemExit(f"Refusing to overwrite existing QC directory: {QC}")
    QC.mkdir(parents=True)
    gate = subprocess.run([str(MEMORY_GATE)], text=True, capture_output=True, check=False)
    gate_line = (gate.stdout or gate.stderr).strip()
    print(f"[memory_gate] {gate_line}", flush=True)
    if gate.returncode != 0:
        raise RuntimeError(gate_line)

    checks: list[dict] = []

    def check(section: str, item: str, passed: bool, observed: object, expected: object) -> None:
        checks.append({
            "section": section,
            "check": item,
            "status": "PASS" if passed else "FAIL",
            "observed": str(observed),
            "expected": str(expected),
        })

    # File-level completion and metadata.
    for section, folder in [("vascular", VASCULAR), ("trajectory", TRAJECTORY)]:
        check(section, "RUN_COMPLETE exists", (folder / "RUN_COMPLETE").is_file(),
              (folder / "RUN_COMPLETE").is_file(), True)
        metadata = json.loads((folder / "run_metadata_v1.json").read_text())
        check(section, "metadata completed", metadata.get("status") == "completed",
              metadata.get("status"), "completed")

    vt = VASCULAR / "tables"
    tt = TRAJECTORY / "tables"

    # Vascular analysis structural checks.
    fingerprints = pd.read_csv(vt / "input_fingerprints_v1.tsv", sep="\t")
    gates = pd.read_csv(vt / "memory_gates_v1.tsv", sep="\t")
    block = pd.read_csv(vt / "sample_spatial_block_tests_v1.tsv", sep="\t")
    cohort = pd.read_csv(vt / "tjrcc_cohort_direction_summary_v1.tsv", sep="\t")
    contrasts = pd.read_csv(vt / "sample_interface_contrasts_v1.tsv", sep="\t")
    interface = pd.read_csv(vt / "tjrcc_interface_direction_summary_v1.tsv", sep="\t")
    primary = pd.read_csv(vt / "primary_5x5_mechanism_summary_v1.tsv", sep="\t")
    profiles = pd.read_csv(vt / "boundary_distance_profiles_v1.tsv", sep="\t")
    gene_qc = pd.read_csv(vt / "signature_gene_qc_v1.tsv", sep="\t")
    clinical = pd.read_csv(vt / "clinical_resistance_label_audit_v1.tsv", sep="\t")
    spot = pd.read_csv(vt / "spot_mechanism_scores_v1.tsv.gz", sep="\t", compression="gzip")

    check("vascular", "14 input samples", len(fingerprints) == 14 and fingerprints["sample"].nunique() == 14,
          f"rows={len(fingerprints)}; samples={fingerprints['sample'].nunique()}", "14;14")
    check("vascular", "spot total", len(spot) == 41688, len(spot), 41688)
    check("vascular", "cohort spot totals", spot.groupby("cohort").size().to_dict() == {"GLMF": 7079, "TJ-RCC": 34609},
          spot.groupby("cohort").size().to_dict(), {"GLMF": 7079, "TJ-RCC": 34609})
    check("vascular", "memory gates", len(gates) == 16 and (gates["exit_code"] == 0).all() and gates["line"].str.contains("PASS").all(),
          f"rows={len(gates)}; pass={(gates['exit_code'] == 0).sum()}", "16;16")
    per_sample_counts = block.groupby(["cohort", "sample"]).size()
    correct_counts = (
        (per_sample_counts.loc["GLMF"] == 14).all()
        and (per_sample_counts.loc["TJ-RCC"] == 18).all()
    )
    check("vascular", "block-test completeness", len(block) == 244 and correct_counts,
          f"rows={len(block)}; per_sample={per_sample_counts.to_dict()}",
          "2 GLMF x 14 plus 12 TJ-RCC x 18 = 244")
    actual_combos = block[["sample", "relation", "route", "grid"]]
    check("vascular", "block-test unique combinations", not actual_combos.duplicated().any() and len(actual_combos) == 244,
          f"rows={len(actual_combos)}; duplicates={actual_combos.duplicated().sum()}", "244;0")
    numeric_block = ["n_spots", "n_blocks", "rho", "direction_aligned_rho", "p_value", "permutations", "p_fdr_bh_within_sample"]
    check("vascular", "finite block statistics", np.isfinite(block[numeric_block].to_numpy(float)).all(),
          int(np.isfinite(block[numeric_block].to_numpy(float)).sum()), block[numeric_block].size)
    check("vascular", "9999 permutations", (block["permutations"] == 9999).all(),
          sorted(block["permutations"].unique()), [9999])
    check("vascular", "signature coverage", len(gene_qc) == 42 and (gene_qc["n_used"] >= 5).all(),
          f"rows={len(gene_qc)}; min_n={gene_qc['n_used'].min()}", "42;>=5")
    check("vascular", "boundary profile completeness", len(profiles) == 14 * 8 * 5,
          len(profiles), 560)
    check("vascular", "clinical outcome unavailable", len(clinical) == 1 and clinical.loc[0, "status"] == "not_estimable" and
          clinical.loc[0, "patients_with_treatment_response"] == 0 and clinical.loc[0, "patients_with_survival"] == 0,
          clinical.to_dict("records"), "12 patients; no response/survival; not_estimable")

    # Recompute all patient-level direction summaries and global BH.
    tj = block.loc[block["cohort"].eq("TJ-RCC")].copy()
    recomputed = []
    for keys, frame in tj.groupby(["relation", "route", "grid", "expected_direction"], sort=True):
        direction = 1 if keys[3] == "+" else -1
        aligned = frame["rho"].to_numpy(float) * direction
        positive = int((aligned > 0).sum())
        recomputed.append({"keys": keys, "n": len(frame), "positive": positive,
                           "median_rho": float(frame["rho"].median()),
                           "median_aligned": float(np.median(aligned)),
                           "p": float(binomtest(positive, len(frame), 0.5).pvalue)})
    p_adj = bh(np.array([row["p"] for row in recomputed]))
    cohort_key = cohort.set_index(["relation", "route", "grid", "expected_direction"])
    cohort_ok = len(recomputed) == len(cohort) == 18
    for row, fdr in zip(recomputed, p_adj):
        observed = cohort_key.loc[row["keys"]]
        cohort_ok &= (
            int(observed["n_patients"]) == row["n"]
            and int(observed["direction_supported_patients"]) == row["positive"]
            and np.isclose(observed["median_rho"], row["median_rho"], atol=1e-12)
            and np.isclose(observed["median_direction_aligned_rho"], row["median_aligned"], atol=1e-12)
            and np.isclose(observed["sign_test_p"], row["p"], atol=1e-12)
            and np.isclose(observed["sign_test_fdr_bh_global"], fdr, atol=1e-12)
        )
    check("vascular", "cohort summary recomputation", cohort_ok, cohort_ok, True)

    # Recompute interface direction summaries.
    tjc = contrasts.loc[contrasts["cohort"].eq("TJ-RCC")].copy()
    expected_direction = {
        "SPP1_monocyte_balance": 1,
        "TYMP_MTA_core_reference_dependent": 1,
        "TYMP_MTA_core_reference_light": 1,
        "T_cell_balance": -1,
        "cytotoxicity_residual": -1,
        "hypoxia_angiogenesis_residual": 1,
        "invasion_EMT_residual": 1,
        "spp1_vascular_niche": 1,
    }
    recomputed_i = []
    for (score, contrast), frame in tjc.groupby(["score", "contrast"], sort=True):
        direction = expected_direction[score]
        aligned = frame["difference_z"].to_numpy(float) * direction
        positive = int((aligned > 0).sum())
        recomputed_i.append({"keys": (score, contrast), "n": len(frame), "positive": positive,
                             "median": float(frame["difference_z"].median()),
                             "median_aligned": float(np.median(aligned)),
                             "p": float(binomtest(positive, len(frame), 0.5).pvalue)})
    p_adj_i = bh(np.array([row["p"] for row in recomputed_i]))
    interface_key = interface.set_index(["score", "contrast"])
    interface_ok = len(recomputed_i) == len(interface) == 16
    for row, fdr in zip(recomputed_i, p_adj_i):
        observed = interface_key.loc[row["keys"]]
        interface_ok &= (
            int(observed["n_patients"]) == row["n"]
            and int(observed["direction_supported_patients"]) == row["positive"]
            and np.isclose(observed["median_difference_z"], row["median"], atol=1e-12)
            and np.isclose(observed["median_direction_aligned_difference_z"], row["median_aligned"], atol=1e-12)
            and np.isclose(observed["sign_test_p"], row["p"], atol=1e-12)
            and np.isclose(observed["sign_test_fdr_bh_global"], fdr, atol=1e-12)
        )
    check("vascular", "interface summary recomputation", interface_ok, interface_ok, True)
    check("vascular", "primary summary is seven 5x5 rows", len(primary) == 7 and primary["grid"].eq("5x5").all(),
          f"rows={len(primary)}; grids={primary['grid'].unique().tolist()}", "7;['5x5']")

    # Verify current input fingerprints. This is deliberately independent and can be I/O intensive.
    fingerprint_ok = True
    fingerprint_rows = []
    for row in fingerprints.itertuples(index=False):
        path = Path(row.path)
        observed_size = path.stat().st_size if path.is_file() else -1
        observed_hash = sha256(path) if path.is_file() else "MISSING"
        ok = observed_size == int(row.size_bytes) and observed_hash == row.sha256
        fingerprint_ok &= ok
        fingerprint_rows.append({"sample": row.sample, "path": row.path, "size": observed_size,
                                 "sha256": observed_hash, "status": "PASS" if ok else "FAIL"})
    pd.DataFrame(fingerprint_rows).to_csv(QC / "input_fingerprint_recheck_v1.tsv", sep="\t", index=False)
    check("vascular", "input fingerprints unchanged", fingerprint_ok,
          f"pass={sum(row['status'] == 'PASS' for row in fingerprint_rows)}/{len(fingerprint_rows)}", "14/14")

    # Trajectory checks and independent statistic recomputation.
    tmeta = json.loads((TRAJECTORY / "run_metadata_v1.json").read_text())
    cells = pd.read_csv(tt / "myeloid_trajectory_cell_scores_v1.tsv.gz", sep="\t", compression="gzip")
    assoc = pd.read_csv(tt / "sample_pseudotime_associations_v1.tsv", sep="\t")
    tcohort = pd.read_csv(tt / "cohort_pseudotime_direction_summary_v1.tsv", sep="\t")
    tsummary = pd.read_csv(tt / "celltype_trajectory_summary_v1.tsv", sep="\t")
    transitions = pd.read_csv(tt / "trajectory_neighbor_transition_matrix_v1.tsv", sep="\t")
    tgates = pd.read_csv(tt / "memory_gates_v1.tsv", sep="\t")
    tgenes = pd.read_csv(tt / "signature_gene_qc_v1.tsv", sep="\t")
    score_names = ["Monocyte_program", "SPP1_state_no_SPP1", "C1QC_state", "Vascular_remodeling",
                   "Immunosuppression", "MTA_context", "SPP1_vs_C1QC"]
    numeric_cells = ["dpt_pseudotime", *score_names]

    check("trajectory", "cell count", len(cells) == 53794, len(cells), 53794)
    check("trajectory", "unique cell IDs", cells["cell_id"].is_unique, cells["cell_id"].nunique(), len(cells))
    check("trajectory", "finite cell scores", np.isfinite(cells[numeric_cells].to_numpy(float)).all(),
          int(np.isfinite(cells[numeric_cells].to_numpy(float)).sum()), cells[numeric_cells].size)
    check("trajectory", "pseudotime bounds", cells["dpt_pseudotime"].between(0, 1).all(),
          f"{cells['dpt_pseudotime'].min():.6g}..{cells['dpt_pseudotime'].max():.6g}", "0..1")
    root_rows = cells.loc[cells["cell_id"].eq(tmeta["root_cell"])]
    check("trajectory", "root cell is monocyte at zero", len(root_rows) == 1 and
          root_rows.iloc[0]["trajectory_celltype"] == "Monocyte" and np.isclose(root_rows.iloc[0]["dpt_pseudotime"], 0),
          root_rows[["cell_id", "trajectory_celltype", "dpt_pseudotime"]].to_dict("records"), "one Monocyte; dpt=0")
    counts = cells.groupby("trajectory_celltype").size().sort_index()
    saved_counts = tsummary.set_index("trajectory_celltype")["n_cells"].sort_index()
    check("trajectory", "cell-type count reconciliation", counts.equals(saved_counts), counts.to_dict(), saved_counts.to_dict())
    check("trajectory", "sample association completeness", len(assoc) == 63 * 7 and assoc["sample_id"].nunique() == 63,
          f"rows={len(assoc)}; samples={assoc['sample_id'].nunique()}", "441;63")
    check("trajectory", "memory gates", len(tgates) == 3 and tgates["line"].str.contains("PASS").all(),
          f"rows={len(tgates)}; pass={tgates['line'].str.contains('PASS').sum()}", "3;3")
    check("trajectory", "signature genes", len(tgenes) == 6 and (tgenes["n_used"] >= 5).sum() == 5 and tgenes["n_used"].min() == 3,
          tgenes.set_index("signature")["n_used"].to_dict(), "six signatures; C1QC=3; others>=5")
    check("trajectory", "transition matrix shape", len(transitions) == 16,
          len(transitions), "4x4=16")
    row_sums = transitions.groupby("source")["source_normalized_fraction"].sum()
    check("trajectory", "transition row normalization", np.allclose(row_sums.to_numpy(), 1, atol=1e-7),
          row_sums.to_dict(), "all 1")

    # Recompute within-sample Spearman and BH from the cell-level table.
    assoc_key = assoc.set_index(["sample_id", "score"])
    assoc_ok = True
    eligible = 0
    for sample, frame in cells.groupby("sample_id", sort=True):
        if len(frame) < 50 or frame["trajectory_celltype"].nunique() < 2:
            continue
        eligible += 1
        pvals = []
        stats = []
        for score in score_names:
            result = spearmanr(frame["dpt_pseudotime"], frame[score], nan_policy="omit")
            stats.append((score, float(result.statistic), float(result.pvalue)))
            pvals.append(float(result.pvalue))
        fdr = bh(np.array(pvals))
        for (score, rho, pvalue), adjusted in zip(stats, fdr):
            observed = assoc_key.loc[(sample, score)]
            assoc_ok &= (
                int(observed["n_cells"]) == len(frame)
                and np.isclose(observed["rho"], rho, atol=1e-12)
                and np.isclose(observed["p_value"], pvalue, atol=1e-12)
                and np.isclose(observed["p_fdr_bh_within_sample"], adjusted, atol=1e-12)
            )
    check("trajectory", "sample Spearman recomputation", assoc_ok and eligible == 63,
          f"eligible={eligible}; exact={assoc_ok}", "eligible=63; exact=True")

    # Recompute cohort direction tests and BH.
    recomputed_t = []
    for score, frame in assoc.groupby("score", sort=True):
        positive = int((frame["rho"] > 0).sum())
        recomputed_t.append({"score": score, "n": len(frame), "positive": positive,
                             "median": float(frame["rho"].median()),
                             "p": float(binomtest(positive, len(frame), 0.5).pvalue)})
    t_fdr = bh(np.array([row["p"] for row in recomputed_t]))
    t_key = tcohort.set_index("score")
    t_ok = len(recomputed_t) == len(tcohort) == 7
    for row, adjusted in zip(recomputed_t, t_fdr):
        observed = t_key.loc[row["score"]]
        t_ok &= (
            int(observed["n_samples"]) == row["n"]
            and int(observed["positive_samples"]) == row["positive"]
            and np.isclose(observed["median_rho"], row["median"], atol=1e-12)
            and np.isclose(observed["sign_test_p"], row["p"], atol=1e-12)
            and np.isclose(observed["sign_test_fdr_bh"], adjusted, atol=1e-12)
        )
    check("trajectory", "cohort sign-test recomputation", t_ok, t_ok, True)
    input_path = Path(tmeta["input"])
    input_hash = sha256(input_path)
    check("trajectory", "reference input fingerprint unchanged", input_hash == tmeta["input_sha256"],
          input_hash, tmeta["input_sha256"])

    # Interpretation guardrails are QC items because they prevent overclaiming.
    tc = tcohort.set_index("score")
    check("interpretation", "SPP1 state is not a universal trajectory", int(tc.loc["SPP1_state_no_SPP1", "positive_samples"]) < 63,
          f"{int(tc.loc['SPP1_state_no_SPP1', 'positive_samples'])}/63 positive",
          "retain heterogeneity; do not call lineage conversion")
    check("interpretation", "clinical resistance remains non-estimable", clinical.loc[0, "status"] == "not_estimable",
          clinical.loc[0, "status"], "not_estimable")

    result = pd.DataFrame(checks)
    result.to_csv(QC / "independent_verification_v2.tsv", sep="\t", index=False)

    # Hash every staged output and QC source, excluding this QC directory itself.
    manifest_rows = []
    for base, label in [(VASCULAR, "vascular"), (TRAJECTORY, "trajectory")]:
        for path in sorted(p for p in base.rglob("*") if p.is_file()):
            manifest_rows.append({"section": label, "path": str(path), "size_bytes": path.stat().st_size,
                                  "sha256": sha256(path)})
    pd.DataFrame(manifest_rows).to_csv(QC / "file_manifest_v2.tsv", sep="\t", index=False)

    n_pass = int(result["status"].eq("PASS").sum())
    n_fail = int(result["status"].eq("FAIL").sum())
    summary = {
        "status": "passed" if n_fail == 0 else "failed",
        "completed_at": now(),
        "memory_gate": gate_line,
        "checks_passed": n_pass,
        "checks_failed": n_fail,
        "files_hashed": len(manifest_rows),
        "important_interpretation": {
            "vascular": "SPP1 vascular niche is inversely associated with T-cell balance in 12/12 TJ-RCC patients at both grids; this is spatial association, not therapy resistance or causality.",
            "trajectory": "DPT is cross-sectional state ordering. C1QC rises in 57/63 samples, whereas the SPP1-state score rises in 41/63; this does not prove monocyte-to-SPP1 TAM lineage conversion.",
        },
    }
    (QC / "independent_verification_summary_v2.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
    if n_fail:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
