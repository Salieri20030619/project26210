#!/usr/bin/env python3
"""Independent structural, numerical, provenance, and figure QC for TJ-RCC advanced_full_v2."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "external/ccRCC_spatial/TJ_RCC_figshare_24599295/advanced_full_v2"
GATE = PROJECT / "check_memory_gate.sh"
SAMPLES = ("H46_T", "J38_T", "R114_T", "R29_T", "R51_T", "S15_T", "X49_T", "X98_T", "Y12_T", "Y27_T", "Y7_T", "Z43_T")
checks: list[dict[str, object]] = []
memory_records: list[str] = []


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def gate(stage: str) -> None:
    result = subprocess.run([str(GATE)], check=True, capture_output=True, text=True)
    line = result.stdout.strip()
    memory_records.append(f"{stage}: {line}")
    print(f"MEMORY_GATE stage={stage} {line}", flush=True)


def record(name: str, passed: bool, detail: str) -> None:
    checks.append({"check": name, "passed": bool(passed), "detail": str(detail)})


def read_json(path: Path) -> dict:
    record(f"exists::{path.relative_to(ROOT)}", path.is_file() and path.stat().st_size > 10, f"bytes={path.stat().st_size if path.exists() else 0}")
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> pd.DataFrame:
    record(f"exists::{path.relative_to(ROOT)}", path.is_file() and path.stat().st_size > 10, f"bytes={path.stat().st_size if path.exists() else 0}")
    return pd.read_csv(path, sep="\t")


def read_optional_tsv(path: Path) -> pd.DataFrame:
    record(f"exists::{path.relative_to(ROOT)}", path.is_file() and path.stat().st_size > 0, f"bytes={path.stat().st_size if path.exists() else 0}")
    try:
        return pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def bools(series: pd.Series) -> np.ndarray:
    if pd.api.types.is_bool_dtype(series):
        return series.to_numpy(bool)
    return series.astype(str).str.lower().map({"true": True, "false": False, "1": True, "0": False}).to_numpy(bool)


def verify_image(path: Path) -> None:
    if not path.is_file():
        record(f"image::{path.relative_to(ROOT)}", False, "missing")
        return
    with Image.open(path) as image:
        array = np.asarray(image.convert("RGB").resize((128, 128)), dtype=float)
        passed = image.width >= 800 and image.height >= 600 and float(array.std()) > 3
        record(f"image::{path.relative_to(ROOT)}", passed, f"{image.width}x{image.height};rgb_sd={array.std():.3f}")


def verify_bayesspace(sample: str) -> int:
    directory = ROOT / "bayesspace" / sample
    status = read_json(directory / "status.json")
    bayes = read_json(directory / "bayesspace_run_metadata.json")
    boundary_meta = read_json(directory / "boundary_run_metadata.json")
    record(f"{sample}::bayes_status", status.get("status") == "completed" and bayes.get("status") == "bayesspace_complete" and boundary_meta.get("status") == "completed", json.dumps(status))
    q = read_tsv(directory / f"{sample}_q_tuning.tsv")
    assignments = read_tsv(directory / f"{sample}_bayesspace_assignments.tsv")
    spots = read_tsv(directory / f"{sample}_computational_boundary_spots.tsv")
    domains = read_tsv(directory / f"{sample}_domain_summary.tsv")
    n_spots = int(bayes["n_spots"])
    selected_q = int(bayes["selected_q"])
    q_ok = len(q) == 6 and q["q"].astype(int).tolist() == list(range(3, 9)) and int(q["selected"].sum()) == 1 and int(q.loc[q["selected"], "q"].iloc[0]) == selected_q
    record(f"{sample}::q_tuning", q_ok, f"selected={selected_q};rows={len(q)}")
    structure_ok = len(assignments) == len(spots) == n_spots and assignments["barcode"].is_unique and spots["barcode"].is_unique and assignments["bayesspace_cluster"].nunique() == selected_q and len(domains) == selected_q
    record(f"{sample}::bayes_structure", structure_ok, f"spots={len(spots)};domains={len(domains)}")
    tumor = bools(spots["computational_tumor"]); boundary = bools(spots["computational_boundary_spot"]); distance = spots["signed_distance_um"].to_numpy(float)
    boundary_ok = tumor.any() and (~tumor).any() and boundary.any() and np.isfinite(distance).all() and (distance[tumor] <= 1e-8).all() and (distance[~tumor] >= -1e-8).all() and int(boundary.sum()) == int(boundary_meta["n_boundary_spots"])
    record(f"{sample}::computational_boundary", boundary_ok, f"tumor={tumor.sum()};nontumor={(~tumor).sum()};boundary={boundary.sum()}")
    for stem in ("q_tuning", "bayesspace_domains", "bayesspace_boundary_overview", "boundary_profile_summary"):
        verify_image(directory / f"{sample}_{stem}.png")
    return n_spots


def verify_spatial_programs(sample: str, expected_spots: int) -> None:
    directory = ROOT / "spatial_programs" / sample
    status = read_json(directory / "status.json")
    spark = read_json(directory / "sparkx_run_metadata.json")
    hot = read_json(directory / "hotspot_run_metadata.json")
    record(f"{sample}::spatial_status", status.get("status") == "completed" and spark.get("status") == "completed" and hot.get("status") == "completed", json.dumps(status))
    genes = read_tsv(directory / f"{sample}_sparkx_hotspot_gene_statistics.tsv")
    memberships = read_tsv(directory / f"{sample}_hotspot_module_membership.tsv")
    scores = read_tsv(directory / f"{sample}_hotspot_module_scores.tsv")
    modules = read_optional_tsv(directory / f"{sample}_hotspot_module_summary.tsv")
    numeric = ["sparkx_combined_p", "sparkx_by_fdr", "hotspot_C", "hotspot_Z", "hotspot_fdr", "consensus_rank"]
    gene_ok = len(genes) >= 1000 and genes["gene"].is_unique and np.isfinite(genes[numeric].to_numpy(float)).all() and set(genes["support_class"]) <= {"both_strict", "sparkx_only", "hotspot_only", "not_strict"}
    gene_ok = gene_ok and len(genes) == int(hot["n_genes_tested"]) == int(spark["n_genes_tested"])
    record(f"{sample}::spatial_gene_table", gene_ok, f"genes={len(genes)};both={(genes.support_class == 'both_strict').sum()}")
    module_ok = memberships["gene"].is_unique and set(memberships["gene"]) <= set(genes["gene"]) and len(scores) == expected_spots and scores["barcode"].is_unique
    module_ok = module_ok and len(modules) == int(hot["n_modules"])
    record(f"{sample}::hotspot_modules", module_ok, f"memberships={len(memberships)};modules={len(modules)};scores={len(scores)}")
    for stem in ("sparkx_hotspot_concordance", "hotspot_local_correlation_modules", "hotspot_module_spatial_maps", "hotspot_module_summary"):
        verify_image(directory / f"{sample}_{stem}.png")


def verify_misty(sample: str, expected_spots: int) -> None:
    directory = ROOT / "misty" / "samples" / sample
    status = read_json(directory / "status.json")
    meta = read_json(directory / "run_metadata_v1.json")
    record(f"{sample}::misty_status", status.get("status") == "complete" and meta.get("status") == "complete" and int(meta.get("n_spots", -1)) == expected_spots and int(meta.get("n_features", -1)) == 20, json.dumps(status))
    improvement = read_tsv(directory / f"{sample}_misty_improvement_stats_v1.tsv")
    contribution = read_tsv(directory / f"{sample}_misty_contribution_stats_v1.tsv")
    importance = read_tsv(directory / f"{sample}_misty_importances_v1.tsv")
    targets = set(meta["feature_ids"])
    improvement_ok = set(improvement["target"]) == targets and {"gain.R2", "gain.RMSE"} <= set(improvement["measure"]) and np.isfinite(improvement["mean"]).all()
    record(f"{sample}::misty_performance", improvement_ok, f"rows={len(improvement)};targets={improvement.target.nunique()}")
    contribution_ok = set(contribution["target"]) == targets and set(contribution["view"]) == {"intra", "juxta.110", "para.200"} and np.isfinite(contribution["fraction"]).all()
    record(f"{sample}::misty_contributions", contribution_ok, f"rows={len(contribution)};views={sorted(contribution.view.unique())}")
    expected_views = {"intra", "juxta.110", "para.200"}
    structural_self = importance["Predictor"].eq(importance["Target"])
    missing_importance = importance["Importance"].isna()
    importance_grid_ok = (
        len(importance) == len(targets) * len(targets) * len(expected_views)
        and set(importance["Target"]) == targets
        and set(importance["Predictor"]) == targets
        and set(importance["view"]) == expected_views
        and not importance.duplicated(["view", "Predictor", "Target"]).any()
    )
    # mistyR reports the self-predictor diagonal as structural NA in every view.
    # All off-diagonal importances must remain finite; no other NA is accepted.
    importance_ok = (
        importance_grid_ok
        and missing_importance.equals(structural_self)
        and np.isfinite(importance.loc[~structural_self, "Importance"].to_numpy(float)).all()
    )
    record(
        f"{sample}::misty_importance",
        importance_ok,
        f"rows={len(importance)};structural_self_na={int((missing_importance & structural_self).sum())};"
        f"unexpected_na={int((missing_importance & ~structural_self).sum())}",
    )
    for stem in ("misty_performance_gains_v1", "misty_view_contributions_v1", "misty_juxta_para_importance_v1"):
        verify_image(directory / f"{sample}_{stem}.png")


def verify_commot(sample: str, expected_spots: int) -> None:
    directory = ROOT / "commot" / "samples" / sample
    status = read_json(directory / "status.json")
    meta = read_json(directory / "run_metadata_v2.json")
    record(f"{sample}::commot_status", status.get("status") == "complete" and meta.get("status") == "complete" and int(meta.get("n_spots", -1)) == expected_spots, json.dumps(status))
    manifest = read_tsv(directory / f"{sample}_commot_candidate_manifest_v2.tsv")
    summary = read_tsv(directory / f"{sample}_commot_transport_summary_v2.tsv")
    edge = read_tsv(directory / f"{sample}_commot_celltype_edge_mass_proxy_v2.tsv")
    hypothesis = read_tsv(directory / f"{sample}_commot_hypothesis_proxy_support_v2.tsv")
    domain = read_tsv(directory / f"{sample}_commot_computational_domain_edge_mass_v2.tsv")
    spot = read_tsv(directory / f"{sample}_commot_spot_sender_receiver_v2.tsv")
    available = bools(manifest["available"])
    n_axes = int(available.sum()); expected_rows = int(meta["n_secreted_axes"]) * 2 + int(meta["n_contact_axes"])
    record(f"{sample}::commot_estimability", len(manifest) == 12 and manifest["axis_id"].is_unique and n_axes == int(meta["n_candidate_axes_modeled"]) and n_axes >= 8, f"available={n_axes}/12;unavailable={';'.join(manifest.loc[~available, 'axis_id'])}")
    finite = np.isfinite(summary[["total_score", "n_edges", "weighted_p90_distance_um"]].to_numpy(float)).all()
    summary_ok = len(summary) == expected_rows and not summary.duplicated(["scale", "axis_id"]).any() and finite and (summary["total_score"] > 0).all()
    record(f"{sample}::commot_summary", summary_ok, f"rows={len(summary)};expected={expected_rows}")
    edge_groups = edge.groupby(["scale", "axis_id"], observed=True)
    edge_ok = len(edge) == expected_rows * 196 and (edge_groups.size() == 196).all() and np.allclose(edge_groups["weighted_mass_fraction"].sum(), 1, atol=1e-9)
    record(f"{sample}::commot_edges", edge_ok, f"rows={len(edge)};groups={len(edge_groups)}")
    domain_groups = domain.groupby(["scale", "axis_id"], observed=True)
    domain_ok = len(domain) == expected_rows * 4 and (domain_groups.size() == 4).all() and np.allclose(domain_groups["edge_mass_fraction"].sum(), 1, atol=1e-9)
    record(f"{sample}::commot_domains", domain_ok, f"rows={len(domain)};groups={len(domain_groups)}")
    spot_ok = len(spot) == expected_spots and spot["barcode"].is_unique and len(hypothesis) == expected_rows
    max_error = 0.0
    for row in summary.itertuples():
        for prefix in ("sender", "receiver"):
            column = f"{prefix}__{row.scale}__{row.axis_id}"
            if column not in spot:
                spot_ok = False; continue
            error = abs(float(spot[column].sum()) - float(row.total_score)) / max(abs(float(row.total_score)), 1e-12)
            max_error = max(max_error, error); spot_ok = spot_ok and error < 1e-8 and (spot[column] >= 0).all()
    record(f"{sample}::commot_transport_balance", spot_ok, f"spot_rows={len(spot)};max_relative_error={max_error:.3g}")
    for stem in ("commot_total_score_v2", "commot_hypothesis_proxy_v2", "commot_selected_spatial_maps_v2"):
        verify_image(directory / f"{sample}_{stem}.png")


def verify_cohort() -> None:
    directory = ROOT / "cohort_synthesis_v1"
    status = read_json(directory / "status.json")
    meta = read_json(directory / "run_metadata_v1.json")
    record("cohort::status", status.get("status") == "completed" and meta.get("status") == "completed" and int(meta.get("n_patients", -1)) == 12, json.dumps(status))
    required = {
        "bayesspace_patient_summary_v1.tsv": 12,
        "bayesspace_patient_program_effects_v1.tsv": 84,
        "candidate_program_gene_replication_v1.tsv": sum(len(v) for v in {
            "a": ("CA9", "EPCAM", "KRT8", "KRT18", "KRT19"), "b": ("POSTN", "COL1A1", "COL1A2", "COL3A1", "COL4A1", "FN1", "SPARC"),
            "c": ("APOE", "C1QA", "C1QB", "C1QC", "TREM2", "LPL"), "d": ("CD8A", "CD8B", "NKG7", "CCL5", "PRF1", "GZMB"),
            "e": ("TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SERPINE1", "CTGF"), "f": ("TYMP", "SPP1", "MIF", "CD74", "CXCR4", "APP")}.values()),
        "commot_axis_estimability_v1.tsv": 12,
    }
    for name, minimum in required.items():
        frame = read_tsv(directory / name)
        record(f"cohort::shape::{name}", len(frame) >= minimum, f"rows={len(frame)};minimum={minimum}")
    for name in (
        "spatial_program_gene_replication_v1.tsv",
        "misty_patient_performance_summary_v1.tsv", "misty_patient_view_contribution_summary_v1.tsv",
        "misty_patient_importance_replication_v1.tsv", "commot_patient_transport_summary_v1.tsv",
        "commot_patient_hypothesis_proxy_summary_v1.tsv", "commot_patient_celltype_proxy_replication_v1.tsv",
        "commot_patient_domain_transport_v1.tsv",
    ):
        frame = read_tsv(directory / name)
        record(f"cohort::nonempty::{name}", not frame.empty, f"rows={len(frame)}")
    module_cohort = read_optional_tsv(directory / "hotspot_patient_module_summary_v1.tsv")
    record("cohort::hotspot_module_table", module_cohort.empty or "sample" in module_cohort, f"rows={len(module_cohort)}")
    for stem in (
        "bayesspace_patient_overview_v1", "bayesspace_program_effects_v1", "spatial_program_candidate_replication_v1",
        "spatial_program_top_recurrent_genes_v1", "misty_patient_performance_v1", "misty_patient_view_contributions_v1",
        "commot_patient_transport_rank_v1", "commot_patient_hypothesis_proxy_v1",
    ):
        verify_image(directory / f"{stem}.png")


def build_manifest(output: Path) -> int:
    rows = []
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or output in path.parents or "/logs/" in str(path):
            continue
        suffix = path.suffix.lower()
        hash_file = suffix not in {".rds", ".h5ad"}
        rows.append({
            "relative_path": str(path.relative_to(ROOT)), "size_bytes": path.stat().st_size,
            "sha256": sha256(path) if hash_file else "",
            "hash_policy": "sha256" if hash_file else "size_only_for_bulk_binary",
        })
    pd.DataFrame(rows).to_csv(output / "file_manifest_v1.tsv", sep="\t", index=False)
    return len(rows)


def main(output: Path) -> int:
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite: {output}")
    output.mkdir(parents=True)
    (output / "status.json").write_text('{"status":"running"}\n')
    gate("tjrcc_advanced_independent_qc_start")
    total_spots = 0
    spot_counts: dict[str, int] = {}
    for sample in SAMPLES:
        gate(f"tjrcc_qc_{sample.lower()}_bayesspace")
        spot_counts[sample] = verify_bayesspace(sample); total_spots += spot_counts[sample]
    gate("tjrcc_advanced_independent_qc_spatial_programs")
    for sample in SAMPLES:
        verify_spatial_programs(sample, spot_counts[sample])
    gate("tjrcc_advanced_independent_qc_misty")
    for sample in SAMPLES:
        verify_misty(sample, spot_counts[sample])
    gate("tjrcc_advanced_independent_qc_commot")
    for sample in SAMPLES:
        verify_commot(sample, spot_counts[sample])
    gate("tjrcc_advanced_independent_qc_cohort")
    verify_cohort()
    failures = sum(not bool(row["passed"]) for row in checks)
    pd.DataFrame(checks).to_csv(output / "verification_checks_v1.tsv", sep="\t", index=False)
    manifest_n = build_manifest(output)
    payload = {
        "status": "passed" if failures == 0 else "failed", "n_checks": len(checks), "n_failed": failures,
        "n_samples": len(SAMPLES), "n_spots": total_spots, "file_manifest_rows": manifest_n,
        "memory_gate_records": memory_records,
        "inference_boundary": "QC verifies computation and internal consistency; it does not validate computational domains against pathology or convert spot-level tests into patient-level causal evidence.",
    }
    (output / "verification_summary_v1.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    (output / "status.json").write_text(json.dumps({"status": payload["status"]}, indent=2) + "\n")
    if failures:
        raise RuntimeError(f"Independent QC failed {failures}/{len(checks)} checks")
    print(f"TJRCC_ADVANCED_INDEPENDENT_QC_ALL_PASS checks={len(checks)} files={manifest_n}", flush=True)
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: verify_tjrcc_advanced_full_v1.py OUTPUT_DIR")
    target = Path(sys.argv[1]).resolve()
    existed = target.exists()
    try:
        raise SystemExit(main(target))
    except Exception as exc:
        if not existed:
            target.mkdir(parents=True, exist_ok=True)
            prior = {}
            if (target / "status.json").exists():
                try:
                    prior = json.loads((target / "status.json").read_text())
                except Exception:
                    prior = {}
            if prior.get("status") != "failed":
                (target / "status.json").write_text(json.dumps({"status": "failed", "error": str(exc)}, indent=2) + "\n")
            (target / "failure_traceback.txt").write_text(traceback.format_exc(), encoding="utf-8")
        raise
