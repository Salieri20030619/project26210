#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from scipy.stats import spearmanr


ROOT = Path("/mnt/BioDisk/ludingyu/project26210/results/part4_spatial")
CROSS = ROOT / "advanced_patch_commot_cross_sample_v3"
SAMPLES = ("GLMF1", "GLMF2")
checks: list[dict[str, object]] = []


def record(name: str, passed: bool, detail: str) -> None:
    checks.append({"check": name, "passed": bool(passed), "detail": detail})


def read_tsv(path: Path) -> pd.DataFrame:
    record(f"exists::{path.name}", path.exists() and path.stat().st_size > 0, str(path))
    return pd.read_csv(path, sep="\t")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_image(path: Path) -> None:
    with Image.open(path) as image:
        array = np.asarray(image.convert("RGB").resize((128, 128)))
        passed = image.width >= 1000 and image.height >= 1000 and float(array.std()) > 5
        record(f"image::{path.name}", passed, f"{image.width}x{image.height}; rgb_sd={array.std():.3f}")


def verify_sample(sample: str) -> None:
    directory = ROOT / f"advanced_patch_commot_{sample.lower()}_v2"
    status = json.loads((directory / "status.json").read_text())
    metadata = json.loads((directory / "run_metadata_v2.json").read_text())
    record(f"{sample}::status", status == {"status": "complete", "sample": sample}, json.dumps(status))
    record(f"{sample}::metadata", metadata.get("status") == "complete" and metadata.get("commot_version") == "0.0.3",
           f"status={metadata.get('status')}; commot={metadata.get('commot_version')}")

    manifest = read_tsv(directory / f"{sample}_commot_candidate_manifest_v2.tsv")
    summary = read_tsv(directory / f"{sample}_commot_transport_summary_v2.tsv")
    edge = read_tsv(directory / f"{sample}_commot_celltype_edge_mass_proxy_v2.tsv")
    hypothesis = read_tsv(directory / f"{sample}_commot_hypothesis_proxy_support_v2.tsv")
    domain = read_tsv(directory / f"{sample}_commot_computational_domain_edge_mass_v2.tsv")
    spot = read_tsv(directory / f"{sample}_commot_spot_sender_receiver_v2.tsv")
    record(f"{sample}::manifest", len(manifest) == 12 and manifest["axis_id"].nunique() == 12 and manifest["available"].all(),
           f"rows={len(manifest)}; available={int(manifest['available'].sum())}")
    scale_counts = summary.groupby("scale").size().to_dict()
    expected_scales = {"secreted_200um": 11, "secreted_500um": 11, "contact_110um": 1}
    record(f"{sample}::summary_shape", len(summary) == 23 and scale_counts == expected_scales and not summary.duplicated(["scale", "axis_id"]).any(),
           f"rows={len(summary)}; scales={scale_counts}")
    numeric_ok = np.isfinite(summary[["total_score", "n_edges", "diagonal_share", "weighted_p90_distance_um"]]).all(axis=None)
    numeric_ok = numeric_ok and (summary["total_score"] > 0).all() and (summary["n_edges"] > 0).all()
    record(f"{sample}::summary_numeric", numeric_ok, f"min_total={summary.total_score.min():.6g}")
    threshold = summary["scale"].map({"secreted_200um": 200.0, "secreted_500um": 500.0, "contact_110um": 110.0})
    record(f"{sample}::distance_threshold", (summary["weighted_p90_distance_um"] <= threshold + 1e-6).all(),
           f"max_excess={(summary['weighted_p90_distance_um'] - threshold).max():.6g}")

    expected_spots = int(metadata["n_spots"])
    sender_columns = [column for column in spot if column.startswith("sender__")]
    receiver_columns = [column for column in spot if column.startswith("receiver__")]
    record(f"{sample}::spot_shape", len(spot) == expected_spots and len(sender_columns) == 23 and len(receiver_columns) == 23,
           f"shape={spot.shape}; sender={len(sender_columns)}; receiver={len(receiver_columns)}")
    balance_ok = True
    max_relative_error = 0.0
    for row in summary.itertuples():
        sender = f"sender__{row.scale}__{row.axis_id}"
        receiver = f"receiver__{row.scale}__{row.axis_id}"
        for column in (sender, receiver):
            observed = float(spot[column].sum())
            relative_error = abs(observed - row.total_score) / max(abs(row.total_score), 1e-12)
            max_relative_error = max(max_relative_error, relative_error)
            balance_ok = balance_ok and relative_error < 1e-9 and (spot[column] >= 0).all()
    record(f"{sample}::spot_transport_balance", balance_ok, f"max_relative_error={max_relative_error:.3g}")

    group = edge.groupby(["scale", "axis_id"], observed=True)
    sizes = group.size()
    mass_sums = group["weighted_mass_fraction"].sum()
    expected_sums = group["marginal_independence_expected_fraction"].sum()
    edge_ok = len(edge) == 4508 and (sizes == 196).all() and np.allclose(mass_sums, 1, atol=1e-10) and np.allclose(expected_sums, 1, atol=1e-10)
    record(f"{sample}::celltype_edge_mass", edge_ok,
           f"rows={len(edge)}; groups={len(sizes)}; mass_sum_range={mass_sums.min():.12g}-{mass_sums.max():.12g}")
    matched = hypothesis.merge(
        edge,
        left_on=["scale", "axis_id", "source_proxy", "target_proxy"],
        right_on=["scale", "axis_id", "source_celltype", "target_celltype"],
        suffixes=("_hyp", "_edge"), validate="one_to_one",
    )
    hypothesis_ok = len(hypothesis) == 23 and len(matched) == 23
    for left, right in [("proxy_log2_enrichment_hyp", "proxy_log2_enrichment_edge"),
                        ("proxy_mass_rank_among_196", "mass_rank_among_196")]:
        hypothesis_ok = hypothesis_ok and np.allclose(matched[left], matched[right], equal_nan=True)
    record(f"{sample}::hypothesis_matches_edge", hypothesis_ok, f"hypotheses={len(hypothesis)}; matched={len(matched)}")

    domain_group = domain.groupby(["scale", "axis_id"], observed=True)
    domain_sums = domain_group["edge_mass_fraction"].sum()
    domain_ok = len(domain) == 92 and (domain_group.size() == 4).all() and np.allclose(domain_sums, 1, atol=1e-10)
    record(f"{sample}::domain_mass", domain_ok, f"rows={len(domain)}; groups={len(domain_sums)}")
    h5ad = directory / f"{sample}_commot_candidate_v2.h5ad"
    record(f"{sample}::h5ad", h5ad.exists() and h5ad.stat().st_size > 1_000_000,
           f"bytes={h5ad.stat().st_size if h5ad.exists() else 0}")
    for stem in ["hypothesis_proxy", "selected_spatial_maps"]:
        verify_image(directory / f"{sample}_commot_{stem}_v2.png")
    verify_image(ROOT / "advanced_patch_commot_figures_v3" / f"{sample}_commot_total_score_v3.png")


def verify_offdiagonal(sample: str) -> None:
    directory = ROOT / "advanced_patch_commot_offdiagonal_v3" / sample
    status = json.loads((directory / "status.json").read_text())
    metadata = json.loads((directory / "run_metadata_v3.json").read_text())
    record(f"{sample}::offdiagonal_status", status.get("status") == "complete" and metadata.get("status") == "complete",
           f"status={status.get('status')}; rows={metadata.get('n_axes_scales')}")
    transport = read_tsv(directory / f"{sample}_commot_offdiagonal_transport_v3.tsv")
    edge = read_tsv(directory / f"{sample}_commot_offdiagonal_celltype_edge_mass_proxy_v3.tsv")
    hypothesis = read_tsv(directory / f"{sample}_commot_offdiagonal_hypothesis_proxy_support_v3.tsv")
    spot = read_tsv(directory / f"{sample}_commot_offdiagonal_spot_sender_receiver_v3.tsv")
    source_summary = pd.read_csv(ROOT / f"advanced_patch_commot_{sample.lower()}_v2" / f"{sample}_commot_transport_summary_v2.tsv", sep="\t")
    joined = transport.merge(source_summary[["scale", "axis_id", "offdiagonal_score", "diagonal_score"]],
                             on=["scale", "axis_id"], suffixes=("_v3", "_v2"), validate="one_to_one")
    transport_ok = len(transport) == 23 and len(joined) == 23
    transport_ok = transport_ok and np.allclose(joined["offdiagonal_score_v3"], joined["offdiagonal_score_v2"], rtol=1e-10, atol=1e-10)
    transport_ok = transport_ok and np.allclose(joined["diagonal_score_v3"], joined["diagonal_score_v2"], rtol=1e-10, atol=1e-10)
    record(f"{sample}::offdiagonal_transport", transport_ok,
           f"rows={len(transport)}; diagonal_share={transport.diagonal_share.min():.3f}-{transport.diagonal_share.max():.3f}")
    groups = edge.groupby(["scale", "axis_id"], observed=True)
    mass_sums = groups["offdiagonal_weighted_mass_fraction"].sum()
    expected_sums = groups["offdiagonal_expected_fraction"].sum()
    edge_ok = len(edge) == 4508 and (groups.size() == 196).all() and np.allclose(mass_sums, 1, atol=1e-10)
    edge_ok = edge_ok and np.allclose(expected_sums, 1, atol=1e-10)
    record(f"{sample}::offdiagonal_edge_mass", edge_ok, f"rows={len(edge)}; groups={len(groups)}")
    matched = hypothesis.merge(
        edge,
        left_on=["scale", "axis_id", "source_proxy", "target_proxy"],
        right_on=["scale", "axis_id", "source_celltype", "target_celltype"],
        suffixes=("_hyp", "_edge"), validate="one_to_one",
    )
    hypothesis_ok = len(hypothesis) == 23 and len(matched) == 23
    hypothesis_ok = hypothesis_ok and np.allclose(
        matched["offdiagonal_proxy_log2_enrichment_hyp"], matched["offdiagonal_proxy_log2_enrichment_edge"], equal_nan=True
    )
    hypothesis_ok = hypothesis_ok and np.array_equal(
        matched["offdiagonal_proxy_mass_rank_among_196"], matched["offdiagonal_mass_rank_among_196"]
    )
    record(f"{sample}::offdiagonal_hypothesis", hypothesis_ok, f"rows={len(hypothesis)}; matched={len(matched)}")
    sender_columns = [column for column in spot if column.startswith("sender__")]
    receiver_columns = [column for column in spot if column.startswith("receiver__")]
    balance_ok = len(sender_columns) == 23 and len(receiver_columns) == 23
    for row in transport.itertuples():
        balance_ok = balance_ok and np.isclose(spot[f"sender__{row.scale}__{row.axis_id}"].sum(), row.offdiagonal_score, rtol=1e-9)
        balance_ok = balance_ok and np.isclose(spot[f"receiver__{row.scale}__{row.axis_id}"].sum(), row.offdiagonal_score, rtol=1e-9)
    record(f"{sample}::offdiagonal_spot_balance", balance_ok, f"shape={spot.shape}")


def verify_failed_v1_and_reproducibility() -> None:
    v1 = ROOT / "advanced_patch_commot_glmf1_v1"
    status = json.loads((v1 / "status.json").read_text())
    audit = v1 / "failure_audit_v1.md"
    record("v1::failure_preserved", status.get("status") == "failed" and audit.exists(),
           f"status={status.get('status')}; audit_sha256={sha256(audit) if audit.exists() else 'missing'}")
    equal = True
    compared = []
    maximum_absolute_difference = 0.0
    maximum_relative_difference = 0.0
    for stem in ["candidate_manifest", "transport_summary", "celltype_edge_mass_proxy", "hypothesis_proxy_support",
                 "computational_domain_edge_mass", "spot_sender_receiver"]:
        left = pd.read_csv(v1 / f"GLMF1_commot_{stem}_v1.tsv", sep="\t")
        right = pd.read_csv(ROOT / "advanced_patch_commot_glmf1_v2" / f"GLMF1_commot_{stem}_v2.tsv", sep="\t")
        numeric = left.select_dtypes(include="number").columns.tolist()
        nonnumeric = [column for column in left.columns if column not in numeric]
        structure_equal = left.columns.tolist() == right.columns.tolist() and left[nonnumeric].equals(right[nonnumeric])
        left_values = left[numeric].to_numpy(dtype=float)
        right_values = right[numeric].to_numpy(dtype=float)
        finite = np.isfinite(left_values) & np.isfinite(right_values)
        differences = np.abs(left_values - right_values)
        relative = differences / np.maximum(np.abs(left_values), 1e-12)
        maximum_absolute_difference = max(maximum_absolute_difference, float(np.nanmax(differences[finite])))
        maximum_relative_difference = max(maximum_relative_difference, float(np.nanmax(relative[finite])))
        table_equal = structure_equal and np.allclose(left_values, right_values, rtol=1e-10, atol=1e-12, equal_nan=True)
        equal = equal and table_equal
        if table_equal:
            compared.append(stem)
    record("v1_v2::glmf1_numeric_reproducibility", equal,
           f"tolerance_equal_tables={len(compared)}/6; max_abs={maximum_absolute_difference:.3g}; "
           f"max_rel={maximum_relative_difference:.3g}; {','.join(compared)}")
    visual_audit = ROOT / "advanced_patch_commot_figures_v3" / "commot_total_score_v2_visual_failure_audit.md"
    render_metadata = ROOT / "advanced_patch_commot_figures_v3" / "render_metadata_v3.json"
    record("v2_v3::total_score_visual_supersession", visual_audit.exists() and render_metadata.exists(),
           f"audit={visual_audit.exists()}; metadata={render_metadata.exists()}")


def verify_cross() -> None:
    status = json.loads((CROSS / "status.json").read_text())
    metadata = json.loads((CROSS / "run_metadata_v3.json").read_text())
    record("cross::status", status.get("status") == "complete" and metadata.get("status") == "complete",
           f"status={status}; metadata_status={metadata.get('status')}")
    transport = read_tsv(CROSS / "commot_cross_sample_transport_summary_v3.tsv")
    concordance = read_tsv(CROSS / "commot_cross_sample_transport_rank_concordance_v3.tsv")
    hypothesis = read_tsv(CROSS / "commot_cross_sample_hypothesis_proxy_v3.tsv")
    primary = read_tsv(CROSS / "commot_cross_sample_primary_hypothesis_support_v3.tsv")
    edge = read_tsv(CROSS / "commot_cross_sample_celltype_edge_proxy_v3.tsv")
    supported = read_tsv(CROSS / "commot_cross_sample_replicated_top_decile_proxy_edges_v3.tsv")
    record("cross::shapes", len(transport) == 23 and len(concordance) == 2 and len(hypothesis) == 23
           and len(primary) == 12 and len(edge) == 4508,
           f"transport={len(transport)}; concordance={len(concordance)}; hypothesis={len(hypothesis)}; primary={len(primary)}; edge={len(edge)}")
    recomputed = ((edge["offdiagonal_mass_rank_among_196_GLMF1"] <= 20)
                  & (edge["offdiagonal_mass_rank_among_196_GLMF2"] <= 20)
                  & (edge["offdiagonal_proxy_log2_enrichment_GLMF1"] > 0)
                  & (edge["offdiagonal_proxy_log2_enrichment_GLMF2"] > 0))
    support_ok = (recomputed == edge["replicated_top_decile_positive_enrichment"]).all() and len(supported) == int(recomputed.sum())
    record("cross::edge_support_rule", support_ok, f"supported={int(recomputed.sum())}; output_rows={len(supported)}")
    primary_rule = ((primary["offdiagonal_proxy_mass_rank_among_196_GLMF1"] <= 20)
                    & (primary["offdiagonal_proxy_mass_rank_among_196_GLMF2"] <= 20)
                    & (primary["offdiagonal_proxy_log2_enrichment_GLMF1"] > 0)
                    & (primary["offdiagonal_proxy_log2_enrichment_GLMF2"] > 0))
    record("cross::primary_support_rule", (primary_rule == primary["replicated_proxy_support"]).all(),
           f"supported={int(primary_rule.sum())}")
    concordance_ok = True
    for row in concordance.itertuples():
        data = transport[transport["scale"] == row.scale]
        rho, pvalue = spearmanr(data["offdiagonal_score_rank_GLMF1"], data["offdiagonal_score_rank_GLMF2"])
        concordance_ok = concordance_ok and np.isclose(rho, row.spearman_rho) and np.isclose(pvalue, row.nominal_two_sided_p)
    record("cross::rank_concordance", concordance_ok, concordance.to_json(orient="records"))
    for stem in ["transport_rank", "hypothesis_proxy", "replicated_proxy_edges"]:
        verify_image(CROSS / f"commot_cross_sample_{stem}_v3.png")


def build_file_manifest() -> int:
    paths: set[Path] = set()
    for directory in [
        ROOT / "advanced_patch_commot_glmf1_v1",
        ROOT / "advanced_patch_commot_glmf1_v2",
        ROOT / "advanced_patch_commot_glmf2_v2",
        ROOT / "advanced_patch_commot_offdiagonal_v3",
        ROOT / "advanced_patch_commot_figures_v3",
        CROSS,
    ]:
        paths.update(path for path in directory.rglob("*") if path.is_file())
    for name in [
        "run_part4_commot_candidate_v1.py", "run_part4_commot_pipeline_v1.sh",
        "run_part4_commot_candidate_v2.py", "run_part4_commot_pipeline_v2.sh",
        "derive_part4_commot_offdiagonal_v3.py", "summarize_part4_commot_cross_sample_v3.py",
        "render_part4_commot_total_score_v3.py", "verify_part4_commot_v3.py",
        "commot_pipeline_v1.log", "commot_pipeline_v2.log",
    ]:
        path = ROOT / name
        if path.exists():
            paths.add(path)
    paths = {path for path in paths if not path.name.startswith("commot_v3_independent_qc")
             and path.name != "commot_v3_file_manifest.tsv"}
    rows = [{"absolute_path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in sorted(paths)]
    frame = pd.DataFrame(rows)
    frame.to_csv(CROSS / "commot_v3_file_manifest.tsv", sep="\t", index=False)
    record("bundle::file_manifest", len(frame) > 50 and (frame["bytes"] > 0).all() and frame["sha256"].str.len().eq(64).all(),
           f"files={len(frame)}")
    return len(frame)


def main() -> None:
    for sample in SAMPLES:
        verify_sample(sample)
        verify_offdiagonal(sample)
    verify_failed_v1_and_reproducibility()
    verify_cross()
    file_count = build_file_manifest()
    frame = pd.DataFrame(checks)
    frame.to_csv(CROSS / "commot_v3_independent_qc_checks.tsv", sep="\t", index=False)
    passed = bool(frame["passed"].all())
    result = {"status": "pass" if passed else "fail", "checks_passed": int(frame["passed"].sum()),
              "checks_total": len(frame), "files_hashed": file_count,
              "failed_checks": frame.loc[~frame["passed"], "check"].tolist()}
    (CROSS / "commot_v3_independent_qc_summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
