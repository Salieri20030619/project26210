#!/usr/bin/env python3
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bh_adjust(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    results = Path(args.results).resolve()
    out = Path(args.out).resolve()
    if out.exists():
        raise RuntimeError(f"Refusing to overwrite {out}")
    out.mkdir(parents=True)

    checks = []

    def check(name, passed, observed, expected):
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "observed": observed, "expected": expected})

    required = [
        "PHASE70_PASS_v1",
        "PHASE70_RUN_STATUS_v1.tsv",
        "section_join_qc_v1.tsv",
        "construct_correlations_by_section_v1.tsv",
        "composition_adjusted_bridges_by_section_v1.tsv",
        "construct_cross_section_meta_v1.tsv",
        "composition_adjusted_bridge_meta_v1.tsv",
        "primary_decisions_v1.tsv",
        "INPUT_MANIFEST_SHA256_v1.tsv",
        "OUTPUT_MANIFEST_SHA256_v1.tsv",
        "PHASE70_COMPOSITION_AWARE_SPATIAL_SUMMARY_v1.md",
        "run_metadata_v1.json",
    ]
    missing = [name for name in required if not (results / name).is_file()]
    check("required_files", not missing, ";".join(missing) if missing else "all present", "all present")
    if missing:
        raise RuntimeError(f"Missing required files: {missing}")

    status = pd.read_csv(results / "PHASE70_RUN_STATUS_v1.tsv", sep="\t")
    qc = pd.read_csv(results / "section_join_qc_v1.tsv", sep="\t")
    constructs = pd.read_csv(results / "construct_correlations_by_section_v1.tsv", sep="\t")
    bridges = pd.read_csv(results / "composition_adjusted_bridges_by_section_v1.tsv", sep="\t")
    construct_meta = pd.read_csv(results / "construct_cross_section_meta_v1.tsv", sep="\t")
    bridge_meta = pd.read_csv(results / "composition_adjusted_bridge_meta_v1.tsv", sep="\t")
    decisions = pd.read_csv(results / "primary_decisions_v1.tsv", sep="\t")
    inputs = pd.read_csv(results / "INPUT_MANIFEST_SHA256_v1.tsv", sep="\t")
    outputs = pd.read_csv(results / "OUTPUT_MANIFEST_SHA256_v1.tsv", sep="\t")
    metadata = json.loads((results / "run_metadata_v1.json").read_text(encoding="utf-8"))

    check("run_gates", bool((status["status"] == "PASS").all()), int((status["status"] == "PASS").sum()), len(status))
    check("section_count", len(qc) == 24, len(qc), 24)
    check("preparation_balance", qc["preparation"].value_counts().to_dict() == {"FFPE": 12, "FROZEN": 12}, str(qc["preparation"].value_counts().to_dict()), "FFPE=12;FROZEN=12")
    check("exact_barcode_match", bool((qc["exact_barcode_match"] == "PASS").all()), int((qc["exact_barcode_match"] == "PASS").sum()), 24)
    check("construct_rows", len(constructs) == 576, len(constructs), 576)
    check("bridge_rows", len(bridges) == 216, len(bridges), 216)
    check("construct_meta_rows", len(construct_meta) == 48, len(construct_meta), 48)
    check("bridge_meta_rows", len(bridge_meta) == 18, len(bridge_meta), 18)
    check("decision_rows", len(decisions) == 9, len(decisions), 9)
    check("finite_primary_coefficients", bool(np.isfinite(bridges["rho"]).all()), int(np.isfinite(bridges["rho"]).sum()), 216)
    check("section_level_meta_n", bool((bridge_meta["n_sections"] == 12).all()), sorted(bridge_meta["n_sections"].unique().tolist()), [12])

    bh_max_error = 0.0
    for meta in [construct_meta, bridge_meta]:
        for _, sub in meta.groupby(["preparation", "analysis_family"], sort=True):
            expected = bh_adjust(sub["wilcoxon_signed_rank_p"].to_numpy(dtype=float))
            bh_max_error = max(bh_max_error, float(np.max(np.abs(expected - sub["wilcoxon_bh_q"].to_numpy(dtype=float)))))
    check("bh_recalculation", bh_max_error <= 1e-12, bh_max_error, "<=1e-12")

    recomputed = []
    for relation_id, sub in bridge_meta.groupby("relation_id", sort=True):
        supported = int((sub["support_status"] == "SUPPORTED").sum())
        recomputed.append(
            "TWO_PREPARATION_SUPPORTED" if supported == 2 else ("ONE_PREPARATION_SUPPORTED" if supported == 1 else "NOT_SUPPORTED")
        )
    observed_decisions = decisions.sort_values("relation_id")["decision"].tolist()
    check("decision_recalculation", observed_decisions == recomputed, str(observed_decisions), str(recomputed))

    bad_inputs = []
    for row in inputs.itertuples(index=False):
        path = Path(row.path)
        if not path.is_file() or path.stat().st_size != int(row.bytes) or sha256(path) != row.sha256:
            bad_inputs.append(str(path))
    check("input_hash_replay", not bad_inputs, len(inputs) - len(bad_inputs), len(inputs))

    bad_outputs = []
    for row in outputs.itertuples(index=False):
        path = results / row.file
        if not path.is_file() or path.stat().st_size != int(row.bytes) or sha256(path) != row.sha256:
            bad_outputs.append(str(path))
    check("output_hash_replay", not bad_outputs, len(outputs) - len(bad_outputs), len(outputs))
    check("metadata_consistency", metadata.get("sections") == 24 and metadata.get("spots") == int(qc["spots"].sum()), f"sections={metadata.get('sections')};spots={metadata.get('spots')}", f"sections=24;spots={int(qc['spots'].sum())}")
    summary = (results / "PHASE70_COMPOSITION_AWARE_SPATIAL_SUMMARY_v1.md").read_text(encoding="utf-8")
    check("interpretation_boundary", "does not alter the frozen response/PFS evidence grade" in summary, "present" if "does not alter the frozen response/PFS evidence grade" in summary else "missing", "present")

    frame = pd.DataFrame(checks)
    frame.to_csv(out / "validation_checks_v1.tsv", sep="\t", index=False)
    passed = int((frame["status"] == "PASS").sum())
    status_json = {"status": "PASS" if passed == len(frame) else "FAIL", "checks_passed": passed, "checks_total": len(frame)}
    (out / "validation_summary_v1.json").write_text(json.dumps(status_json, indent=2) + "\n", encoding="utf-8")
    if passed != len(frame):
        raise RuntimeError(frame.loc[frame["status"] != "PASS"].to_string(index=False))
    (out / "PHASE70_INDEPENDENT_VALIDATION_PASS_v1").write_text("status=PASS\n", encoding="utf-8")


if __name__ == "__main__":
    main()
