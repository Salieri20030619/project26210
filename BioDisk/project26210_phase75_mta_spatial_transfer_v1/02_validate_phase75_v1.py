#!/usr/bin/env python3
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


def bh(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(len(values))
    out[order] = np.minimum(ranked, 1.0)
    return out


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def check(condition, name, observed, expected, rows):
    rows.append({"check": name, "status": "PASS" if condition else "FAIL", "observed": observed, "expected": expected})


def replay_meta(meta, family_col):
    max_error = 0.0
    decision_mismatch = 0
    for _, idx in meta.groupby(["preparation", family_col], sort=True).groups.items():
        idx = list(idx)
        expected_q = bh(meta.loc[idx, "wilcoxon_signed_rank_p"].to_numpy(dtype=float))
        max_error = max(max_error, float(np.max(np.abs(expected_q - meta.loc[idx, "wilcoxon_bh_q"].to_numpy(dtype=float)))))
    expected_status = np.where(
        (meta["wilcoxon_bh_q"] < 0.05)
        & (meta["median_rho"] * meta["expected_direction"] > 0)
        & (meta["expected_direction_fraction"] >= 8 / 12),
        "SUPPORTED", "NOT_SUPPORTED"
    )
    decision_mismatch = int(np.sum(expected_status != meta["support_status"].to_numpy()))
    return max_error, decision_mismatch


def replay_cross(meta, decisions, id_col, family_col, family_value):
    mismatch = 0
    sub = meta[meta[family_col] == family_value]
    for _, row in decisions.iterrows():
        x = sub[sub[id_col] == row[id_col]]
        status = dict(zip(x["preparation"], x["support_status"]))
        n = sum(status.get(p) == "SUPPORTED" for p in ["FFPE", "FROZEN"])
        expected = "TWO_PREPARATION_SUPPORTED" if n == 2 else ("ONE_PREPARATION_SUPPORTED" if n == 1 else "NOT_SUPPORTED")
        mismatch += int(row["ffpe_status"] != status.get("FFPE", "MISSING"))
        mismatch += int(row["frozen_status"] != status.get("FROZEN", "MISSING"))
        mismatch += int(int(row["preparations_supported"]) != n)
        mismatch += int(row["decision"] != expected)
    return mismatch


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    results = Path(args.results).resolve()
    out = Path(args.out).resolve()
    part = Path(str(out) + ".part")
    if out.exists() or part.exists():
        raise RuntimeError("Refusing existing validation output")
    part.mkdir(parents=True)

    section = pd.read_csv(results / "section_qc_v1.tsv", sep="\t")
    coverage = pd.read_csv(results / "program_gene_coverage_v1.tsv", sep="\t")
    relations = pd.read_csv(results / "relations_by_section_v1.tsv", sep="\t")
    rmeta = pd.read_csv(results / "relation_meta_v1.tsv", sep="\t")
    pdec = pd.read_csv(results / "primary_relation_decisions_v1.tsv", sep="\t")
    ndec = pd.read_csv(results / "neighbor_relation_decisions_v1.tsv", sep="\t")
    chains = pd.read_csv(results / "chain_products_by_section_v1.tsv", sep="\t")
    cmeta = pd.read_csv(results / "chain_product_meta_v1.tsv", sep="\t")
    cdec = pd.read_csv(results / "chain_product_decisions_v1.tsv", sep="\t")

    rows = []
    check(len(section) == 24, "sections_24", len(section), 24, rows)
    check((section["preparation"] == "FFPE").sum() == 12, "ffpe_12", int((section["preparation"] == "FFPE").sum()), 12, rows)
    check((section["preparation"] == "FROZEN").sum() == 12, "frozen_12", int((section["preparation"] == "FROZEN").sum()), 12, rows)
    check((section["exact_barcode_match"] == "PASS").all(), "exact_barcode_all", int((section["exact_barcode_match"] == "PASS").sum()), 24, rows)
    check((section["variable_response_genes"] >= 2).all(), "response_variance_gate", int((section["variable_response_genes"] >= 2).sum()), 24, rows)
    check(len(coverage) == 216, "coverage_rows", len(coverage), 216, rows)
    check(len(relations) == 288, "relation_rows", len(relations), 288, rows)
    check(np.isfinite(relations["rho"]).all(), "finite_relation_rho", int(np.isfinite(relations["rho"]).sum()), 288, rows)
    check(len(rmeta) == 24, "relation_meta_rows", len(rmeta), 24, rows)
    check(len(pdec) == 6 and len(ndec) == 6, "relation_decision_rows", f"{len(pdec)}+{len(ndec)}", "6+6", rows)
    check(len(chains) == 96 and np.isfinite(chains["rho"]).all(), "chain_rows_finite", len(chains), 96, rows)
    check(len(cmeta) == 8 and len(cdec) == 4, "chain_meta_decisions", f"{len(cmeta)}+{len(cdec)}", "8+4", rows)

    rerr, rmismatch = replay_meta(rmeta, "analysis_family")
    cerr, cmismatch = replay_meta(cmeta, "analysis_family")
    check(rerr < 1e-12, "relation_bh_replay", rerr, "<1e-12", rows)
    check(cerr < 1e-12, "chain_bh_replay", cerr, "<1e-12", rows)
    check(rmismatch == 0, "relation_support_replay", rmismatch, 0, rows)
    check(cmismatch == 0, "chain_support_replay", cmismatch, 0, rows)
    cross_mismatch = replay_cross(rmeta, pdec, "relation_id", "analysis_family", "PRIMARY_SAME_SPOT")
    cross_mismatch += replay_cross(rmeta, ndec, "relation_id", "analysis_family", "SUPPORTING_K6_NEIGHBOR")
    cross_mismatch += replay_cross(cmeta, cdec, "chain_id", "analysis_family", "PRIMARY_CHAIN_PRODUCT")
    check(cross_mismatch == 0, "cross_preparation_decision_replay", cross_mismatch, 0, rows)

    checks = pd.DataFrame(rows)
    checks.to_csv(part / "validation_checks_v1.tsv", sep="\t", index=False)
    status = "PASS" if (checks["status"] == "PASS").all() else "FAIL"
    summary = {"status": status, "checks": len(checks), "passed": int((checks["status"] == "PASS").sum()), "failed": int((checks["status"] == "FAIL").sum())}
    (part / "validation_summary_v1.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    core = sorted(p for p in results.iterdir() if p.is_file())
    pd.DataFrame([{"file": p.name, "bytes": p.stat().st_size, "sha256": sha256(p)} for p in core]).to_csv(part / "CORE_OUTPUT_SHA256_v1.tsv", sep="\t", index=False)
    if status != "PASS":
        raise RuntimeError(f"Phase75 independent validation failed: {summary}")
    (part / "INDEPENDENT_VALIDATION_PASS_v1").write_text("status=PASS\n", encoding="utf-8")
    os.replace(part, out)
    print(json.dumps(summary))


if __name__ == "__main__":
    main()

