#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.mkdir(parents=True)
    required = [
        "tcga_kirc_natural_perturbation_cohort_v1.tsv",
        "tcga_kirc_interaction_primary_v1.tsv",
        "tcga_kirc_interaction_sensitivity_v1.tsv",
        "external_scrna_interaction_primary_v1.tsv",
        "external_scrna_leave_one_cohort_out_v1.tsv",
        "mtap_molecular_validation_gate_v1.json",
        "PHASE79_DECISION_v1.json",
        "PHASE79_PASS_v1",
    ]
    missing = [name for name in required if not (args.results / name).is_file()]
    if missing:
        raise RuntimeError(f"Missing outputs: {missing}")
    cohort = pd.read_csv(
        args.results / "tcga_kirc_natural_perturbation_cohort_v1.tsv", sep="\t"
    )
    primary = pd.read_csv(
        args.results / "tcga_kirc_interaction_primary_v1.tsv", sep="\t"
    )
    external = pd.read_csv(
        args.results / "external_scrna_interaction_primary_v1.tsv", sep="\t"
    )
    loo = pd.read_csv(
        args.results / "external_scrna_leave_one_cohort_out_v1.tsv", sep="\t"
    )
    gate = json.loads(
        (args.results / "mtap_molecular_validation_gate_v1.json").read_text()
    )
    decision = json.loads(
        (args.results / "PHASE79_DECISION_v1.json").read_text()
    )
    target = primary.loc[primary["outcome"] == "spp1_state"].iloc[0]
    external_target = external.loc[external["outcome"] == "spp1_state"].iloc[0]
    loo_target = loo.loc[loo["outcome"] == "spp1_state"]
    genomic = bool(gate["pass"] and target["beta"] > 0 and target["q"] < 0.05)
    external_support = bool(
        external_target["beta"] > 0
        and external_target["q"] < 0.05
        and (loo_target["beta"] > 0).all()
    )
    checks = pd.DataFrame([
        {"check": "tcga_unique_patients", "value": cohort["patient"].nunique(), "pass": cohort["patient"].nunique() == len(cohort)},
        {"check": "tcga_at_least_500", "value": len(cohort), "pass": len(cohort) >= 500},
        {"check": "primary_three_plus_control", "value": len(primary), "pass": len(primary) == 4},
        {"check": "external_three_endpoints", "value": len(external), "pass": len(external) == 3},
        {"check": "loo_two_cohorts_each", "value": len(loo), "pass": len(loo) == 6},
        {"check": "molecular_gate_true", "value": gate["pass"], "pass": gate["pass"] is True},
        {"check": "decision_genomic_recomputed", "value": decision["genomic_support"], "pass": decision["genomic_support"] == genomic},
        {"check": "decision_external_recomputed", "value": decision["external_scrna_support"], "pass": decision["external_scrna_support"] == external_support},
        {"check": "primary_q_valid", "value": int(np.isfinite(primary.loc[primary["outcome"] != "c1qc_specificity_control", "q"]).sum()), "pass": primary.loc[primary["outcome"] != "c1qc_specificity_control", "q"].between(0, 1).all()},
    ])
    checks.to_csv(args.output / "phase79_validation_checks_v1.tsv", sep="\t", index=False)
    manifest = pd.DataFrame([
        {"file": name, "sha256": sha256(args.results / name)}
        for name in required
    ])
    manifest.to_csv(args.output / "phase79_output_manifest_sha256_v1.tsv", sep="\t", index=False)
    if not checks["pass"].all():
        raise RuntimeError("Phase79 validation failed")
    (args.output / "PHASE79_VALIDATION_PASS_v1").write_text(
        "status=PASS\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
