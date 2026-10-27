#!/usr/bin/env python3
"""Consolidate Phase 53 scoring-method rows into dataset-by-relation findings."""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "final_evidence_local_v1" / "final_bridge_evidence_v1.tsv"
FINAL = HERE / "consolidated_findings_local_v1"
PART = HERE / "consolidated_findings_local_v1.part"
REFERENCE_LAYER = "expanded_any_positive_cnv"


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    staging = path.with_name(path.name + ".part")
    frame.to_csv(staging, sep="\t", index=False, na_rep="NA", encoding="utf-8", lineterminator="\n")
    os.replace(staging, path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    if FINAL.exists() or PART.exists():
        raise FileExistsError("Refusing to overwrite an existing consolidated-findings directory")
    PART.mkdir(parents=True)
    source = pd.read_csv(SOURCE, sep="\t")
    source["reportable_v2"] = source["reportable_v2"].astype(str).str.lower() == "true"
    source["method_preference"] = source["method"].map({
        "primary_gene_z_mean": 0,
        "sensitivity_mean_log_normalized": 1,
    }).fillna(9)
    source["selection_p"] = pd.to_numeric(source["p_value_for_reporting"], errors="coerce").fillna(float("inf"))

    keys = ["dataset", "evidence_layer", "relation"]
    rows = []
    for _, group in source.groupby(keys, sort=True, dropna=False):
        reportable = group.loc[group["reportable_v2"]]
        candidates = reportable if not reportable.empty else group
        best = candidates.sort_values(["selection_p", "method_preference", "method"]).iloc[0]
        tested_rhos = pd.to_numeric(group["spearman_rho"], errors="coerce").dropna()
        signs = set((tested_rhos > 0).astype(int) - (tested_rhos < 0).astype(int))
        row = best.drop(labels=["method_preference", "selection_p"]).to_dict()
        row.update({
            "finding_reportable": not reportable.empty,
            "selected_method": best["method"],
            "methods_tested": ";".join(sorted(group["method"].astype(str).unique())),
            "n_methods_tested": int(group["method"].nunique()),
            "reportable_methods": ";".join(sorted(reportable["method"].astype(str).unique())) if not reportable.empty else "none",
            "n_methods_reportable": int(reportable["method"].nunique()),
            "method_direction_concordance": "concordant" if len(signs) <= 1 else "discordant",
            "finding_unit": "dataset_x_layer_x_relation; methods_are_robustness_checks_not_independent_replications",
        })
        rows.append(row)
    findings = pd.DataFrame(rows)
    reportable_findings = findings.loc[findings["finding_reportable"]].copy()
    independent_findings = reportable_findings.loc[
        reportable_findings["overlap_flag"] == "INDEPENDENT_GENE_SETS"
    ].copy()

    order = ["tumor_to_caf", "caf_to_tam", "tam_to_tcell", "tam_to_vascular", "tumor_to_tam", "tam_state_internal"]
    status_map = {
        "tumor_to_caf": "NO_REPORTABLE_SUPPORT",
        "caf_to_tam": "SMALL_N_INDEPENDENT_SUPPORT",
        "tam_to_tcell": "NO_REPORTABLE_SUPPORT",
        "tam_to_vascular": "NO_REPORTABLE_SUPPORT",
        "tumor_to_tam": "PROVISIONAL_LAYER_ONLY_SUPPORT",
        "tam_state_internal": "PARTIAL_INTERNAL_STATE_SUPPORT",
    }
    interpretations = {
        "tumor_to_caf": "F1/TYMP-to-POSTN CAF was not reportable; the main limitation is sparse patient overlap.",
        "caf_to_tam": "One independent POSTN CAF-to-APOE TAM finding was obtained in GSE264586 (n=5, exact p=0.0167).",
        "tam_to_tcell": "No TAM-to-T-cell dysfunction/cytotoxicity finding survived the final exact/FDR rule.",
        "tam_to_vascular": "No TAM-to-vascular finding survived; the candidate signatures also share LGALS3.",
        "tumor_to_tam": "Two findings occur only under the provisional tumour-like layer and disappear under expanded CNV-positive malignancy.",
        "tam_state_internal": "Two independent SPP1-to-suppression findings support an immunosuppressive TAM state; five APOE findings are state-consistency evidence because the signatures overlap.",
    }
    decision_rows = []
    for edge in order:
        rep = reportable_findings.loc[reportable_findings["edge_group"] == edge]
        indep = independent_findings.loc[independent_findings["edge_group"] == edge]
        caution = rep.loc[rep["overlap_flag"] != "INDEPENDENT_GENE_SETS"]
        decision_rows.append({
            "edge_group": edge,
            "reportable_findings": len(rep),
            "independent_reportable_findings": len(indep),
            "overlap_caution_findings": len(caution),
            "independent_datasets": ";".join(sorted(indep["dataset"].unique())) or "none",
            "expanded_layer_independent_findings": len(indep.loc[indep["evidence_layer"] == REFERENCE_LAYER]),
            "provisional_layer_independent_findings": len(indep.loc[indep["evidence_layer"] == "provisional_tumor_like"]),
            "decision": status_map[edge],
            "interpretation": interpretations[edge],
        })
    decisions = pd.DataFrame(decision_rows)

    checks = []
    def check(name: str, passed: bool, observed: object, expected: object) -> None:
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "observed": observed, "expected": expected})
    check("source_row_count", len(source) == 380, len(source), 380)
    check("finding_keys_unique", not findings.duplicated(keys).any(), int(findings.duplicated(keys).sum()), 0)
    check("consolidated_finding_count", len(findings) == 190, len(findings), 190)
    check("reportable_finding_count", len(reportable_findings) == 10, len(reportable_findings), 10)
    check("independent_reportable_finding_count", len(independent_findings) == 5, len(independent_findings), 5)
    check("all_reportable_methods_direction_concordant", (reportable_findings["method_direction_concordance"] == "concordant").all(), int((reportable_findings["method_direction_concordance"] == "concordant").sum()), len(reportable_findings))
    check("no_reportable_opposing_direction", not (reportable_findings["direction_vs_reference"] == "opposing").any(), int((reportable_findings["direction_vs_reference"] == "opposing").sum()), 0)
    check("decision_groups_complete", list(decisions["edge_group"]) == order, ";".join(decisions["edge_group"]), ";".join(order))
    validation = pd.DataFrame(checks)
    if (validation["status"] != "PASS").any():
        atomic_tsv(validation, PART / "validation_checks_v1.tsv")
        raise RuntimeError("Consolidated finding validation failed")

    atomic_tsv(findings, PART / "consolidated_bridge_findings_v1.tsv")
    atomic_tsv(reportable_findings, PART / "reportable_bridge_findings_v1.tsv")
    atomic_tsv(independent_findings, PART / "independent_reportable_bridge_findings_v1.tsv")
    atomic_tsv(decisions, PART / "edge_decisions_by_finding_v1.tsv")
    atomic_tsv(validation, PART / "validation_checks_v1.tsv")
    metadata = {
        "status": "PASS",
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "source_method_rows": len(source),
        "consolidated_findings": len(findings),
        "reportable_findings": len(reportable_findings),
        "independent_reportable_findings": len(independent_findings),
        "finding_unit": "dataset x evidence layer x relation",
        "method_rule": "Primary and sensitivity scores are robustness checks and are not counted as independent replication.",
    }
    (PART / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    (PART / "PHASE53_CONSOLIDATED_FINDINGS_PASS_v1").write_text(
        "status=PASS\nvalidation=8/8\nconsolidated_findings=190\nreportable_findings=10\nindependent_reportable_findings=5\n",
        encoding="utf-8",
    )
    manifest_rows = []
    for path in sorted(PART.iterdir()):
        if path.is_file():
            manifest_rows.append({"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)})
    atomic_tsv(pd.DataFrame(manifest_rows), PART / "FILE_MANIFEST_SHA256_v1.tsv")
    os.replace(PART, FINAL)
    print(f"PASS findings={len(findings)} reportable={len(reportable_findings)} independent={len(independent_findings)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
