#!/usr/bin/env python3
"""Build the deduplicated, overlap-aware final evidence tables for Phase 53."""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
ALL_RESULTS = HERE / "analysis_local_v2" / "bridge_screen_all_exact_v2.tsv"
OVERLAP = HERE / "signature_overlap_audit_local_v2" / "bridge_relation_signature_overlap_v1.tsv"
FINAL = HERE / "final_evidence_local_v1"
PART = HERE / "final_evidence_local_v1.part"
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
        raise FileExistsError("Refusing to overwrite an existing final-evidence directory")
    PART.mkdir(parents=True)

    results = pd.read_csv(ALL_RESULTS, sep="\t")
    overlap = pd.read_csv(OVERLAP, sep="\t")
    merged = results.merge(
        overlap[["relation", "overlap_count", "overlap_genes", "jaccard", "overlap_flag", "cross_compartment"]],
        on="relation",
        how="left",
        validate="many_to_one",
    )
    if merged["overlap_flag"].isna().any():
        raise RuntimeError("One or more relations are missing from the signature-overlap audit")

    malignant_dependent = (
        merged["x_compartment"].astype(str).str.startswith("malignant")
        | merged["y_compartment"].astype(str).str.startswith("malignant")
    )
    merged["malignancy_layer_relevant"] = malignant_dependent
    keep = malignant_dependent | (merged["evidence_layer"] == REFERENCE_LAYER)
    final = merged.loc[keep].copy()
    final["deduplication_rule"] = "retain_each_malignancy_layer"
    final.loc[~final["malignancy_layer_relevant"], "deduplication_rule"] = (
        "nonmalignant_edge_once_in_expanded_reference_layer"
    )
    final["evidence_weight"] = "independent_gene_sets"
    final.loc[final["overlap_flag"] != "INDEPENDENT_GENE_SETS", "evidence_weight"] = (
        "state_consistency_only_due_to_signature_overlap"
    )

    reportable = final.loc[final["reportable_v2"].astype(str).str.lower() == "true"].copy()
    independent = reportable.loc[reportable["overlap_flag"] == "INDEPENDENT_GENE_SETS"].copy()

    order = [
        "tumor_to_caf",
        "caf_to_tam",
        "tam_to_tcell",
        "tam_to_vascular",
        "tumor_to_tam",
        "tam_state_internal",
    ]
    decision_rows = []
    decision_text = {
        "tumor_to_caf": (
            "NO_REPORTABLE_SUPPORT",
            "No reportable F1/TYMP-to-POSTN CAF association; patient overlap is sparse in every cohort.",
        ),
        "caf_to_tam": (
            "SMALL_N_INDEPENDENT_SUPPORT",
            "POSTN CAF-to-APOE TAM is positive in GSE264586 (n=5) by exact permutation; external replication is absent.",
        ),
        "tam_to_tcell": (
            "NO_REPORTABLE_SUPPORT",
            "No exact/FDR-reportable TAM-to-T-cell dysfunction or cytotoxicity relation in the full-cohort screen.",
        ),
        "tam_to_vascular": (
            "NO_REPORTABLE_SUPPORT",
            "No reportable TAM-to-vascular relation; the tested APOE/SPP1 vascular pairs also share LGALS3.",
        ),
        "tumor_to_tam": (
            "PROVISIONAL_LAYER_ONLY_SUPPORT",
            "TYMP-to-SPP1 in GSE178481 and F1-to-C1QC in GSE207493 appear only in the broad provisional tumour-like layer, not the expanded CNV-positive layer.",
        ),
        "tam_state_internal": (
            "PARTIAL_INTERNAL_STATE_SUPPORT",
            "SPP1-to-myeloid-suppression has independent-gene-set support in two cohorts; APOE-to-suppression is strong but partly circular because four genes overlap.",
        ),
    }
    for edge in order:
        subset = final.loc[final["edge_group"] == edge]
        rep = reportable.loc[reportable["edge_group"] == edge]
        indep = independent.loc[independent["edge_group"] == edge]
        status, interpretation = decision_text[edge]
        decision_rows.append(
            {
                "edge_group": edge,
                "tested_rows_after_deduplication": len(subset),
                "reportable_rows": len(rep),
                "independent_reportable_rows": len(indep),
                "independent_reportable_datasets": ";".join(sorted(indep["dataset"].unique())) or "none",
                "expanded_layer_independent_reportable_rows": len(
                    indep.loc[indep["evidence_layer"] == REFERENCE_LAYER]
                ),
                "provisional_layer_independent_reportable_rows": len(
                    indep.loc[indep["evidence_layer"] == "provisional_tumor_like"]
                ),
                "decision": status,
                "interpretation": interpretation,
            }
        )
    decisions = pd.DataFrame(decision_rows)

    checks = []
    def check(name: str, passed: bool, observed: object, expected: object) -> None:
        checks.append({
            "check": name,
            "status": "PASS" if passed else "FAIL",
            "observed": observed,
            "expected": expected,
        })

    nonmal = final.loc[~final["malignancy_layer_relevant"]]
    nonmal_keys = ["dataset", "method", "relation"]
    check("source_row_count", len(results) == 660, len(results), 660)
    check("deduplicated_row_count", len(final) == 380, len(final), 380)
    check("nonmalignant_edges_unique_after_deduplication", not nonmal.duplicated(nonmal_keys).any(), int(nonmal.duplicated(nonmal_keys).sum()), 0)
    check("nonmalignant_edges_use_reference_layer", (nonmal["evidence_layer"] == REFERENCE_LAYER).all(), ";".join(sorted(nonmal["evidence_layer"].unique())), REFERENCE_LAYER)
    check("reportable_row_count", len(reportable) == 15, len(reportable), 15)
    check("independent_reportable_row_count", len(independent) == 6, len(independent), 6)
    check("small_n_reportable_has_exact_p", reportable.loc[reportable["n_patients"].between(5, 9), "small_n_exact_permutation_p"].notna().all(), int(reportable.loc[reportable["n_patients"].between(5, 9), "small_n_exact_permutation_p"].notna().sum()), int(reportable["n_patients"].between(5, 9).sum()))
    check("no_reportable_opposing_direction", not (reportable["direction_vs_reference"] == "opposing").any(), int((reportable["direction_vs_reference"] == "opposing").sum()), 0)
    check("expanded_tumor_to_tam_independent_reportable", len(independent.loc[(independent["edge_group"] == "tumor_to_tam") & (independent["evidence_layer"] == REFERENCE_LAYER)]) == 0, len(independent.loc[(independent["edge_group"] == "tumor_to_tam") & (independent["evidence_layer"] == REFERENCE_LAYER)]), 0)
    check("decision_groups_complete", list(decisions["edge_group"]) == order, ";".join(decisions["edge_group"]), ";".join(order))
    validation = pd.DataFrame(checks)
    if (validation["status"] != "PASS").any():
        atomic_tsv(validation, PART / "validation_checks_v1.tsv")
        raise RuntimeError("Final evidence validation failed")

    atomic_tsv(final, PART / "final_bridge_evidence_v1.tsv")
    atomic_tsv(reportable, PART / "final_reportable_bridge_evidence_v1.tsv")
    atomic_tsv(independent, PART / "final_independent_reportable_evidence_v1.tsv")
    atomic_tsv(decisions, PART / "final_edge_decisions_v1.tsv")
    atomic_tsv(validation, PART / "validation_checks_v1.tsv")

    metadata = {
        "status": "PASS",
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "source_rows": len(results),
        "deduplicated_rows": len(final),
        "reportable_rows": len(reportable),
        "independent_reportable_rows": len(independent),
        "reference_layer_for_nonmalignant_edges": REFERENCE_LAYER,
        "small_n_rule": "n=5-9: full-permutation exact two-sided Spearman p < 0.05; n<5 direction only; n>=10 BH FDR < 0.05",
    }
    (PART / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    (PART / "PHASE53_FINAL_EVIDENCE_PASS_v1").write_text(
        "status=PASS\nvalidation=10/10\nsource_rows=660\ndeduplicated_rows=380\nreportable_rows=15\nindependent_reportable_rows=6\n",
        encoding="utf-8",
    )

    manifest_rows = []
    for path in sorted(PART.iterdir()):
        if path.is_file():
            manifest_rows.append({"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)})
    atomic_tsv(pd.DataFrame(manifest_rows), PART / "FILE_MANIFEST_SHA256_v1.tsv")
    os.replace(PART, FINAL)
    print(f"PASS source={len(results)} final={len(final)} reportable={len(reportable)} independent={len(independent)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
