#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results_local_v1"
VALIDATION = ROOT / "validation_local_v1"
OUT = ROOT / "terminal_validation_local_v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if OUT.exists():
        raise RuntimeError(f"Refusing to overwrite existing terminal directory: {OUT}")
    OUT.mkdir(parents=False)
    analysis = json.loads((RESULTS / "status.json").read_text(encoding="utf-8"))
    validation = json.loads((VALIDATION / "status.json").read_text(encoding="utf-8"))
    patients = pd.read_csv(RESULTS / "m03_patient_tree_summary_v1.tsv", sep="\t")
    nodes = pd.read_csv(RESULTS / "m03_tree_nodes_cnmf_annotation_v1.tsv", sep="\t")
    clone_nodes = nodes.loc[nodes["node_type"].eq("formal_primary_inferred_cnv_clone")]
    multi_factor_patients = int((patients["n_distinct_clone_dominant_factors"] >= 2).sum())
    checks = [
        ("analysis_complete", analysis["status"] == "PASS_M03_ANALYSIS_TABLES_COMPLETE_PRESENTATION_PENDING", analysis["status"]),
        ("independent_validation", validation["status"] == "PASS_INDEPENDENT_VALIDATION", validation["status"]),
        ("independent_checks", validation["checks_passed"] == 20 and validation["checks_total"] == 20, f"{validation['checks_passed']}/{validation['checks_total']}"),
        ("profile_reproduction", analysis["n_profile_reproductions"] == analysis["n_profile_hash_matches"] == 227, "227/227"),
        ("tree_scope", analysis["n_evaluable_patients"] == 14 and analysis["n_primary_clone_nodes"] == 43 and analysis["n_tree_edges"] == 43, "14 patients; 43 clones; 43 edges"),
        ("multi_factor_patient_sentinel", multi_factor_patients == 4, str(multi_factor_patients)),
        ("clone_factor_scope", clone_nodes["dominant_cnmf_factor"].notna().all(), str(len(clone_nodes))),
        ("interpretation_boundary", "not mutation-validated temporal ancestry" in analysis["interpretation_boundary"], analysis["interpretation_boundary"]),
        ("failed_attempt_preserved", (ROOT / "validation_failed_attempt1_assumed_48_nodes" / "status.json").is_file(), "validator-only fixed-count assumption"),
        ("no_figures", analysis["figures_created"] is False and not any(ROOT.rglob("*.png")) and not any(ROOT.rglob("*.pdf")), "no PNG/PDF"),
    ]
    failed = [name for name, passed, _ in checks if not passed]
    path = OUT / "terminal_checks_v1.tsv.part"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["check", "pass", "detail"])
        writer.writerows((name, str(passed), detail) for name, passed, detail in checks)
    path.replace(OUT / "terminal_checks_v1.tsv")

    files = []
    for path in sorted(list(RESULTS.glob("*")) + list(VALIDATION.glob("*"))):
        if path.is_file():
            files.append({"relative_path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "sha256": sha256(path)})
    manifest = OUT / "validated_file_manifest_v1.tsv.part"
    pd.DataFrame(files).to_csv(manifest, sep="\t", index=False, lineterminator="\n")
    manifest.replace(OUT / "validated_file_manifest_v1.tsv")
    payload = {
        "status": "PASS_PHASE32_M03_TERMINAL_VALIDATION" if not failed else "FAIL_PHASE32_M03_TERMINAL_VALIDATION",
        "timestamp_local": datetime.now().astimezone().isoformat(timespec="seconds"),
        "checks_passed": len(checks) - len(failed),
        "checks_total": len(checks),
        "independent_component_checks": "20/20",
        "failed_checks": failed,
        "presentation_status": "TREE_FIGURE_PENDING_USER_REVIEW",
    }
    status = OUT / "status.json.part"
    status.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    status.replace(OUT / "status.json")
    if failed:
        raise RuntimeError(f"Terminal validation failed: {failed}")


if __name__ == "__main__":
    main()
