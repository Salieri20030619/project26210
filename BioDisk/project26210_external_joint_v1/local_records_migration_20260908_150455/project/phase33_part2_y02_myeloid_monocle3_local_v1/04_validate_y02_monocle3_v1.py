#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(); parser.add_argument("--analysis", type=Path, required=True); parser.add_argument("--output", type=Path, required=True); args = parser.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    args.output.mkdir(parents=True)
    status = json.loads((args.analysis / "status.json").read_text())
    cell = pd.read_csv(args.analysis / "y02_monocle3_dpt_cell_scores_v1.tsv.gz", sep="\t")
    effects = pd.read_csv(args.analysis / "y02_monocle3_patient_effects_v1.tsv", sep="\t")
    summary = pd.read_csv(args.analysis / "y02_monocle3_primary_summary_v1.tsv", sep="\t")
    conc = pd.read_csv(args.analysis / "y02_monocle3_dpt_patient_concordance_v1.tsv", sep="\t")
    state = pd.read_csv(args.analysis / "y02_monocle3_state_summary_v1.tsv", sep="\t")
    checks = [
        ("status", status["status"] == "PASS_Y02_MONOCLE3_ANALYSIS_COMPLETE_PRESENTATION_PENDING", status["status"]),
        ("cell_scope", len(cell) == 50761 and cell["cell_id"].nunique() == 50761, str(len(cell))),
        ("study_scope", cell["study_v1"].nunique() == 5, str(cell["study_v1"].nunique())),
        ("patient_scope", cell["patient"].nunique() == 40, str(cell["patient"].nunique())),
        ("finite_pseudotime", np.isfinite(cell["monocle3_pseudotime_v1"]).all() and cell["monocle3_pseudotime_v1"].between(0, 1).all(), "all finite in [0,1]"),
        ("patient_endpoint_grid", not effects.duplicated(["study", "patient", "endpoint"]).any() and set(effects["endpoint"]) == {"monocyte_program", "spp1_state_no_spp1", "spp1_vs_c1qc", "SPP1_minus_Monocyte_pseudotime"}, str(len(effects))),
        ("summary_grid", len(summary) == 20 and not summary.duplicated(["study", "endpoint"]).any(), str(len(summary))),
        ("concordance_units", not conc.duplicated(["study", "patient"]).any(), str(len(conc))),
        ("state_grid", len(state) == 20, str(len(state))),
        ("two_sided_no_direction", status["two_sided_testing"] is True and status["direction_prespecified"] is False, "two-sided; no direction prespecified"),
        ("reporting_policy", set(summary["reporting_tier"]).issubset({"FDR_AT_LEAST_10", "NOMINAL_P_3_TO_9", "DIRECTION_ONLY_LT3"}), ",".join(sorted(summary["reporting_tier"].unique()))),
        ("interpretation_guard", "not lineage tracing" in status["interpretation"], status["interpretation"]),
        ("no_figures", status["figures_created"] is False and not any(args.analysis.rglob("*.png")) and not any(args.analysis.rglob("*.pdf")), "no figures"),
    ]
    failed = [name for name, ok, _ in checks if not ok]
    with (args.output / "verification_checks_v1.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n"); writer.writerow(["check", "pass", "detail"]); writer.writerows((name, str(ok), detail) for name, ok, detail in checks)
    payload = {"status": "PASS_INDEPENDENT_VALIDATION" if not failed else "FAIL_INDEPENDENT_VALIDATION", "checks_passed": len(checks)-len(failed), "checks_total": len(checks), "failed_checks": failed, "analysis_file_sha256": {p.name: sha256(p) for p in sorted(args.analysis.iterdir()) if p.is_file()}}
    (args.output / "status.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    if failed: raise RuntimeError(f"Validation failed: {failed}")


if __name__ == "__main__": main()
