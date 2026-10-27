#!/usr/bin/env python3
"""Independent structural and numerical verification for TYMP-SPP1 state v1."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from scipy.stats import spearmanr


EXPECTED_SAMPLES = {"H46_T", "J38_T", "R114_T", "R29_T", "R51_T", "S15_T", "X49_T", "X98_T", "Y12_T", "Y27_T", "Y7_T", "Z43_T"}


def sha256(path: Path, chunk: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while data := handle.read(chunk):
            digest.update(data)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    root = args.output
    checks = []

    def check(name: str, passed: bool, detail: str = "") -> None:
        checks.append({"check": name, "passed": bool(passed), "detail": detail})

    metadata = json.loads((root / "run_metadata_v1.json").read_text())
    check("metadata_completed", metadata.get("status") == "completed", str(metadata.get("status")))
    check("metadata_patients", metadata.get("n_patients") == 12, str(metadata.get("n_patients")))
    check("metadata_spots", metadata.get("n_spots") == 34609, str(metadata.get("n_spots")))
    check("metadata_samples", set(metadata.get("samples", [])) == EXPECTED_SAMPLES)

    fingerprints = pd.read_csv(root / "tables/input_fingerprints_v1.tsv", sep="\t")
    check("fingerprint_samples", set(fingerprints["sample"]) == EXPECTED_SAMPLES)
    for row in fingerprints.itertuples():
        path = Path(row.path)
        check(f"input_exists_{row.sample}", path.is_file())
        check(f"input_size_{row.sample}", path.stat().st_size == row.size_bytes)
        check(f"input_hash_{row.sample}", sha256(path) == row.sha256)

    scores = pd.read_csv(root / "tables/patient_spot_state_scores_v1.tsv.gz", sep="\t")
    check("score_rows", len(scores) == 34609, str(len(scores)))
    check("score_samples", set(scores["sample"]) == EXPECTED_SAMPLES)
    check("score_unique_ids", not scores.duplicated(["sample", "barcode"]).any())
    numeric = scores.select_dtypes(include=[np.number])
    check("score_finite", np.isfinite(numeric.to_numpy()).all())
    check("macro_flags_balanced_refdep", scores.groupby("sample").macro_enriched_reference_dependent.mean().between(0.49, 0.51).all())
    check("macro_flags_balanced_reflight", scores.groupby("sample").macro_enriched_reference_light.mean().between(0.49, 0.51).all())

    tests = pd.read_csv(root / "tables/patient_state_associations_v1.tsv", sep="\t")
    check("association_rows", len(tests) == 12 * 2 * 2 * 2 * 3 * 2, str(len(tests)))
    check("association_samples", set(tests["sample"]) == EXPECTED_SAMPLES)
    check("association_finite", np.isfinite(tests[["rho", "p_value", "within_patient_fdr_bh"]].to_numpy()).all())
    check("association_ranges", tests.rho.between(-1, 1).all() and tests.p_value.between(0, 1).all() and tests.within_patient_fdr_bh.between(0, 1).all())

    primary = tests[(tests.route == "reference_dependent") & (tests.adjustment == "composition_qc") &
                    (tests.subset == "macrophage_enriched") & (tests.outcome == "SPP1_state_no_SPP1")]
    check("primary_24_rows", len(primary) == 24)
    sample = "H46_T"; grid = "4x4"
    score = scores[scores["sample"] == sample]
    score = score[score.macro_enriched_reference_dependent]
    block = score.groupby("block_4x4", as_index=False).agg(
        x=("neighbor_THM_refdep", "median"),
        y=("SPP1_state_no_SPP1_resid_reference_dependent_composition_qc", "median"),
    )
    recomputed = float(spearmanr(block.x, block.y).statistic)
    recorded = float(primary[(primary["sample"] == sample) & (primary.grid == grid)].rho.iloc[0])
    check("recompute_primary_rho", np.isclose(recomputed, recorded, atol=1e-12), f"{recomputed} vs {recorded}")

    summary = pd.read_csv(root / "tables/cohort_state_direction_summary_v1.tsv", sep="\t")
    check("summary_rows", len(summary) == 2 * 2 * 2 * 3 * 2, str(len(summary)))
    check("summary_patient_counts", (summary.patients_tested == 12).all())
    check("summary_direction_counts", ((summary.positive_patients + summary.negative_patients) == 12).all())
    check("summary_fdr_range", summary.sign_test_fdr_bh.between(0, 1).all())

    mediation = pd.read_csv(root / "tables/patient_exploratory_mediation_v1.tsv", sep="\t")
    check("mediation_rows", len(mediation) == 48, str(len(mediation)))
    check("mediation_finite", np.isfinite(mediation.select_dtypes(include=[np.number]).to_numpy()).all())
    downstream = pd.read_csv(root / "tables/patient_interface_downstream_tests_v1.tsv", sep="\t")
    check("downstream_rows", len(downstream) == 12 * 5 * 2, str(len(downstream)))
    check("downstream_finite", np.isfinite(downstream[["rho", "p_value", "within_patient_fdr_bh"]].to_numpy()).all())

    memory = pd.read_csv(root / "memory_gates_v1.tsv", sep="\t")
    check("memory_gates_all_pass", (memory.exit_code == 0).all(), f"n={len(memory)}")
    check("memory_gate_per_sample", set(memory.loc[memory.stage == "sample_analysis", "sample"]) == EXPECTED_SAMPLES)

    expected_png = 12 + 4
    pngs = sorted((root / "figures").rglob("*.png"))
    pdfs = sorted((root / "figures").rglob("*.pdf"))
    check("png_count", len(pngs) == expected_png, str(len(pngs)))
    check("pdf_count", len(pdfs) == expected_png, str(len(pdfs)))
    for path in pngs:
        with Image.open(path) as image:
            check(f"png_valid_{path.name}", image.width >= 1000 and image.height >= 800, f"{image.size}")

    manifest = pd.read_csv(root / "file_manifest_v1.tsv", sep="\t")
    check("file_manifest_unique", manifest.path.is_unique)
    for row in manifest.itertuples():
        path = root / row.path
        check(f"manifest_exists_{row.Index}", path.is_file())
        check(f"manifest_size_{row.Index}", path.stat().st_size == row.size_bytes)
        check(f"manifest_hash_{row.Index}", sha256(path) == row.sha256)

    qc = pd.DataFrame(checks)
    qc.to_csv(root / "verification_checks_v1.tsv", sep="\t", index=False)
    metadata_out = {"status": "passed" if qc.passed.all() else "failed", "checks": len(qc),
                    "passed": int(qc.passed.sum()), "failed": int((~qc.passed).sum())}
    (root / "verification_metadata_v1.json").write_text(json.dumps(metadata_out, indent=2), encoding="utf-8")
    if not qc.passed.all():
        print(qc[~qc.passed].to_string(index=False))
        raise SystemExit(1)
    print(f"TYMP_SPP1_STATE_VERIFICATION_ALL_PASS checks={len(qc)}")


if __name__ == "__main__":
    main()
