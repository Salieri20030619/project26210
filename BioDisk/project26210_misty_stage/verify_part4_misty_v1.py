#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from scipy.spatial import cKDTree


PROJECT_ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
RESULT_ROOT = PROJECT_ROOT / "results/part4_spatial"
INPUT_DIR = RESULT_ROOT / "advanced_patch_misty_inputs_v1"
CROSS_V1 = RESULT_ROOT / "advanced_patch_misty_cross_sample_v1"
CROSS_V2 = RESULT_ROOT / "advanced_patch_misty_cross_sample_v2"
CROSS_V3 = RESULT_ROOT / "advanced_patch_misty_cross_sample_v3"
QC_DIR = RESULT_ROOT / "advanced_patch_misty_qc_v1"
MANUAL_REVIEW = QC_DIR / "manual_figure_review_v1.tsv"
SAMPLES = {"GLMF1": 3199, "GLMF2": 3880}
VIEWS = {"intra", "juxta.110", "para.200"}


class Verifier:
    def __init__(self) -> None:
        self.rows: list[dict[str, str]] = []

    def check(self, name: str, condition: bool, detail: str) -> None:
        self.rows.append(
            {"check_id": f"MISTY_QC_{len(self.rows) + 1:03d}", "check": name,
             "status": "PASS" if condition else "FAIL", "detail": detail}
        )
        if not condition:
            raise AssertionError(f"{name}: {detail}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sample_dir(sample: str) -> Path:
    return RESULT_ROOT / f"advanced_patch_misty_{sample.lower()}_v1"


def verify_inputs(verifier: Verifier) -> None:
    manifest = pd.read_csv(INPUT_DIR / "misty_feature_manifest_v1.tsv", sep="\t")
    verifier.check("feature_manifest_rows", manifest.shape[0] == 20, f"rows={manifest.shape[0]}")
    verifier.check("feature_ids_unique", manifest["feature_id"].is_unique, "20 unique feature IDs")
    verifier.check(
        "feature_classes",
        manifest["feature_class"].value_counts().to_dict()
        == {"cell2location_q05_abundance": 14, "prespecified_gene_program": 5, "replicated_hotspot_program": 1},
        str(manifest["feature_class"].value_counts().to_dict()),
    )
    stromal = manifest.loc[manifest["feature_id"] == "program_shared_stromal", "genes"].iloc[0].split(";")
    verifier.check("shared_stromal_genes", len(stromal) == len(set(stromal)) == 37, f"genes={len(stromal)}")
    metadata = json.loads((INPUT_DIR / "misty_input_metadata_v1.json").read_text())
    verifier.check("input_metadata_status", metadata["status"] == "complete", str(metadata["status"]))
    verifier.check(
        "spatial_parameters",
        metadata["juxtaview_threshold_um"] == 110.0
        and metadata["paraview_l_um"] == 200.0
        and metadata["paraview_zoi_um"] == 110.0
        and metadata["paraview_family"] == "gaussian",
        "juxta=110 um; para Gaussian l=200 um zoi=110 um",
    )
    for sample, expected_spots in SAMPLES.items():
        measurements = pd.read_csv(INPUT_DIR / f"{sample}_misty_measurements_v1.tsv", sep="\t")
        positions = pd.read_csv(INPUT_DIR / f"{sample}_misty_positions_v1.tsv", sep="\t")
        verifier.check(
            f"{sample}_input_shape",
            measurements.shape == (expected_spots, 21) and positions.shape == (expected_spots, 3),
            f"measurements={measurements.shape}; positions={positions.shape}",
        )
        verifier.check(
            f"{sample}_barcode_alignment",
            measurements["barcode"].is_unique and measurements["barcode"].equals(positions["barcode"]),
            f"barcodes={expected_spots}",
        )
        matrix = measurements.iloc[:, 1:].to_numpy(dtype=float)
        verifier.check(f"{sample}_finite_features", np.isfinite(matrix).all(), "all finite")
        max_mean = float(np.max(np.abs(matrix.mean(axis=0))))
        max_sd_error = float(np.max(np.abs(matrix.std(axis=0, ddof=0) - 1)))
        verifier.check(
            f"{sample}_standardized_features",
            max_mean < 1e-8 and max_sd_error < 1e-8,
            f"max_abs_mean={max_mean:.3g}; max_population_sd_error={max_sd_error:.3g}",
        )
        coords = positions[["x_um", "y_um"]].to_numpy(dtype=float)
        tree = cKDTree(coords)
        counts = np.array([len(x) - 1 for x in tree.query_ball_point(coords, r=110.0)])
        verifier.check(
            f"{sample}_six_neighbor_geometry",
            np.median(counts) == 6 and counts.max() <= 6,
            f"min={counts.min()}; median={np.median(counts):.0f}; max={counts.max()}",
        )


def verify_sample_results(verifier: Verifier) -> None:
    for sample, expected_spots in SAMPLES.items():
        directory = sample_dir(sample)
        status = json.loads((directory / "status.json").read_text())
        metadata = json.loads((directory / "run_metadata_v1.json").read_text())
        verifier.check(f"{sample}_status", status["status"] == "complete", status["status"])
        verifier.check(
            f"{sample}_run_metadata",
            metadata["status"] == "complete"
            and metadata["n_spots"] == expected_spots
            and metadata["n_features"] == 20
            and metadata["cv_folds"] == 10
            and metadata["random_forest"]["num_trees"] == 100,
            f"spots={metadata['n_spots']}; features={metadata['n_features']}; folds={metadata['cv_folds']}",
        )
        raw_performance = pd.read_csv(directory / "raw_misty/performance.txt", sep=r"\s+")
        raw_coefficients = pd.read_csv(directory / "raw_misty/coefficients.txt", sep=r"\s+")
        raw_importances = list((directory / "raw_misty").glob("importances_*.txt"))
        verifier.check(
            f"{sample}_raw_outputs",
            raw_performance.shape == (20, 7) and raw_coefficients.shape == (20, 9) and len(raw_importances) == 60,
            f"performance={raw_performance.shape}; coefficients={raw_coefficients.shape}; importance_files={len(raw_importances)}",
        )
        tables = {
            "improvements": (f"{sample}_misty_improvements_v1.tsv", (160, 4)),
            "contributions": (f"{sample}_misty_contributions_v1.tsv", (160, 4)),
            "importances": (f"{sample}_misty_importances_v1.tsv", (1200, 5)),
            "improvement_stats": (f"{sample}_misty_improvement_stats_v1.tsv", (120, 5)),
            "contribution_stats": (f"{sample}_misty_contribution_stats_v1.tsv", (60, 6)),
            "importance_aggregated": (f"{sample}_misty_importances_aggregated_v1.tsv", (1200, 6)),
        }
        loaded = {}
        for label, (filename, expected_shape) in tables.items():
            frame = pd.read_csv(directory / filename, sep="\t")
            loaded[label] = frame
            verifier.check(f"{sample}_{label}_shape", frame.shape == expected_shape, f"shape={frame.shape}")
        verifier.check(
            f"{sample}_view_sets",
            set(loaded["contribution_stats"]["view"]) == VIEWS
            and set(loaded["importances"]["view"]) == VIEWS,
            f"views={sorted(set(loaded['importances']['view']))}",
        )
        self_rows = loaded["importances"]["Predictor"] == loaded["importances"]["Target"]
        verifier.check(
            f"{sample}_importance_nan_pattern",
            loaded["importances"].loc[self_rows, "Importance"].isna().all()
            and loaded["importances"].loc[~self_rows, "Importance"].notna().all()
            and int(self_rows.sum()) == 60,
            f"self_rows={int(self_rows.sum())}",
        )
        gains = loaded["improvement_stats"].pivot(index="target", columns="measure", values="mean")
        positive = int(((gains["gain.R2"] > 0) & (gains["gain.RMSE"] > 0)).sum())
        expected_positive = 18 if sample == "GLMF1" else 19
        verifier.check(f"{sample}_positive_gain_count", positive == expected_positive, f"positive={positive}/20")
        for stem in ("performance_gains", "view_contributions", "juxta_para_importance"):
            for suffix in ("png", "pdf"):
                path = directory / f"{sample}_misty_{stem}_v1.{suffix}"
                verifier.check(f"{sample}_{stem}_{suffix}", path.is_file() and path.stat().st_size > 1000, f"bytes={path.stat().st_size}")


def verify_cross_results(verifier: Verifier) -> None:
    metadata = json.loads((CROSS_V3 / "cross_sample_run_metadata_v3.json").read_text())
    verifier.check("cross_status", metadata["status"] == "complete", metadata["status"])
    verifier.check("cross_positive_gain_count", metadata["positive_gain_both_sections"] == 17, "positive=17/20")
    verifier.check(
        "cross_replicated_counts",
        metadata["replicated_top_decile_counts"] == {"intra": 20, "juxta.110": 15, "para.200": 13},
        str(metadata["replicated_top_decile_counts"]),
    )
    files = {
        "performance": ("cross_sample_misty_performance_v3.tsv", (20, 20)),
        "contributions": ("cross_sample_misty_view_contributions_v3.tsv", (60, 12)),
        "importance": ("cross_sample_misty_importance_support_v3.tsv", (1140, 12)),
        "replicated": ("cross_sample_misty_replicated_top_decile_v3.tsv", (48, 12)),
    }
    frames = {}
    for label, (filename, expected_shape) in files.items():
        frame = pd.read_csv(CROSS_V3 / filename, sep="\t")
        frames[label] = frame
        verifier.check(f"cross_{label}_shape", frame.shape == expected_shape, f"shape={frame.shape}")
    recomputed = frames["importance"]["replicated_top_decile"].astype(bool)
    verifier.check(
        "cross_replicated_table_exact",
        int(recomputed.sum()) == 48
        and frames["replicated"][["view", "Predictor", "Target"]].apply(tuple, axis=1).is_unique,
        f"replicated={int(recomputed.sum())}",
    )
    for base in ("performance", "view_contributions", "importance_support", "replicated_top_decile"):
        hashes = [
            sha256(CROSS_V1 / f"cross_sample_misty_{base}_v1.tsv"),
            sha256(CROSS_V2 / f"cross_sample_misty_{base}_v2.tsv"),
            sha256(CROSS_V3 / f"cross_sample_misty_{base}_v3.tsv"),
        ]
        verifier.check(f"cross_{base}_v1_v2_v3_identical", len(set(hashes)) == 1, hashes[0])
    for stem in ("performance_gains", "view_contributions", "importance_replication"):
        for suffix in ("png", "pdf"):
            path = CROSS_V3 / f"cross_sample_misty_{stem}_v3.{suffix}"
            verifier.check(f"cross_{stem}_{suffix}", path.is_file() and path.stat().st_size > 1000, f"bytes={path.stat().st_size}")


def verify_manual_figures(verifier: Verifier) -> None:
    review = pd.read_csv(MANUAL_REVIEW, sep="\t")
    final = review.loc[review["artifact_scope"] == "formal_final"]
    verifier.check("manual_review_final_rows", final.shape[0] == 9, f"rows={final.shape[0]}")
    verifier.check("manual_review_final_pass", (final["review_status"] == "pass").all(), "final PNG pass=9/9")
    verifier.check(
        "manual_review_v2_failure_recorded",
        ((review["artifact_scope"] == "superseded_audit") & (review["review_status"] == "fail")).sum() == 1,
        "one superseded v2 label-overlap failure",
    )
    for relative in final["figure_path"]:
        path = PROJECT_ROOT / relative
        with Image.open(path) as image:
            verifier.check(
                f"png_dimensions_{path.name}",
                image.width >= 1800 and image.height >= 1400,
                f"{image.width}x{image.height}",
            )


def write_manifest() -> int:
    roots = [INPUT_DIR, sample_dir("GLMF1"), sample_dir("GLMF2"), CROSS_V1, CROSS_V2, CROSS_V3, QC_DIR]
    standalone = [
        RESULT_ROOT / "prepare_part4_misty_inputs_v1.py",
        RESULT_ROOT / "run_part4_misty_sample_v1.R",
        RESULT_ROOT / "run_part4_misty_pipeline_v1.sh",
        RESULT_ROOT / "summarize_part4_misty_cross_sample_v1.py",
        RESULT_ROOT / "summarize_part4_misty_cross_sample_v2.py",
        RESULT_ROOT / "summarize_part4_misty_cross_sample_v3.py",
        RESULT_ROOT / "verify_part4_misty_v1.py",
        RESULT_ROOT / "misty_pipeline_v1.log",
    ]
    files: set[Path] = set()
    for root in roots:
        if root.exists():
            files.update(path for path in root.rglob("*") if path.is_file())
    files.update(path for path in standalone if path.is_file())
    rows = []
    for path in sorted(files):
        if path.name == "file_manifest_v1.tsv":
            continue
        rows.append(
            {"relative_path": str(path.relative_to(PROJECT_ROOT)), "size_bytes": path.stat().st_size, "sha256": sha256(path)}
        )
    pd.DataFrame(rows).to_csv(QC_DIR / "file_manifest_v1.tsv", sep="\t", index=False)
    return len(rows)


def main() -> None:
    QC_DIR.mkdir(parents=True, exist_ok=True)
    verifier = Verifier()
    verify_inputs(verifier)
    verify_sample_results(verifier)
    verify_cross_results(verifier)
    verify_manual_figures(verifier)
    qc_path = QC_DIR / "verification_qc_v1.tsv"
    pd.DataFrame(verifier.rows).to_csv(qc_path, sep="\t", index=False)
    file_count = write_manifest()
    summary = {
        "status": "complete",
        "checks_passed": len(verifier.rows),
        "checks_failed": 0,
        "manifest_files": file_count,
        "formal_png_manual_review_passed": 9,
        "formal_png_manual_review_failed": 0,
        "superseded_visual_failures_recorded": 1,
        "cross_sample_positive_multiview_gain_targets": 17,
        "replicated_top_decile_pairs": {"intra": 20, "juxta.110": 15, "para.200": 13},
        "interpretation_limit": "predictive multiview associations only; importance is neither signed effect nor causality",
    }
    (QC_DIR / "verification_summary_v1.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"MISTY_VERIFICATION_ALL_PASS checks={len(verifier.rows)} files={file_count} final_png=9")


if __name__ == "__main__":
    main()
