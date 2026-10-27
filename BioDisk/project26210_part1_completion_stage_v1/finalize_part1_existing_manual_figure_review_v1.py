#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path

from PIL import Image
import pandas as pd


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def inspect(path: Path, status: str, selected: bool, note: str) -> dict[str, object]:
    with Image.open(path) as image:
        width, height = image.size
        image.verify()
    return {
        "figure": path.name,
        "relative_path": str(path.relative_to(PROJECT)),
        "status": status,
        "selected_for_use": selected,
        "width_px": width,
        "height_px": height,
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
        "manual_checks": "content;orientation;color_scale;labels;legend;clipping;readability",
        "note": note,
    }


def update_json(path: Path, **updates: object) -> None:
    data = json.loads(path.read_text(encoding="utf-8"))
    data.update(updates)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main() -> None:
    completed_at = datetime.now().astimezone().isoformat(timespec="seconds")
    clone_dir = ROOT / "04_subclones_tymp/cnv_clone_tymp_cnmf_association_v1"
    clone_qc = ROOT / "04_subclones_tymp/cnv_clone_tymp_cnmf_association_verification_v1"
    clone_notes = {
        "factor_contrast_threshold_sensitivity_v1.png": "PASS: all seven factors and three thresholds visible; zero line and legend readable.",
        "formal_primary_clone_tymp_landscape_v1.png": "PASS: patient rows, clone ranks, point sizes and TYMP color scale visible without clipping.",
        "selected_F1_patient_tymp_deltas_v1.png": "PASS: patient labels, zero line and threshold-sensitive red bars readable.",
        "selected_F1_threshold_patient_sensitivity_v1.png": "PASS: all thresholds, points, medians and zero line visible.",
    }
    clone_rows = [inspect(clone_dir / name, "PASS", True, note) for name, note in clone_notes.items()]
    pd.DataFrame(clone_rows).to_csv(clone_qc / "manual_figure_review_v1.tsv", sep="\t", index=False)
    (clone_qc / "MANUAL_FIGURE_REVIEW_COMPLETE").write_text(completed_at + "\n", encoding="utf-8")
    update_json(
        clone_qc / "verification_summary_v1.json",
        manual_figure_review_pending=False,
        manual_figure_review_completed_at=completed_at,
        manual_figure_review_file="manual_figure_review_v1.tsv",
        n_manual_figures_pass=4,
        n_manual_figures_fail=0,
    )

    traj = ROOT / "09_metabolism_trajectory/trajectory_v1"
    traj_notes = {
        "Monocle3_pseudotime_v1.png": "PASS: complete embedding and labeled pseudotime color scale; no clipping.",
        "CytoTRACE2_potency_v1.png": "PASS: complete embedding, axes and color scale; no clipping.",
        "scTour_pseudotime_v1.png": "PASS: complete embedding, axes and color scale; no clipping.",
        "CellRank_fate_confidence_v1.png": "PASS: complete embedding, axes and probability color scale; no clipping.",
        "CellRank_fate_Cycle_v1.png": "PASS: complete embedding, axes and probability color scale; no clipping.",
        "CellRank_fate_Hypoxia_v1.png": "PASS: complete embedding, axes and probability color scale; no clipping.",
        "CellRank_fate_TYMP_state_v1.png": "PASS: complete embedding, axes and probability color scale; no clipping.",
    }
    traj_rows = [inspect(traj / "figures" / name, "PASS", True, note) for name, note in traj_notes.items()]
    traj_rows.append(
        inspect(
            traj / "figures/Slingshot_pseudotime_v1.png",
            "FAIL",
            False,
            "FAIL: pseudotime color scale is absent; retained as immutable visual-failure audit.",
        )
    )
    traj_rows.append(
        inspect(
            traj / "figures/Slingshot_pseudotime_v2.png",
            "PASS",
            True,
            "PASS: v2 adds labeled pseudotime color scale and explicit UMAP axes; values were not recomputed.",
        )
    )
    pd.DataFrame(traj_rows).to_csv(traj / "qc/manual_figure_review_v1.tsv", sep="\t", index=False)
    (traj / "qc/MANUAL_FIGURE_REVIEW_COMPLETE").write_text(completed_at + "\n", encoding="utf-8")
    update_json(
        traj / "qc/trajectory_verification_metadata_v1.json",
        manual_figure_review_pending=False,
        manual_figure_review_completed_at=completed_at,
        manual_figure_review_file="qc/manual_figure_review_v1.tsv",
        n_selected_manual_figures_pass=8,
        n_archived_manual_figures_fail=1,
        selected_slingshot_figure="figures/Slingshot_pseudotime_v2.png",
    )
    print("PART1_EXISTING_MANUAL_FIGURE_REVIEW_V1_PASS clone=4 trajectory_selected=8 archived_fail=1")


if __name__ == "__main__":
    main()
