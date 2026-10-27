#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import pandas as pd
from PIL import Image


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
STAGE = Path("/mnt/BioDisk/wangmingyuan/project26210_tymp_overnight_stage_v1")
FIG = STAGE / "joint_figures_v4"
QC = STAGE / "joint_qc_v2"
GATE = ROOT / "check_memory_gate.sh"
EXPECTED = [
    "vascular_niche_patient_effects_v1",
    "vascular_niche_cohort_summary_v1",
    "myeloid_cross_sectional_state_overview_v2",
    "myeloid_state_patient_directions_v1",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    output_tsv = QC / "joint_figure_verification_v1.tsv"
    output_json = QC / "joint_figure_verification_summary_v1.json"
    if output_tsv.exists() or output_json.exists():
        raise SystemExit("Refusing to overwrite final figure verification")
    gate = subprocess.run([str(GATE)], text=True, capture_output=True, check=False)
    gate_line = (gate.stdout or gate.stderr).strip()
    if gate.returncode != 0:
        raise RuntimeError(gate_line)
    prior = json.loads((QC / "independent_verification_summary_v2.json").read_text())
    if prior["status"] != "passed" or prior["checks_failed"] != 0:
        raise ValueError("Independent numeric verification did not pass")
    metadata = json.loads((FIG / "figure_generation_metadata_v1.json").read_text())
    rows = []
    all_pass = True
    for stem in EXPECTED:
        png = FIG / f"{stem}.png"
        pdf = FIG / f"{stem}.pdf"
        png_ok = png.is_file() and png.stat().st_size > 100_000
        pdf_ok = pdf.is_file() and pdf.stat().st_size > 10_000 and pdf.read_bytes()[:5] == b"%PDF-"
        width = height = 0
        if png_ok:
            with Image.open(png) as image:
                image.verify()
            with Image.open(png) as image:
                width, height = image.size
            png_ok &= width >= 2500 and height >= 1500
        passed = png_ok and pdf_ok
        all_pass &= passed
        rows.append({
            "figure": stem,
            "png": str(png),
            "png_size_bytes": png.stat().st_size if png.is_file() else -1,
            "width_px": width,
            "height_px": height,
            "png_sha256": sha256(png) if png.is_file() else "MISSING",
            "pdf": str(pdf),
            "pdf_size_bytes": pdf.stat().st_size if pdf.is_file() else -1,
            "pdf_sha256": sha256(pdf) if pdf.is_file() else "MISSING",
            "status": "PASS" if passed else "FAIL",
        })
    expected_names = sorted([f"{stem}.{suffix}" for stem in EXPECTED for suffix in ["png", "pdf"]])
    metadata_ok = metadata.get("status") == "completed" and sorted(metadata.get("files", [])) == expected_names
    all_pass &= metadata_ok
    pd.DataFrame(rows).to_csv(output_tsv, sep="\t", index=False)
    summary = {
        "status": "passed" if all_pass else "failed",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "memory_gate": gate_line,
        "numeric_verification": "35/35 PASS",
        "figures_passed": sum(row["status"] == "PASS" for row in rows),
        "figures_failed": sum(row["status"] == "FAIL" for row in rows),
        "metadata_ok": metadata_ok,
        "manual_review": str(QC / "joint_figure_manual_review_v1.tsv"),
    }
    output_json.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    if not all_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
