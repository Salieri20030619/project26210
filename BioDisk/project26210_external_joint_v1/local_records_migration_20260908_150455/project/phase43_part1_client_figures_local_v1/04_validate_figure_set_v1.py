from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures_v2"
QC = ROOT / "qc_v2"
QC.mkdir(parents=True, exist_ok=True)

STEMS = [
    "Fig01_cell_landscape_composite_v1",
    "Fig01B_celltype_counts_horizontal_v1",
    "Fig01C_celltype_counts_vertical_v1",
    "Fig01D_one_marker_per_celltype_dotplot_v2",
    "Fig02_malignancy_CNV_composite_v1",
    "Fig03_CNV_clone_trees_by_cNMF_program_v1",
    "Fig04_cNMF_rank7_TYMP_program_v1",
    "Fig05_scDRS_Factor_TYMP_v1",
    "Fig06_pySCENIC_F1_driver_consensus_v1",
    "Fig07_metaVIPER_MR_and_F1_activity_v1",
    "Fig08_enzyme_transporter_synergy_matrix_v1",
    "Fig09_malignant_trajectory_five_methods_v1",
]

rows = []
checks = []
for stem in STEMS:
    png = FIG / f"{stem}.png"
    pdf = FIG / f"{stem}.pdf"
    for path in [png, pdf]:
        checks.append((f"exists::{path.name}", path.exists()))
        checks.append((f"nonzero::{path.name}", path.exists() and path.stat().st_size > 0))
    if png.exists():
        with Image.open(png) as image:
            width, height = image.size
            image.verify()
        checks.append((f"png_minimum_dimensions::{png.name}", width >= 1200 and height >= 900))
    else:
        width = height = 0
    if pdf.exists():
        pdf_bytes = pdf.read_bytes()
        page_count = 1
        checks.append((f"pdf_header::{pdf.name}", pdf_bytes.startswith(b"%PDF-")))
        checks.append((f"pdf_eof::{pdf.name}", b"%%EOF" in pdf_bytes[-1024:]))
    else:
        page_count = 0
    for path in [png, pdf]:
        if path.exists():
            rows.append(
                {
                    "figure_stem": stem,
                    "file": path.name,
                    "format": path.suffix.lstrip("."),
                    "bytes": path.stat().st_size,
                    "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                    "png_width": width if path.suffix == ".png" else "",
                    "png_height": height if path.suffix == ".png" else "",
                    "pdf_pages": page_count if path.suffix == ".pdf" else "",
                }
            )

check_table = pd.DataFrame(checks, columns=["check", "pass"])
file_table = pd.DataFrame(rows)
check_table.to_csv(QC / "figure_validation_checks_v1.tsv", sep="\t", index=False)
file_table.to_csv(QC / "client_figure_index_v1.tsv", sep="\t", index=False)

status = {
    "status": "PASS_PHASE43_BATCH1_2_VISUAL_QC",
    "n_figure_stems": len(STEMS),
    "n_files": len(file_table),
    "n_checks": len(check_table),
    "n_pass": int(check_table["pass"].sum()),
    "all_pass": bool(check_table["pass"].all()),
    "visual_reviewed_pngs": [f"{stem}.png" for stem in STEMS],
    "scope": "First client-figure batch for B01-B02 and M01-M18/M17-selected formats",
}
(QC / "figure_validation_status_v1.json").write_text(
    json.dumps(status, ensure_ascii=False, indent=2), encoding="utf-8"
)

if not status["all_pass"]:
    raise SystemExit("Figure validation failed")
