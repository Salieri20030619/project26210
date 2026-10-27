from __future__ import annotations

import hashlib
from pathlib import Path

from PIL import Image


PHASE = Path(__file__).resolve().parent
FIG = PHASE / "figures_complete_v2"
QC = PHASE / "qc_v1"
QC.mkdir(parents=True, exist_ok=True)

expected_stems = [
    "Fig01B_celltype_counts_horizontal_v1",
    "Fig01C_celltype_counts_vertical_v1",
    "Fig01D_canonical_marker_dotplot_v1",
    "Fig01D_one_marker_per_celltype_dotplot_v2",
    "Fig01_cell_landscape_composite_v1",
    "Fig02B_CNV_burden_violin_v1",
    "Fig02_malignancy_CNV_composite_v1",
    "Fig03_CNV_clone_trees_by_cNMF_program_v1",
    "Fig04_cNMF_rank7_TYMP_program_v1",
    "Fig05_scDRS_Factor_TYMP_v1",
    "Fig06_pySCENIC_F1_driver_consensus_v1",
    "Fig07_metaVIPER_MR_and_F1_activity_v1",
    "Fig08_enzyme_transporter_synergy_matrix_v1",
    "Fig09_malignant_trajectory_five_methods_v1",
    "Fig10_pySCENIC_selected_regulon_networks_v1",
    "Fig11_metaVIPER_candidate_activity_umaps_v1",
    "Fig12_metabolic_pathway_activity_and_contrasts_v1",
    "Fig13_CellChat_communication_landscape_v1",
    "Fig14_NicheNet_ligand_receptor_evidence_v1",
    "Fig15_TCGA_TYMP_PRTG_IGDCC4_survival_v1",
    "FigS01_TCGA_TYMP_receptor_survival_all_axes_v1",
]

rows = []
for stem in expected_stems:
    for ext in ("png", "pdf"):
        path = FIG / f"{stem}.{ext}"
        exists = path.exists()
        rows.append((stem, ext, "exists", exists, str(path)))
        if not exists:
            continue
        size = path.stat().st_size
        rows.append((stem, ext, "nonempty", size > 1000, str(size)))
        if ext == "png":
            with Image.open(path) as im:
                w, h = im.size
                rows.append((stem, ext, "png_dimensions", w >= 1200 and h >= 700, f"{w}x{h}"))
                rows.append((stem, ext, "png_mode", im.mode in {"RGB", "RGBA"}, im.mode))
        else:
            raw = path.read_bytes()
            rows.append((stem, ext, "pdf_header", raw.startswith(b"%PDF-"), raw[:8].decode("latin1", errors="replace")))
            rows.append((stem, ext, "pdf_eof", b"%%EOF" in raw[-2048:], "EOF_in_tail"))
        sha = hashlib.sha256(path.read_bytes()).hexdigest()
        rows.append((stem, ext, "sha256", True, sha))

out = QC / "FIGURE_VALIDATION_CHECKS_v2.tsv"
with out.open("w", encoding="utf-8", newline="\n") as handle:
    handle.write("figure_stem\textension\tcheck\tpass\tdetail\n")
    for row in rows:
        handle.write("\t".join(map(str, row)) + "\n")

failures = [r for r in rows if r[2] != "sha256" and not r[3]]
summary = QC / "FIGURE_VALIDATION_SUMMARY_v2.txt"
summary.write_text(
    f"expected_stems\t{len(expected_stems)}\n"
    f"expected_files\t{len(expected_stems) * 2}\n"
    f"checks_total\t{len(rows)}\n"
    f"checks_failed\t{len(failures)}\n"
    f"status\t{'PASS' if not failures else 'FAIL'}\n",
    encoding="utf-8",
)
print(summary.read_text(encoding="utf-8"), end="")
if failures:
    raise SystemExit(1)
