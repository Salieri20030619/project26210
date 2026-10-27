#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import pandas as pd
from PIL import Image


ROOT = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_three_figures_stage_v1")
FIGDIR = ROOT / "figures"
META = ROOT / "figure_generation_metadata_v4.json"
MANIFEST = ROOT / "panel_source_manifest_v4.tsv"
LEGENDS = ROOT / "figure_legends_v4.md"
REPORT = ROOT / "figure_qc_report_v4.tsv"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def pdf_pages(path: Path) -> int:
    text = subprocess.run(["pdfinfo", str(path)], check=True, text=True,
                          stdout=subprocess.PIPE).stdout
    for line in text.splitlines():
        if line.startswith("Pages:"):
            return int(line.split(":", 1)[1].strip())
    raise RuntimeError(f"No page count in pdfinfo output: {path}")


def main() -> None:
    meta = json.loads(META.read_text(encoding="utf-8"))
    rows: list[dict[str, str]] = []

    def check(name: str, passed: bool, evidence: str) -> None:
        rows.append({"check": name, "status": "PASS" if passed else "FAIL", "evidence": evidence})

    for key in ("figure1", "figure2", "figure3"):
        output = meta["outputs"][key]
        png, pdf = Path(output["png"]), Path(output["pdf"])
        check(f"{key}_files_exist", png.is_file() and pdf.is_file(), f"{png}; {pdf}")
        with Image.open(png) as image:
            dims = image.size
        check(f"{key}_png_dimensions", dims[0] >= 8_000 and dims[1] >= 4_500,
              f"{dims[0]}x{dims[1]}; 450 dpi export with tight bounding box")
        check(f"{key}_pdf_single_page", pdf_pages(pdf) == 1, f"pages={pdf_pages(pdf)}")
        check(f"{key}_png_hash", sha256(png) == output["png_sha256"], sha256(png))
        check(f"{key}_pdf_hash", sha256(pdf) == output["pdf_sha256"], sha256(pdf))
        check(f"{key}_file_size", png.stat().st_size > 1_000_000 and pdf.stat().st_size > 50_000,
              f"png={png.stat().st_size}; pdf={pdf.stat().st_size}")

    panel = pd.read_csv(MANIFEST, sep="\t", dtype=str)
    counts = panel.groupby("figure")["panel"].nunique().to_dict()
    check("panel_manifest_counts", counts == {"1": 5, "2": 6, "3": 6}, json.dumps(counts, sort_keys=True))
    check("panel_manifest_sources", panel["source"].notna().all() and (panel["source"].str.len() > 0).all(),
          f"rows={len(panel)}")

    legend_text = LEGENDS.read_text(encoding="utf-8")
    check("three_complete_legends", all(f"## Figure {i}." in legend_text for i in (1, 2, 3)),
          f"characters={len(legend_text)}")
    check("inference_boundaries", all(term in legend_text for term in
          ("do not independently prove", "not a pathology-defined", "No lineage direction")),
          "mechanistic, pathology-boundary, and lineage caveats present")

    summaries = meta["summaries"]
    check("figure1_statistics", summaries["figure1"]["glmf_tymp_mta_fdr_pass"] == 4
          and summaries["figure1"]["reference_light_positive_patients"] == 10
          and abs(summaries["figure1"]["reference_light_global_fdr"] - 0.154296875) < 1e-12,
          json.dumps(summaries["figure1"], sort_keys=True))
    check("figure2_statistics", summaries["figure2"]["representative_patient"] == "Z43_T"
          and summaries["figure2"]["supported_relations"] == 2,
          json.dumps(summaries["figure2"], sort_keys=True))
    check("figure3_statistics", summaries["figure3"]["cells"] == 53794
          and summaries["figure3"]["samples"] == 63,
          json.dumps(summaries["figure3"], sort_keys=True))

    result = pd.DataFrame(rows)
    result.to_csv(REPORT, sep="\t", index=False)
    if not result["status"].eq("PASS").all():
        raise SystemExit(result[result["status"].eq("FAIL")].to_string(index=False))
    print(json.dumps({"status": "PASS", "checks": len(result), "report": str(REPORT)}, indent=2))


if __name__ == "__main__":
    main()
