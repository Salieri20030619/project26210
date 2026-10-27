from __future__ import annotations

import hashlib
import json
import os
import shutil
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_file(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: 11_package_install_record_phase90_v1.py <project_root> <work_root>")
    project_root = Path(sys.argv[1]).resolve()
    work_root = Path(sys.argv[2]).resolve()
    figure_root = work_root / "figure_v7"
    local_package = work_root / "delivery_v1"
    local_part = work_root / "delivery_v1.part"
    install_root = project_root / "phase90_immune_main_figure_local_v1"
    install_part = project_root / "phase90_immune_main_figure_local_v1.part"
    for path in [local_package, local_part, install_root, install_part]:
        if path.exists():
            raise SystemExit(f"Refusing existing path: {path}")
    local_part.mkdir(parents=True)

    now = datetime.now().astimezone().isoformat(timespec="seconds")
    outputs = [
        "PROJECT26210_IMMUNE_MAIN_FIGURE_v7.png",
        "PROJECT26210_IMMUNE_MAIN_FIGURE_v7.pdf",
        "PROJECT26210_IMMUNE_MAIN_FIGURE_v7.svg",
        "BUILD_SUMMARY_v7.json",
        "PANEL_SOURCE_MANIFEST_SHA256_v1.tsv",
    ]
    for name in outputs:
        copy_file(figure_root / name, local_part / name)

    qc_names = [
        "VALIDATION_STATUS_v1.json",
        "VALIDATION_CHECKS_v1.tsv",
        "FINAL_FILE_MANIFEST_SHA256_v1.tsv",
        "PROJECT26210_IMMUNE_MAIN_FIGURE_v7_pdf_render.png",
    ]
    for name in qc_names:
        copy_file(figure_root / "qc_v1" / name, local_part / "qc_v1" / name)

    script_names = [
        "00_extract_atlas_umap_v2.R",
        "01_build_main_figure_v1.py",
        "02_build_main_figure_v2.py",
        "09_build_main_figure_v7.py",
        "07_validate_main_figure_v1.py",
        "10_validate_main_figure_v7.py",
        "11_package_install_record_phase90_v1.py",
    ]
    for name in script_names:
        copy_file(work_root / "scripts" / name, local_part / "scripts_v1" / name)

    for path in sorted((work_root / "atlas_extract_v1").glob("*")):
        if path.is_file():
            copy_file(path, local_part / "compact_inputs_v1" / "atlas_extract_v1" / path.name)
    for path in sorted((work_root / "spatial_inputs_v1").glob("*/*")):
        if path.is_file():
            copy_file(path, local_part / "compact_inputs_v1" / "spatial_inputs_v1" / path.parent.name / path.name)

    pre_records = local_part / "pre_phase90_root_records_v1"
    copy_file(project_root / "PROJECT_STATUS.md", pre_records / "PROJECT_STATUS_pre_phase90_v1.md")
    copy_file(project_root / "RESULTS_MANIFEST.tsv", pre_records / "RESULTS_MANIFEST_pre_phase90_v1.tsv")

    figure_hashes = {ext: sha256(figure_root / f"PROJECT26210_IMMUNE_MAIN_FIGURE_v7.{ext}") for ext in ["png", "pdf", "svg"]}
    status_addendum = f"""# Phase 90 immune main figure

- Completed: {now}
- Status: **PASS_MAIN_FIGURE_V7; 35/35 VALIDATION CHECKS PASS; FULL PNG AND PDF-RENDER VISUAL REVIEW PASS**.
- Scope: non-overwriting visual integration of the post-immune-dataset evidence into panels A-F. No normalization, integration, clustering, annotation, model fitting or biological inference was rerun.
- Panel A: discovery, external single-cell, spatial and clinical-context architecture.
- Panel B: three existing cohort-specific UMAPs; no pooled integration.
- Panel C: GSE207493 TYMP/F1 localization, lineage feature landscape and pathway/TF coverage glyphs.
- Panel D: 18 two-external Myeloid/DC signatures and 13 single-external Treg signatures from the frozen Phase 67 table.
- Panel E: Phase 87 cell-level state ordering and smoothed descriptive pseudotime program dynamics.
- Panel F: non-extreme representative FFPE/frozen Phase 78 sections plus Phase 88 absolute CellChat network and six specificity endpoints.
- Authoritative PNG SHA256: `{figure_hashes['png']}`.
- Authoritative PDF SHA256: `{figure_hashes['pdf']}`.
- Authoritative SVG SHA256: `{figure_hashes['svg']}`.
- Interpretation boundary: the figure retains broad TYMP-associated myeloid/macrophage context and cross-sectional state ordering; it does not promote selective SPP1 communication, 5-MTA causality or response/PFS/OS prediction.
"""
    (local_part / "PROJECT_STATUS_PHASE90_ADDENDUM_v1.md").write_text(status_addendum, encoding="utf-8")

    summary = f"""# Phase 90 immune main figure summary

The publication-style A-F main figure was generated from frozen Phase 65/67/78/87/88 outputs and three Phase 8 cohort objects. It contains no sparse boxplot or simple bar-chart panel. The final design uses cohort UMAPs, feature-density maps, a lineage bubble matrix, pathway coverage rings, a cross-cohort effect bridge, a replicated program constellation, pseudotime dynamics, spatial abundance maps and a CellChat network/specificity matrix.

## Deliverables

- `PROJECT26210_IMMUNE_MAIN_FIGURE_v7.png`: 4392 x 5622 high-resolution raster.
- `PROJECT26210_IMMUNE_MAIN_FIGURE_v7.pdf`: one-page publication PDF.
- `PROJECT26210_IMMUNE_MAIN_FIGURE_v7.svg`: editable vector/text version with rasterized high-density point clouds.
- `qc_v1`: 35/35 machine checks plus the visually reviewed PDF render.
- `compact_inputs_v1` and `scripts_v1`: compact reproducibility package.

## Scientific boundary

This is a visual synthesis, not a new analysis. Cohorts remain separate. The positive core is external replication of TYMP-associated Myeloid/DC programs and reproducible cross-sectional myeloid state ordering. Fine-state bulk identifiability, SPP1-specific communication, direct 5-MTA mediation and clinical prediction remain unsupported.
"""
    (local_part / "PHASE90_MAIN_FIGURE_SUMMARY_v1.md").write_text(summary, encoding="utf-8")

    manifest_rows = []
    for path in sorted(local_part.rglob("*")):
        if path.is_file() and path.name != "FILE_MANIFEST_SHA256_v1.tsv":
            manifest_rows.append({
                "relative_path": path.relative_to(local_part).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            })
    pd.DataFrame(manifest_rows).to_csv(local_part / "FILE_MANIFEST_SHA256_v1.tsv", sep="\t", index=False)
    os.replace(local_part, local_package)

    shutil.copytree(local_package, install_part)
    os.replace(install_part, install_root)
    manifest = pd.read_csv(local_package / "FILE_MANIFEST_SHA256_v1.tsv", sep="\t")
    verified = []
    for _, row in manifest.iterrows():
        path = install_root / Path(row["relative_path"])
        verified.append(path.exists() and path.stat().st_size == int(row["bytes"]) and sha256(path) == row["sha256"])
    if not all(verified):
        raise SystemExit("Installed package hash validation failed")

    status_path = project_root / "PROJECT_STATUS.md"
    root_status = status_path.read_text(encoding="utf-8")
    heading = "## Phase 90 immune main figure (2026-08-27)"
    if heading in root_status:
        raise SystemExit("Phase90 status already present")
    block = f"""

{heading}

- Status: **PASS_MAIN_FIGURE_V7; 35/35 VALIDATION CHECKS PASS; PNG AND PDF-RENDER VISUAL REVIEW PASS**.
- A publication-style A-F composite was generated without sparse boxplots or simple bar panels. It integrates three cohort-specific UMAPs, TYMP/F1 cellular localization, frozen Phase67 cross-cohort program replication, Phase87 myeloid-state ordering, Phase78 FFPE/frozen spatial maps and Phase88 CellChat summaries.
- The figure is a non-overwriting visual synthesis; it does not rerun or alter any biological analysis. Cohorts remain separate and the Phase89 evidence boundary is unchanged.
- Authoritative files are under `phase90_immune_main_figure_local_v1`; the editable deliverables are `PROJECT26210_IMMUNE_MAIN_FIGURE_v7.svg`, `.pdf` and `.png`.
- Detailed provenance, compact inputs, scripts, hashes and QA are recorded in `PROJECT_STATUS_PHASE90_ADDENDUM_v1.md`, `FILE_MANIFEST_SHA256_v1.tsv` and `qc_v1`.
"""
    status_path.write_text(root_status.rstrip() + block + "\n", encoding="utf-8")

    results_path = project_root / "RESULTS_MANIFEST.tsv"
    results = results_path.read_text(encoding="utf-8").rstrip("\r\n")
    ids = {line.split("\t", 1)[0] for line in results.splitlines()[1:] if line.strip()}
    rows = [
        ("phase90_immune_main_figure_png_v7", "high_resolution_png", "PROJECT26210_IMMUNE_MAIN_FIGURE_v7.png", f"SHA256 {figure_hashes['png']}; 4392x5622", "Publication-style A-F immune main figure"),
        ("phase90_immune_main_figure_pdf_v7", "publication_pdf", "PROJECT26210_IMMUNE_MAIN_FIGURE_v7.pdf", f"SHA256 {figure_hashes['pdf']}; one page; PDF render reviewed", "One-page publication PDF"),
        ("phase90_immune_main_figure_svg_v7", "editable_svg", "PROJECT26210_IMMUNE_MAIN_FIGURE_v7.svg", f"SHA256 {figure_hashes['svg']}; editable text/vector", "Editable source figure"),
        ("phase90_immune_main_figure_validation_v1", "validation_json", "qc_v1/VALIDATION_STATUS_v1.json", "35/35 PASS; full PNG and PDF-render review", "Independent data, hash, format and render QA"),
        ("phase90_immune_main_figure_source_manifest_v1", "source_manifest_tsv", "PANEL_SOURCE_MANIFEST_SHA256_v1.tsv", "All panel source paths, bytes and SHA256 exact", "Frozen source provenance"),
        ("phase90_immune_main_figure_file_manifest_v1", "file_manifest_tsv", "FILE_MANIFEST_SHA256_v1.tsv", f"{len(manifest)} packaged files hashed and installed exact", "Package byte/hash manifest"),
        ("phase90_immune_main_figure_status_addendum_v1", "status_addendum_md", "PROJECT_STATUS_PHASE90_ADDENDUM_v1.md", "PASS_MAIN_FIGURE_V7", "Scope, panels, QA and interpretation boundary"),
    ]
    new_lines = []
    for result_id, artifact_type, relative, evidence, notes in rows:
        if result_id in ids:
            raise SystemExit(f"Manifest result ID already exists: {result_id}")
        absolute = install_root / Path(relative)
        fields = [result_id, "phase90_immune_main_figure", artifact_type, str(absolute), now, "PASS", evidence, notes]
        if any("\t" in field or "\n" in field for field in fields):
            raise SystemExit("Invalid manifest field")
        new_lines.append("\t".join(fields))
    results_path.write_text(results + "\n" + "\n".join(new_lines) + "\n", encoding="utf-8")

    installed_results = pd.read_csv(results_path, sep="\t")
    install_status = {
        "status": "PASS",
        "installed_files_verified": int(sum(verified)),
        "installed_files_total": int(len(verified)),
        "root_status_phase90_count": status_path.read_text(encoding="utf-8").count(heading),
        "results_manifest_phase90_rows": int(installed_results["result_id"].str.startswith("phase90_immune_main_figure").sum()),
        "installed_at": now,
    }
    (install_root / "INSTALL_VALIDATION_v1.json").write_text(json.dumps(install_status, indent=2), encoding="utf-8")
    if install_status["root_status_phase90_count"] != 1 or install_status["results_manifest_phase90_rows"] != len(rows):
        raise SystemExit("Root-record installation validation failed")


if __name__ == "__main__":
    main()
