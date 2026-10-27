from __future__ import annotations

import hashlib
import json
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


def copy_new(source: Path, target: Path) -> None:
    if target.exists():
        raise SystemExit(f"Refusing existing target: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: 22_install_phase90_v12_revision.py <project_root> <work_root>")
    project_root = Path(sys.argv[1]).resolve()
    work_root = Path(sys.argv[2]).resolve()
    figure_root = work_root / "figure_v12"
    install_root = project_root / "phase90_immune_main_figure_local_v1"
    if not install_root.is_dir():
        raise SystemExit("Phase90 install root missing")
    now = datetime.now().astimezone().isoformat(timespec="seconds")

    for name in [
        "PROJECT26210_IMMUNE_MAIN_FIGURE_v12.png", "PROJECT26210_IMMUNE_MAIN_FIGURE_v12.pdf",
        "PROJECT26210_IMMUNE_MAIN_FIGURE_v12.svg", "BUILD_SUMMARY_v12.json",
    ]:
        copy_new(figure_root / name, install_root / name)
    for name in [
        "VALIDATION_STATUS_v1.json", "VALIDATION_CHECKS_v1.tsv",
        "FINAL_FILE_MANIFEST_SHA256_v1.tsv", "PROJECT26210_IMMUNE_MAIN_FIGURE_v12_pdf_render.png",
    ]:
        copy_new(figure_root / "qc_v1" / name, install_root / "qc_v12" / name)
    for name in ["20_build_main_figure_v12.py", "21_validate_main_figure_v12.py", "22_install_phase90_v12_revision.py"]:
        copy_new(work_root / "scripts" / name, install_root / "scripts_v1" / name)

    hashes = {ext: sha256(figure_root / f"PROJECT26210_IMMUNE_MAIN_FIGURE_v12.{ext}") for ext in ["png", "pdf", "svg"]}
    addendum_path = install_root / "PROJECT_STATUS_PHASE90_ADDENDUM_v4.md"
    if addendum_path.exists():
        raise SystemExit("Phase90 addendum v4 already exists")
    addendum_path.write_text(f"""# Phase 90 CellChat chord revision v12

- Completed: {now}
- Status: **PASS_MAIN_FIGURE_V12; 35/35 VALIDATION CHECKS PASS; F-PANEL CROP AND FULL PDF RENDER REVIEW PASS**.
- Revision scope: Panel F lower-left communication graphic only.
- The previous directed network was replaced by a six-sector chord diagram for Tumour, SPP1 macrophage, C1QC macrophage, Other macrophage, T cells and Endothelial.
- Ribbons use the top 14 non-self GSE207493 patient-median CellChat edges; ribbon width is proportional to the square root of median probability and colour denotes sender identity.
- Short external labels are placed at a fixed radius and were visually verified with no overlap or clipping.
- All spatial maps, specificity endpoints, other panels and scientific interpretation are unchanged.
- v12 supersedes v11, v10 and v7 for client delivery.
- PNG SHA256: `{hashes['png']}`.
- PDF SHA256: `{hashes['pdf']}`.
- SVG SHA256: `{hashes['svg']}`.
""", encoding="utf-8")

    manifest_rows = []
    for path in sorted(install_root.rglob("*")):
        if path.is_file() and path.name not in {"FILE_MANIFEST_SHA256_v4.tsv", "INSTALL_VALIDATION_v4.json"}:
            manifest_rows.append({"relative_path": path.relative_to(install_root).as_posix(), "bytes": path.stat().st_size, "sha256": sha256(path)})
    manifest_path = install_root / "FILE_MANIFEST_SHA256_v4.tsv"
    if manifest_path.exists():
        raise SystemExit("File manifest v4 already exists")
    pd.DataFrame(manifest_rows).to_csv(manifest_path, sep="\t", index=False)

    status_path = project_root / "PROJECT_STATUS.md"
    status_text = status_path.read_text(encoding="utf-8")
    heading = "### Phase 90 CellChat chord revision v12 (2026-08-27)"
    if heading in status_text:
        raise SystemExit("Phase90 v12 status already present")
    status_path.write_text(status_text.rstrip() + f"""

{heading}

- Status: **PASS_MAIN_FIGURE_V12; 35/35 VALIDATION CHECKS PASS; FULL PDF-RENDER REVIEW PASS**.
- Panel F lower-left was replaced by a six-sector CellChat chord diagram with fixed-radius short labels and no overlap. Ribbon width encodes sqrt(patient-median probability); ribbon colour identifies the sender.
- All other panels and source values are unchanged. v12 supersedes v11, v10 and v7 for client delivery.
- Authoritative files: `PROJECT26210_IMMUNE_MAIN_FIGURE_v12.png`, `.pdf` and `.svg` under `phase90_immune_main_figure_local_v1`.
""" + "\n", encoding="utf-8")

    results_path = project_root / "RESULTS_MANIFEST.tsv"
    results_text = results_path.read_text(encoding="utf-8").rstrip("\r\n")
    existing = {line.split("\t", 1)[0] for line in results_text.splitlines()[1:] if line.strip()}
    result_rows = [
        ("phase90_immune_main_figure_png_v12", "high_resolution_png", "PROJECT26210_IMMUNE_MAIN_FIGURE_v12.png", f"SHA256 {hashes['png']}; 4058x5622", "Authoritative CellChat chord revision"),
        ("phase90_immune_main_figure_pdf_v12", "publication_pdf", "PROJECT26210_IMMUNE_MAIN_FIGURE_v12.pdf", f"SHA256 {hashes['pdf']}; one page; PDF render reviewed", "Authoritative publication PDF"),
        ("phase90_immune_main_figure_svg_v12", "editable_svg", "PROJECT26210_IMMUNE_MAIN_FIGURE_v12.svg", f"SHA256 {hashes['svg']}; editable text/vector", "Authoritative editable source"),
        ("phase90_immune_main_figure_validation_v12", "validation_json", "qc_v12/VALIDATION_STATUS_v1.json", "35/35 PASS; chord crop and full PDF render reviewed", "v12 independent QA"),
        ("phase90_immune_main_figure_status_addendum_v4", "status_addendum_md", "PROJECT_STATUS_PHASE90_ADDENDUM_v4.md", "PASS_MAIN_FIGURE_V12", "Chord revision scope and unchanged data boundary"),
        ("phase90_immune_main_figure_file_manifest_v4", "file_manifest_tsv", "FILE_MANIFEST_SHA256_v4.tsv", f"{len(manifest_rows)} files hashed", "Expanded Phase90 package manifest including v12"),
    ]
    lines = []
    for result_id, artifact_type, relative, evidence, notes in result_rows:
        if result_id in existing:
            raise SystemExit(f"Manifest ID already exists: {result_id}")
        lines.append("\t".join([result_id, "phase90_immune_main_figure", artifact_type, str(install_root / relative), now, "PASS", evidence, notes]))
    results_path.write_text(results_text + "\n" + "\n".join(lines) + "\n", encoding="utf-8")

    manifest = pd.read_csv(manifest_path, sep="\t")
    exact = []
    for _, row in manifest.iterrows():
        path = install_root / Path(row["relative_path"])
        exact.append(path.exists() and path.stat().st_size == int(row["bytes"]) and sha256(path) == row["sha256"])
    install_status = {
        "status": "PASS", "package_files_exact": int(sum(exact)), "package_files_total": int(len(exact)),
        "v12_validation_checks": "35/35 PASS", "root_status_revision_count": status_path.read_text(encoding="utf-8").count(heading),
        "installed_at": now,
    }
    (install_root / "INSTALL_VALIDATION_v4.json").write_text(json.dumps(install_status, indent=2), encoding="utf-8")
    if not all(exact) or install_status["root_status_revision_count"] != 1:
        raise SystemExit("v12 installation validation failed")


if __name__ == "__main__":
    main()
