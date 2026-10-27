#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from zipfile import ZipFile


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
REPORT = ROOT / "docs/client_reports/part4_spatial/project26210_TYMP_5MTA_SPP1机制空间验证精选版_v2.docx"
FULL_REPORT = ROOT / "docs/client_reports/part4_spatial/project26210_第四部分空间转录组阶段性报告_客户样式综合图版_POSTN_APOE论文对标增强版_v2.docx"
RESULT = ROOT / "results/part4_spatial/tymp_curated_validation_report_v2"
FIG1 = RESULT / "Figure1_TYMP_MTA_SPP1_spatial_niche_8panel_v2.png"
FIG2 = RESULT / "Figure2_TYMP_SPP1_patient_spatial_atlas_6panel_v2.png"
PDF = RESULT / "project26210_TYMP_5MTA_SPP1_mechanism_spatial_atlas_curated_report_v2.pdf"
OUT = RESULT / "installation_verification_v2.tsv"
SUMMARY = RESULT / "installation_verification_summary_v2.json"
W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
A = "http://schemas.openxmlformats.org/drawingml/2006/main"


def q(ns: str, tag: str) -> str:
    return f"{{{ns}}}{tag}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


checks: list[tuple[str, bool, str]] = []


def check(name: str, condition: bool, detail: str) -> None:
    checks.append((name, bool(condition), detail))


expected = {
    REPORT: "55b95f7a15b3f33ba9e3d10f65bb88a76d9104882bbc26585508790f5bc55f39",
    PDF: "0edf5100e08495e02ac152d12f6c05cd181d4e391117e7dad56999c809ad806c",
    FIG1: "7c70c5dfb4692e0efe5be1201be7f96849e78ec64bae78f82151e08ce890bf76",
    FIG2: "8d3180730ff099d94db0b5d328e99bc55506a9eb80ffbda852e731e31ec4794e",
    FULL_REPORT: "4c830304b52b7e7ad9c4e40902cc81558e5b50690e4089ac4ee52ffd0b4db69b",
}
for path, digest in expected.items():
    check(f"exists_{path.name}", path.is_file(), str(path))
    if path.is_file():
        check(f"hash_{path.name}", sha256(path) == digest, sha256(path))

with ZipFile(REPORT) as archive:
    check("report_zip_crc", archive.testzip() is None, str(archive.testzip()))
    media = [name for name in archive.namelist() if name.startswith("word/media/")]
    media_hashes = {hashlib.sha256(archive.read(name)).hexdigest() for name in media}
    check("report_two_media", len(media) == 2, str(media))
    check("report_embedded_master_figures", media_hashes == {expected[FIG1], expected[FIG2]}, str(media_hashes))
    document = ET.fromstring(archive.read("word/document.xml"))
    text = "".join(item.text or "" for item in document.iter(q(W, "t")))
    check("report_five_tables", len(list(document.iter(q(W, "tbl")))) == 5, str(len(list(document.iter(q(W, "tbl"))))))
    check("report_two_drawings", len(list(document.iter(q(W, "drawing")))) == 2, str(len(list(document.iter(q(W, "drawing"))))))
    check("report_required_text", all(marker in text for marker in ["执行摘要", "10/12", "推荐客户表述", "精选版交付说明"]), "required markers")

independent = json.loads((RESULT / "curated_report_independent_verification_summary_v2.json").read_text(encoding="utf-8"))
render_qc = json.loads((RESULT / "curated_report_render_qc_v2.json").read_text(encoding="utf-8"))
install = json.loads((RESULT / "installation_metadata_v2.json").read_text(encoding="utf-8"))
check("independent_58_of_58", independent.get("status") == "passed" and independent.get("passed") == independent.get("total") == 58, json.dumps(independent, ensure_ascii=False))
check("render_7_of_7", render_qc.get("status") == "passed" and render_qc.get("pdf", {}).get("pages") == 7, json.dumps(render_qc, ensure_ascii=False))
check("install_candidate_status", install.get("status") == "installed" and install.get("report_status") == "candidate", json.dumps(install, ensure_ascii=False))

status_text = (ROOT / "PROJECT_STATUS.md").read_text(encoding="utf-8")
handoff_text = (ROOT / "CODEX_HANDOFF.md").read_text(encoding="utf-8")
manifest_text = (ROOT / "RESULTS_MANIFEST.tsv").read_text(encoding="utf-8")
check("project_status_registered", str(REPORT) in status_text and "58/58" in status_text and "7/7" in status_text, "PROJECT_STATUS.md")
check("handoff_registered", str(REPORT) in handoff_text and str(RESULT) in handoff_text, "CODEX_HANDOFF.md")
for result_id in [
    "p4_tymp_curated_report_v2_20260803",
    "p4_tymp_curated_figure1_v2_20260803",
    "p4_tymp_curated_figure2_v2_20260803",
    "p4_tymp_curated_report_verification_v2_atlas_20260803",
    "p4_tymp_curated_report_manual_v2_atlas_20260803",
]:
    check(f"manifest_{result_id}", manifest_text.count(result_id) == 1, str(manifest_text.count(result_id)))
backup = RESULT / "project_records_before_install"
check("records_backup_complete", backup.is_dir() and all((backup / name).is_file() for name in ["PROJECT_STATUS.md", "CODEX_HANDOFF.md", "RESULTS_MANIFEST.tsv"]), str(backup))

passed = sum(ok for _, ok, _ in checks)
OUT.write_text(
    "check\tstatus\tdetail\n" + "".join(f"{name}\t{'PASS' if ok else 'FAIL'}\t{detail}\n" for name, ok, detail in checks),
    encoding="utf-8",
)
summary = {
    "status": "passed" if passed == len(checks) else "failed",
    "passed": passed,
    "total": len(checks),
    "failed": [name for name, ok, _ in checks if not ok],
    "report": str(REPORT),
    "report_sha256": sha256(REPORT),
    "full_report_unchanged_sha256": sha256(FULL_REPORT),
    "verification": str(OUT),
}
SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps(summary, ensure_ascii=False, indent=2))
if passed != len(checks):
    raise SystemExit(1)
