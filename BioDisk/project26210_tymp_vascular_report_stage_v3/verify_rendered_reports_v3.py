#!/usr/bin/env python3
"""Verify rendered PDFs and the recorded manual review for Part 4 v3 reports."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import subprocess
from datetime import datetime
from pathlib import Path


STAGE = Path("/mnt/BioDisk/wangmingyuan/project26210_tymp_vascular_report_stage_v3")
FULL_PDF = STAGE / "rendered_full_v3" / "project26210_第四部分空间转录组阶段性报告_客户样式综合图版_POSTN_APOE论文对标增强版_v3_candidate.pdf"
CURATED_PDF = STAGE / "rendered_curated_v3" / "project26210_TYMP_5MTA_SPP1机制空间验证精选版_v3_candidate.pdf"
MANUAL = STAGE / "report_render_manual_review_v3.tsv"
OUT_TSV = STAGE / "render_verification_v3.tsv"
OUT_JSON = STAGE / "render_verification_summary_v3.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def command(*args: str) -> str:
    return subprocess.run(args, check=True, text=True, capture_output=True).stdout


checks: list[dict[str, str]] = []


def add(name: str, passed: bool, observed: str, expected: str) -> None:
    checks.append(
        {
            "check": name,
            "status": "PASS" if passed else "FAIL",
            "observed": observed,
            "expected": expected,
        }
    )


for label, path, pages, page_size in (
    ("full", FULL_PDF, 76, "612 x 792 pts (letter)"),
    ("curated", CURATED_PDF, 11, "595.304 x 841.89 pts (A4)"),
):
    add(f"{label}_pdf_exists", path.is_file(), str(path), "file exists")
    info = command("pdfinfo", str(path))
    page_match = re.search(r"^Pages:\s+(\d+)$", info, flags=re.MULTILINE)
    size_match = re.search(r"^Page size:\s+(.+)$", info, flags=re.MULTILINE)
    observed_pages = int(page_match.group(1)) if page_match else -1
    observed_size = size_match.group(1).strip() if size_match else "missing"
    add(f"{label}_page_count", observed_pages == pages, str(observed_pages), str(pages))
    add(f"{label}_page_size", observed_size == page_size, observed_size, page_size)
    add(f"{label}_pdf_nonempty", path.stat().st_size > 100_000, str(path.stat().st_size), ">100000 bytes")

full_text = command("pdftotext", "-layout", str(FULL_PDF), "-")
curated_text = command("pdftotext", "-layout", str(CURATED_PDF), "-")

for marker in (
    "4.15 【血管生态位与髓系横断面状态排序扩展】",
    "图 4.50",
    "图 4.51",
    "图 4.52",
    "图 4.53",
    "35/35 独立检查",
    "DPT 不是 RNA velocity 或谱系示踪",
):
    add(f"full_marker_{marker}", marker in full_text, str(marker in full_text), "True")

for marker in (
    "3.2 血管生态位与髓系横断面状态排序",
    "图 3",
    "图 4",
    "横断面 DPT",
    "精选版不替代完整技术报告",
):
    add(f"curated_marker_{marker}", marker in curated_text, str(marker in curated_text), "True")

page_markers = (
    ("full", FULL_PDF, 68, "4.15.1 TYMP–SPP1 血管生态位的患者级空间关系"),
    ("full", FULL_PDF, 70, "4.15.2 髓系横断面状态排序"),
    ("full", FULL_PDF, 72, "4.15.3 QA 与证据边界"),
    ("curated", CURATED_PDF, 9, "3.2 血管生态位与髓系横断面状态排序"),
    ("curated", CURATED_PDF, 10, "4. 面向客户方向的结论"),
    ("curated", CURATED_PDF, 11, "5. 精选版交付说明"),
)
for label, path, page, marker in page_markers:
    text = command("pdftotext", "-f", str(page), "-l", str(page), "-layout", str(path), "-")
    add(f"{label}_page_{page}_marker", marker in text, str(marker in text), marker)

with MANUAL.open(encoding="utf-8", newline="") as handle:
    manual_rows = list(csv.DictReader(handle, delimiter="\t"))
add("manual_review_row_count", len(manual_rows) == 11, str(len(manual_rows)), "11")
add(
    "manual_review_all_pass",
    bool(manual_rows) and all(row["status"] == "PASS" for row in manual_rows),
    ",".join(row["status"] for row in manual_rows),
    "all PASS",
)
add(
    "manual_review_both_reports",
    {row["report"] for row in manual_rows} == {"full_v3", "curated_v3"},
    ",".join(sorted({row["report"] for row in manual_rows})),
    "curated_v3,full_v3",
)

with OUT_TSV.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["check", "status", "observed", "expected"], delimiter="\t")
    writer.writeheader()
    writer.writerows(checks)

failed = [row for row in checks if row["status"] != "PASS"]
summary = {
    "timestamp": datetime.now().astimezone().isoformat(timespec="seconds"),
    "status": "passed" if not failed else "failed",
    "checks_total": len(checks),
    "checks_passed": len(checks) - len(failed),
    "checks_failed": len(failed),
    "full_pdf": str(FULL_PDF),
    "full_pdf_sha256": sha256(FULL_PDF),
    "full_pdf_bytes": FULL_PDF.stat().st_size,
    "curated_pdf": str(CURATED_PDF),
    "curated_pdf_sha256": sha256(CURATED_PDF),
    "curated_pdf_bytes": CURATED_PDF.stat().st_size,
    "manual_review": str(MANUAL),
    "failed_checks": failed,
}
OUT_JSON.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps(summary, ensure_ascii=False, indent=2))
raise SystemExit(1 if failed else 0)
