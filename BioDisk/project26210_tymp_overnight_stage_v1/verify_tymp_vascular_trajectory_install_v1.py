#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import pandas as pd


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
STAGE_V = Path("/mnt/BioDisk/wangmingyuan/project26210_tymp_vascular_niche_stage_v1")
STAGE_T = Path("/mnt/BioDisk/wangmingyuan/project26210_tymp_overnight_stage_v1")
DEST_V = ROOT / "results/part4_spatial/tymp_vascular_niche_v1"
DEST_T = ROOT / "results/part4_spatial/tymp_myeloid_state_ordering_v1"
DEST_S = ROOT / "results/part4_spatial/tymp_vascular_trajectory_synthesis_v1"
OUT = STAGE_T / "installation_verification_v1.tsv"
SUMMARY = STAGE_T / "installation_verification_summary_v1.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if OUT.exists() or SUMMARY.exists():
        raise SystemExit("Refusing to overwrite installation verification")
    gate = subprocess.run([str(ROOT / "check_memory_gate.sh")], text=True, capture_output=True, check=False)
    gate_line = (gate.stdout or gate.stderr).strip()
    if gate.returncode != 0:
        raise RuntimeError(gate_line)
    rows = []

    def add(check: str, passed: bool, observed: object, expected: object) -> None:
        rows.append({"check": check, "status": "PASS" if passed else "FAIL",
                     "observed": str(observed), "expected": str(expected)})

    for label, source, destination in [
        ("vascular source tree", STAGE_V / "output_v1", DEST_V),
        ("trajectory source tree", STAGE_T / "trajectory_v1", DEST_T),
        ("formal figure tree", STAGE_T / "joint_figures_v4", DEST_S / "figures"),
        ("formal QC tree", STAGE_T / "joint_qc_v2", DEST_S / "qc"),
    ]:
        source_files = sorted(path.relative_to(source) for path in source.rglob("*") if path.is_file())
        matches = []
        for relative in source_files:
            target = destination / relative
            matches.append(target.is_file() and target.stat().st_size == (source / relative).stat().st_size
                           and sha256(target) == sha256(source / relative))
        add(label, bool(source_files) and all(matches), f"{sum(matches)}/{len(matches)}", f"{len(matches)}/{len(matches)}")

    installed_manifest = pd.read_csv(DEST_S / "installed_file_manifest_v1.tsv", sep="\t")
    manifest_matches = []
    for row in installed_manifest.itertuples(index=False):
        path = Path(row.path)
        manifest_matches.append(path.is_file() and path.stat().st_size == int(row.size_bytes) and sha256(path) == row.sha256)
    add("installed file manifest", all(manifest_matches), f"{sum(manifest_matches)}/{len(manifest_matches)}",
        f"{len(manifest_matches)}/{len(manifest_matches)}")

    result_manifest = ROOT / "RESULTS_MANIFEST.tsv"
    with result_manifest.open(newline="", encoding="utf-8") as handle:
        parsed = list(csv.reader(handle, delimiter="\t"))
    widths = {len(row) for row in parsed}
    ids = [row[0] for row in parsed[1:]]
    new_ids = [value for value in ids if value.endswith("_20260804") and value.startswith("p4_")]
    add("RESULTS_MANIFEST column count", widths == {18}, widths, {18})
    add("RESULTS_MANIFEST unique IDs", len(ids) == len(set(ids)), f"{len(ids)}/{len(set(ids))}", "equal")
    add("new Part4 manifest rows", len(new_ids) == 10, len(new_ids), 10)

    status = (ROOT / "PROJECT_STATUS.md").read_text(encoding="utf-8")
    handoff = (ROOT / "CODEX_HANDOFF.md").read_text(encoding="utf-8")
    add("PROJECT_STATUS synchronized", "Part 4 TYMP血管生态位与髓系横断面状态排序已完成" in status and str(DEST_S.relative_to(ROOT)) in status,
        "found" if "Part 4 TYMP血管生态位与髓系横断面状态排序已完成" in status else "missing", "found")
    add("CODEX_HANDOFF synchronized", "2026-08-04 Part 4 TYMP血管生态位与髓系状态排序终态" in handoff and str(DEST_S) in handoff,
        "found" if "2026-08-04 Part 4 TYMP血管生态位与髓系状态排序终态" in handoff else "missing", "found")

    reports = {
        ROOT / "docs/client_reports/part4_spatial/project26210_第四部分空间转录组阶段性报告_客户样式综合图版_POSTN_APOE论文对标增强版_v2.docx":
            "4c830304b52b7e7ad9c4e40902cc81558e5b50690e4089ac4ee52ffd0b4db69b",
        ROOT / "docs/client_reports/part4_spatial/project26210_TYMP_5MTA_SPP1机制空间验证精选版_v1.docx":
            "aa274dac1601d0c0d263de1f32e1a5c296f5ad3980e4eb04b7a0987d2a3ae230",
        ROOT / "docs/client_reports/part4_spatial/project26210_TYMP_5MTA_SPP1机制空间验证精选版_v2.docx":
            "55b95f7a15b3f33ba9e3d10f65bb88a76d9104882bbc26585508790f5bc55f39",
    }
    report_matches = [path.is_file() and sha256(path) == expected for path, expected in reports.items()]
    add("client DOCX hashes unchanged", all(report_matches), f"{sum(report_matches)}/3", "3/3")

    frame = pd.DataFrame(rows)
    frame.to_csv(OUT, sep="\t", index=False)
    passed = int(frame["status"].eq("PASS").sum())
    failed = int(frame["status"].eq("FAIL").sum())
    summary = {
        "status": "passed" if failed == 0 else "failed",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "memory_gate": gate_line,
        "checks_passed": passed,
        "checks_failed": failed,
        "client_docx_unchanged": all(report_matches),
    }
    SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(frame.to_string(index=False))
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
