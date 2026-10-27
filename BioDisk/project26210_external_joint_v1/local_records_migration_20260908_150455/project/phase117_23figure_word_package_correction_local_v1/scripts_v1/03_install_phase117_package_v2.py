from __future__ import annotations

import csv
import hashlib
import json
import shutil
import zipfile
from datetime import datetime, timezone, timedelta
from pathlib import Path


WORK = Path(r"C:\Users\31937\project26210_phase117_23figure_word_package_work_v1")
STAGING = WORK / "package_v2"
ROOT = Path(r"D:\project26210_external_joint_v1")
PHASE = ROOT / "phase117_23figure_word_package_correction_local_v1"
DELIVERABLE = PHASE / "deliverable_v2"
DOWNLOAD = Path(r"C:\Users\31937\Downloads\项目26210_2.1-2.2图件证据包_23图完整版_20260902_v2")
ZIP_PATH = DOWNLOAD.with_suffix(".zip")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest().upper()


def files_under(root: Path) -> list[Path]:
    return sorted((p for p in root.rglob("*") if p.is_file()), key=lambda p: str(p.relative_to(root)).lower())


for required in [STAGING, WORK / "validation_v2", WORK / "render_word_v1"]:
    if not required.exists():
        raise FileNotFoundError(required)
for target in [PHASE, DOWNLOAD, ZIP_PATH]:
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite: {target}")

PHASE.mkdir(parents=True)
shutil.copytree(STAGING, DELIVERABLE)
shutil.copytree(WORK / "validation_v2", PHASE / "validation_v2")
shutil.copytree(WORK / "render_word_v1", PHASE / "word_render_v2")
(PHASE / "scripts_v1").mkdir()
for name in [
    "01_build_23figure_word_v2.py",
    "02_validate_23figure_package_v2.py",
    "03_install_phase117_package_v2.py",
]:
    shutil.copy2(WORK / name, PHASE / "scripts_v1" / name)

shutil.copytree(DELIVERABLE, DOWNLOAD)
with zipfile.ZipFile(ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
    for path in files_under(DOWNLOAD):
        zf.write(path, arcname=str(Path(DOWNLOAD.name) / path.relative_to(DOWNLOAD)))

auth = {str(p.relative_to(DELIVERABLE)): sha256(p) for p in files_under(DELIVERABLE)}
down = {str(p.relative_to(DOWNLOAD)): sha256(p) for p in files_under(DOWNLOAD)}
if auth != down:
    raise RuntimeError("Downloads copy does not exactly match authoritative deliverable")

png_count = sum(p.suffix.lower() == ".png" for p in files_under(DELIVERABLE / "01_图件PNG_23张"))
pdf_backup_count = sum(p.suffix.lower() == ".pdf" for p in files_under(DELIVERABLE / "02_矢量PDF备份_18张"))
forbidden = [p for p in files_under(DELIVERABLE) if p.suffix.lower() in {".tsv", ".json"}]
if png_count != 23 or pdf_backup_count != 18 or forbidden:
    raise RuntimeError(f"Unexpected package counts: PNG={png_count}, PDF backups={pdf_backup_count}, forbidden={forbidden}")

with zipfile.ZipFile(ZIP_PATH) as zf:
    zip_count = len([n for n in zf.namelist() if not n.endswith("/")])
if zip_count != len(auth):
    raise RuntimeError(f"ZIP entry mismatch: {zip_count} != {len(auth)}")

manifest = PHASE / "FILE_MANIFEST_SHA256_v2.tsv"
with manifest.open("w", encoding="utf-8-sig", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["relative_path", "bytes", "sha256"])
    for rel in sorted(auth):
        p = DELIVERABLE / rel
        writer.writerow([rel, p.stat().st_size, auth[rel]])

tz = timezone(timedelta(hours=8))
summary = {
    "status": "PASS",
    "generated_at": datetime.now(tz).isoformat(timespec="seconds"),
    "authoritative_phase": str(PHASE),
    "deliverable": str(DELIVERABLE),
    "downloads": str(DOWNLOAD),
    "zip": str(ZIP_PATH),
    "deliverable_file_count": len(auth),
    "downloads_exact_hash_matches": len(auth),
    "zip_file_entries": zip_count,
    "png_count": png_count,
    "pdf_backup_count": pdf_backup_count,
    "forbidden_tsv_json_count": len(forbidden),
    "docx_sha256": sha256(DELIVERABLE / "项目26210_2.1-2.2图件与统计说明_23图完整版_v2.docx"),
    "pdf_sha256": sha256(DELIVERABLE / "项目26210_2.1-2.2图件与统计说明_23图完整版_v2.pdf"),
    "zip_sha256": sha256(ZIP_PATH),
}
(PHASE / "INSTALL_AND_DELIVERY_VALIDATION_v2.json").write_text(
    json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
)
print(json.dumps(summary, ensure_ascii=False, indent=2))
