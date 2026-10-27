#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

BASE = Path("/mnt/BioDisk/wangmingyuan/project26210_immunotherapy_scrna_download_precheck_v1")
FILELIST = BASE / "metadata" / "GSE294109_filelist.txt"
FILES = BASE / "downloads" / "GSE294109" / "files"
STAGING = BASE / "downloads" / "GSE294109" / "range_staging_v1"
REPORT = BASE / "precheck" / "GSE294109_range_repair_manifest_v1.tsv"
CHUNK_SIZE = 16 * 1024 * 1024
MAX_WORKERS = 12


@dataclass(frozen=True)
class Job:
    filename: str
    url: str
    start: int
    end: int
    part: Path


def official_files() -> dict[str, int]:
    result: dict[str, int] = {}
    with FILELIST.open(encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if row and row[0] == "File":
                result[row[1]] = int(row[3])
    return result


def sample_url(filename: str) -> str:
    gsm = filename.split("_", 1)[0]
    group = f"{gsm[:-3]}nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{group}/{gsm}/suppl/{filename}"


def download(job: Job) -> tuple[str, int, int]:
    expected = job.end - job.start + 1
    if job.part.exists() and job.part.stat().st_size == expected:
        return job.filename, job.start, expected
    tmp = job.part.with_suffix(job.part.suffix + ".tmp")
    if tmp.exists() and tmp.stat().st_size > expected:
        tmp.unlink()
    for attempt in range(1, 51):
        current = tmp.stat().st_size if tmp.exists() else 0
        if current == expected:
            break
        request_start = job.start + current
        with tmp.open("ab") as target:
            result = subprocess.run(
                [
                    "curl",
                    "-fsSL",
                    "--connect-timeout",
                    "30",
                    "--max-time",
                    "180",
                    "--range",
                    f"{request_start}-{job.end}",
                    job.url,
                ],
                stdout=target,
            )
        if result.returncode != 0:
            time.sleep(min(5 * attempt, 30))
    else:
        raise RuntimeError(f"Range retry limit reached for {job.filename} {job.start}-{job.end}")
    actual = tmp.stat().st_size
    if actual != expected:
        raise RuntimeError(f"Range size mismatch for {job.filename}: expected {expected}, got {actual}")
    tmp.replace(job.part)
    return job.filename, job.start, actual


def main() -> None:
    FILES.mkdir(parents=True, exist_ok=True)
    STAGING.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    expected_sizes = official_files()
    initial_sizes = {name: (FILES / name).stat().st_size if (FILES / name).exists() else 0 for name in expected_sizes}

    jobs: list[Job] = []
    for filename, expected in expected_sizes.items():
        actual = initial_sizes[filename]
        if actual > expected:
            raise RuntimeError(f"Oversize file requires manual review: {filename} {actual}>{expected}")
        for start in range(actual, expected, CHUNK_SIZE):
            end = min(start + CHUNK_SIZE - 1, expected - 1)
            part = STAGING / f"{filename}.{start:012d}-{end:012d}.part"
            jobs.append(Job(filename, sample_url(filename), start, end, part))

    print(f"incomplete_files={sum(initial_sizes[n] != s for n, s in expected_sizes.items())}", flush=True)
    print(f"range_jobs={len(jobs)} workers={MAX_WORKERS}", flush=True)
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures = [pool.submit(download, job) for job in jobs]
        for index, future in enumerate(as_completed(futures), start=1):
            filename, start, size = future.result()
            print(f"completed={index}/{len(jobs)} file={filename} start={start} bytes={size}", flush=True)

    jobs_by_file: dict[str, list[Job]] = {}
    for job in jobs:
        jobs_by_file.setdefault(job.filename, []).append(job)

    for filename, expected in expected_sizes.items():
        original = FILES / filename
        actual = initial_sizes[filename]
        if actual == expected:
            continue
        assembled = STAGING / f"{filename}.assembled"
        with assembled.open("wb") as target:
            if actual:
                with original.open("rb") as source:
                    shutil.copyfileobj(source, target, length=8 * 1024 * 1024)
            for job in sorted(jobs_by_file[filename], key=lambda item: item.start):
                with job.part.open("rb") as source:
                    shutil.copyfileobj(source, target, length=8 * 1024 * 1024)
        if assembled.stat().st_size != expected:
            raise RuntimeError(f"Assembled size mismatch: {filename}")
        os.replace(assembled, original)

    with REPORT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["filename", "official_bytes", "initial_bytes", "final_bytes", "status"])
        for filename, expected in expected_sizes.items():
            final = (FILES / filename).stat().st_size if (FILES / filename).exists() else 0
            writer.writerow([filename, expected, initial_sizes[filename], final, "PASS" if final == expected else "FAIL"])

    failures = [name for name, expected in expected_sizes.items() if (FILES / name).stat().st_size != expected]
    if failures:
        raise RuntimeError(f"Final byte-size failures: {failures}")
    print(f"all_files_match=PASS file_count={len(expected_sizes)}", flush=True)


if __name__ == "__main__":
    main()
