#!/usr/bin/env python3
from __future__ import annotations

import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

BASE = Path("/mnt/BioDisk/wangmingyuan/project26210_immunotherapy_scrna_download_precheck_v1")
STAGING = BASE / "downloads" / "GSE294109" / "range_staging_v1"
NAME = "GSM8898823_sc48T.matrix.mtx.gz"
RANGE_START = 4_999_025
RANGE_END = 21_776_240
PART = STAGING / f"{NAME}.{RANGE_START:012d}-{RANGE_END:012d}.part"
TMP = PART.with_suffix(PART.suffix + ".tmp")
URL = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8898nnn/GSM8898823/suppl/{NAME}"


def fetch(item: tuple[int, int, Path]) -> Path:
    start, end, output = item
    expected = end - start + 1
    subprocess.run(
        [
            "curl",
            "-fsSL",
            "--retry",
            "30",
            "--retry-all-errors",
            "--retry-delay",
            "5",
            "--connect-timeout",
            "30",
            "--max-time",
            "240",
            "--range",
            f"{start}-{end}",
            "--output",
            str(output),
            URL,
        ],
        check=True,
    )
    if output.stat().st_size != expected:
        raise RuntimeError(f"Tail range mismatch: {output}")
    return output


def main() -> None:
    expected_part = RANGE_END - RANGE_START + 1
    current = TMP.stat().st_size
    if current >= expected_part:
        if current == expected_part:
            TMP.replace(PART)
            return
        raise RuntimeError("Temporary part is oversize")

    remaining_start = RANGE_START + current
    remaining = RANGE_END - remaining_start + 1
    workers = 8
    unit = (remaining + workers - 1) // workers
    jobs: list[tuple[int, int, Path]] = []
    for index, start in enumerate(range(remaining_start, RANGE_END + 1, unit)):
        end = min(start + unit - 1, RANGE_END)
        output = STAGING / f"{NAME}.tail{index:02d}.{start:012d}-{end:012d}"
        jobs.append((start, end, output))

    print(f"current_bytes={current} remaining_bytes={remaining} tail_jobs={len(jobs)}", flush=True)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        outputs = list(pool.map(fetch, jobs))
    with TMP.open("ab") as target:
        for output in outputs:
            with output.open("rb") as source:
                target.write(source.read())
    if TMP.stat().st_size != expected_part:
        raise RuntimeError("Completed last range has wrong size")
    TMP.replace(PART)
    print(f"last_range_complete=PASS bytes={expected_part}", flush=True)


if __name__ == "__main__":
    main()
