#!/usr/bin/env python3
"""Recover a valid first gzip member while preserving any trailing bytes."""

from __future__ import annotations

import argparse
import shutil
import zlib
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("repaired", type=Path)
    parser.add_argument("trailing", type=Path)
    args = parser.parse_args()

    source_size = args.source.stat().st_size
    decoder = zlib.decompressobj(wbits=31)
    bytes_fed = 0
    gzip_end = None
    with args.source.open("rb") as handle:
        while True:
            chunk = handle.read(8 * 1024 * 1024)
            if not chunk:
                break
            bytes_fed += len(chunk)
            decoder.decompress(chunk)
            if decoder.eof:
                gzip_end = bytes_fed - len(decoder.unused_data)
                break

    if gzip_end is None:
        raise RuntimeError("no complete gzip member found")
    trailing_size = source_size - gzip_end
    if trailing_size <= 0:
        raise RuntimeError("source has no trailing bytes; repair is unnecessary")

    args.repaired.parent.mkdir(parents=True, exist_ok=True)
    args.trailing.parent.mkdir(parents=True, exist_ok=True)
    with args.source.open("rb") as source, args.repaired.open("wb") as repaired:
        remaining = gzip_end
        while remaining:
            block = source.read(min(8 * 1024 * 1024, remaining))
            if not block:
                raise RuntimeError("unexpected EOF while copying valid gzip member")
            repaired.write(block)
            remaining -= len(block)
        with args.trailing.open("wb") as trailing:
            shutil.copyfileobj(source, trailing, length=8 * 1024 * 1024)

    print(f"source_bytes={source_size}")
    print(f"gzip_member_bytes={gzip_end}")
    print(f"trailing_bytes={trailing_size}")


if __name__ == "__main__":
    main()
