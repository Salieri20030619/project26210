#!/usr/bin/env python3
"""Safely inventory and extract the official GSE207493 scRNA archive."""

from __future__ import annotations

import argparse
import csv
import os
import tarfile
import zipfile
from pathlib import Path


def inside(child: Path, parent: Path) -> bool:
    try:
        child.relative_to(parent)
        return True
    except ValueError:
        return False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("target", type=Path)
    parser.add_argument("inventory", type=Path)
    args = parser.parse_args()

    archive = args.archive.resolve(strict=True)
    target = args.target.resolve()
    target.mkdir(parents=True, exist_ok=True)
    args.inventory.parent.mkdir(parents=True, exist_ok=True)

    rows: list[tuple[str, str, str, int]] = []
    with tarfile.open(archive, "r:gz") as tf:
        members = tf.getmembers()
        for member in members:
            destination = (target / member.name).resolve()
            if not inside(destination, target):
                raise RuntimeError(f"unsafe archive path: {member.name!r}")
            if member.issym() or member.islnk():
                raise RuntimeError(f"links are not allowed: {member.name!r}")
            if not (member.isfile() or member.isdir()):
                raise RuntimeError(f"unsupported archive entry: {member.name!r}")
            kind = "file" if member.isfile() else "directory"
            rows.append(("outer_tar", member.name, kind, member.size))

        for member in members:
            tf.extract(member, target)

    for nested_zip in sorted(target.rglob("*.zip")):
        nested_target = target / nested_zip.stem
        nested_target.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(nested_zip) as zf:
            bad_member = zf.testzip()
            if bad_member is not None:
                raise RuntimeError(f"nested zip CRC failure: {bad_member!r}")
            infos = zf.infolist()
            for info in infos:
                destination = (nested_target / info.filename).resolve()
                if not inside(destination, nested_target):
                    raise RuntimeError(f"unsafe nested zip path: {info.filename!r}")
                unix_mode = (info.external_attr >> 16) & 0o170000
                if unix_mode == 0o120000:
                    raise RuntimeError(f"nested zip symlink is not allowed: {info.filename!r}")
                if info.flag_bits & 0x1:
                    raise RuntimeError(f"encrypted nested zip entry is not allowed: {info.filename!r}")
                kind = "directory" if info.is_dir() else "file"
                rows.append((nested_zip.name, info.filename, kind, info.file_size))
            zf.extractall(nested_target)

    with args.inventory.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["source", "member", "type", "size_bytes"])
        writer.writerows(rows)

    marker = target / "EXTRACTION_COMPLETE.txt"
    marker.write_text(
        f"archive={archive}\nentries={len(rows)}\nfiles={sum(r[2] == 'file' for r in rows)}\n",
        encoding="utf-8",
    )
    print(f"extracted {len(rows)} entries to {target}")


if __name__ == "__main__":
    main()
