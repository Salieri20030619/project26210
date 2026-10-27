#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import pickle
from datetime import datetime
from pathlib import Path


FROZEN_TFS = {
    "myeloid_dc_tymp": ["IRF1", "SPI1", "RELB", "STAT6", "NFATC1", "ELF4", "CUX1", "IRF7", "STAT2", "ELF2", "CEBPA", "ZBTB7A", "NFKB1", "CEBPB", "RELA"],
    "treg_tymp": ["GATA3", "IRF2", "NFKB2", "REL", "TCF3", "RELB", "STAT1", "TFEB"],
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


parser = argparse.ArgumentParser()
parser.add_argument("--myeloid", required=True, type=Path)
parser.add_argument("--treg", required=True, type=Path)
parser.add_argument("--output", required=True, type=Path)
args = parser.parse_args()

part = args.output.with_name(args.output.name + ".part")
if args.output.exists() or part.exists():
    raise SystemExit("Refusing existing output directory")
part.mkdir(parents=True)

sources = {"myeloid_dc_tymp": args.myeloid, "treg_tymp": args.treg}
rows: list[dict] = []
summary: list[dict] = []
for endpoint, path in sources.items():
    with path.open("rb") as handle:
        regulons = pickle.load(handle)
    by_tf = {}
    for regulon in regulons:
        name = str(regulon.name)
        tf = name.rsplit("(", 1)[0]
        if tf in FROZEN_TFS[endpoint] and name.endswith("(+)"):
            if tf in by_tf:
                raise RuntimeError(f"Duplicate positive regulon for {endpoint} {tf}")
            by_tf[tf] = regulon
    missing = sorted(set(FROZEN_TFS[endpoint]) - set(by_tf))
    if missing:
        raise RuntimeError(f"Missing frozen regulons for {endpoint}: {missing}")
    for order, tf in enumerate(FROZEN_TFS[endpoint], start=1):
        regulon = by_tf[tf]
        weights = dict(regulon.gene2weight)
        if len(weights) < 5:
            raise RuntimeError(f"Too few targets for {endpoint} {tf}: {len(weights)}")
        for gene, weight in sorted(weights.items()):
            rows.append({
                "endpoint": endpoint,
                "tf_order": order,
                "tf": tf,
                "regulon_name": str(regulon.name),
                "target_gene": str(gene).upper(),
                "weight": float(weight),
            })
        summary.append({
            "endpoint": endpoint,
            "tf_order": order,
            "tf": tf,
            "regulon_name": str(regulon.name),
            "n_targets": len(weights),
        })

write_tsv(part / "phase64_frozen_regulon_targets_v1.tsv", rows, ["endpoint", "tf_order", "tf", "regulon_name", "target_gene", "weight"])
write_tsv(part / "phase64_frozen_regulon_summary_v1.tsv", summary, ["endpoint", "tf_order", "tf", "regulon_name", "n_targets"])
metadata = {
    "status": "PASS_FROZEN_REGULON_EXTRACTION_V1",
    "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    "myeloid_tfs": len(FROZEN_TFS["myeloid_dc_tymp"]),
    "treg_tfs": len(FROZEN_TFS["treg_tymp"]),
    "target_rows": len(rows),
    "source_sha256": {endpoint: sha256(path) for endpoint, path in sources.items()},
}
(part / "extraction_metadata_v1.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
manifest = []
for path in sorted(part.iterdir()):
    manifest.append({"file": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path)})
write_tsv(part / "FILE_MANIFEST_SHA256_v1.tsv", manifest, ["file", "size_bytes", "sha256"])
(part / "EXTRACTION_PASS_v1").write_text("status=PASS\n", encoding="utf-8")
os.replace(part, args.output)
print(json.dumps(metadata, indent=2))
