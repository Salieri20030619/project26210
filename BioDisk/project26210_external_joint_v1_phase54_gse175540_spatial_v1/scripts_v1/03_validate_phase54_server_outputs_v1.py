#!/usr/bin/env python3
"""Independent post-run validation for the Phase 54 server outputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime
from pathlib import Path

import anndata as ad
import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    temporary = path.with_name(path.name + ".part")
    frame.to_csv(temporary, sep="\t", index=False, encoding="utf-8", lineterminator="\n")
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    run = args.root / "run_v1"
    final = args.root / "postrun_validation_v1"
    part = args.root / "postrun_validation_v1.part"
    if final.exists() or part.exists():
        raise FileExistsError("Refusing to overwrite post-run validation")
    part.mkdir(parents=True)

    status = pd.read_csv(run / "section_status_v1.tsv", sep="\t")
    failures = pd.read_csv(run / "section_failures_v1.tsv", sep="\t")
    relations = pd.read_csv(run / "all_section_spatial_bridge_relations_v1.tsv", sep="\t")
    moran = pd.read_csv(run / "all_section_spatial_autocorrelation_v1.tsv", sep="\t")
    tls = pd.read_csv(run / "all_section_tls_signature_comparison_v1.tsv", sep="\t")
    meta = pd.read_csv(run / "cross_section_relation_meta_v1.tsv", sep="\t")

    object_rows = []
    for sample in sorted(status["sample_id"].astype(str)):
        section = run / "sections_v1" / sample
        h5ad_path = section / "processed_section_v1.h5ad"
        backed = ad.read_h5ad(h5ad_path, backed="r")
        object_rows.append({
            "sample_id": sample,
            "readable": True,
            "n_obs": backed.n_obs,
            "n_vars": backed.n_vars,
            "counts_layer_present": "counts" in backed.layers,
            "spatial_present": "spatial" in backed.obsm,
            "pca_present": "X_pca" in backed.obsm,
            "umap_present": "X_umap" in backed.obsm,
            "cluster_present": "leiden_0_5" in backed.obs,
            "h5ad_bytes": h5ad_path.stat().st_size,
            "h5ad_sha256": sha256(h5ad_path),
            "section_pass_present": (section / "SECTION_PASS_v1").exists(),
            "spatial_plot_present": (section / "target_signature_spatial_maps_v1.png").exists(),
            "cluster_plot_present": (section / "embedding_and_cluster_maps_v1.png").exists(),
        })
        backed.file.close()
    objects = pd.DataFrame(object_rows)
    merged = objects.merge(status[["sample_id", "qc_pass_spots"]], on="sample_id", validate="one_to_one")

    rows = []
    def add(name, passed, observed, expected):
        rows.append({"check": name, "status": "PASS" if passed else "FAIL", "observed": observed, "expected": expected})
    add("section_status_rows", len(status) == 24, len(status), 24)
    add("failure_rows", len(failures) == 0, len(failures), 0)
    add("h5ad_readable_count", len(objects) == 24 and objects["readable"].all(), int(objects["readable"].sum()), 24)
    add("h5ad_spot_counts_match_status", (merged["n_obs"] == merged["qc_pass_spots"]).all(), int((merged["n_obs"] == merged["qc_pass_spots"]).sum()), 24)
    add("required_h5ad_components", objects[["counts_layer_present", "spatial_present", "pca_present", "umap_present", "cluster_present"]].all(axis=None), int(objects[["counts_layer_present", "spatial_present", "pca_present", "umap_present", "cluster_present"]].sum().sum()), 120)
    add("section_pass_markers", objects["section_pass_present"].all(), int(objects["section_pass_present"].sum()), 24)
    add("section_plot_pairs", objects[["spatial_plot_present", "cluster_plot_present"]].all(axis=None), int(objects[["spatial_plot_present", "cluster_plot_present"]].sum().sum()), 48)
    add("bridge_relation_rows", len(relations) == 1056, len(relations), 1056)
    add("spatial_autocorrelation_rows", len(moran) == 864, len(moran), 864)
    add("tls_comparison_rows", len(tls) == 864, len(tls), 864)
    add("cross_section_meta_rows", len(meta) == 176, len(meta), 176)
    add("background_pass_marker", (args.root / "PHASE54_BACKGROUND_PASS_v1").exists(), (args.root / "PHASE54_BACKGROUND_PASS_v1").exists(), True)
    checks = pd.DataFrame(rows)
    atomic_tsv(objects, part / "h5ad_readability_and_hashes_v1.tsv")
    atomic_tsv(checks, part / "postrun_validation_checks_v1.tsv")
    if (checks["status"] != "PASS").any():
        raise RuntimeError("Post-run validation failed")

    deliverables = [
        run / "section_status_v1.tsv", run / "section_failures_v1.tsv",
        run / "all_section_spatial_bridge_relations_v1.tsv",
        run / "all_section_spatial_autocorrelation_v1.tsv",
        run / "all_section_tls_signature_comparison_v1.tsv",
        run / "cross_section_relation_meta_v1.tsv", run / "validation_checks_v1.tsv",
        run / "run_metadata_v1.json", run / "PHASE54_SPATIAL_PROCESSING_PASS_v1",
        args.root / "PHASE54_BACKGROUND_PASS_v1",
    ]
    manifest = pd.DataFrame([
        {"relative_path": str(path.relative_to(args.root)), "bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in deliverables
    ])
    atomic_tsv(manifest, part / "compact_deliverables_manifest_v1.tsv")
    metadata = {
        "status": "PASS", "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "checks_passed": int((checks["status"] == "PASS").sum()), "checks_total": len(checks),
        "h5ad_objects": len(objects), "total_h5ad_bytes": int(objects["h5ad_bytes"].sum()),
    }
    (part / "POSTRUN_VALIDATION_STATUS_v1.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    (part / "POSTRUN_VALIDATION_PASS_v1").write_text(
        f"status=PASS\nchecks={len(checks)}/{len(checks)}\nh5ad_objects={len(objects)}\n",
        encoding="utf-8",
    )
    os.replace(part, final)
    print(f"PASS checks={len(checks)}/{len(checks)} h5ad={len(objects)} bytes={metadata['total_h5ad_bytes']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

