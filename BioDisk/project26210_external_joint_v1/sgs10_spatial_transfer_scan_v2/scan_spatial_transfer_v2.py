from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import socket
import sys
import traceback
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

FORMAL_ROOT = Path(sys.argv[1]).resolve()
PART4_ROOT = Path(sys.argv[2]).resolve()
OUT_ROOT = Path(sys.argv[3]).resolve()

EXPECTED_H5AD = {
    "GLMF1": PART4_ROOT / "analysis/GLMF1_analyzed.h5ad",
    "GLMF2": PART4_ROOT / "analysis/GLMF2_analyzed.h5ad",
}

CONTROL_FILES = [
    FORMAL_ROOT / "AGENTS.md",
    FORMAL_ROOT / "PROJECT_STATUS.md",
    FORMAL_ROOT / "RESULTS_MANIFEST.tsv",
]

TABULAR_SUFFIXES = {".csv", ".tsv", ".txt", ".json", ".parquet", ".feather"}
SCRIPT_SUFFIXES = {".py", ".r", ".sh", ".R"}
IMAGE_SUFFIXES = {".png", ".jpg", ".jpeg", ".tif", ".tiff", ".svs"}
DATA_SUFFIXES = TABULAR_SUFFIXES | SCRIPT_SUFFIXES | IMAGE_SUFFIXES | {
    ".h5ad", ".h5", ".hdf5", ".mtx", ".gz", ".npz", ".npy", ".rds"
}

PRUNE_DIRS = {
    ".git", "__pycache__", ".ipynb_checkpoints", ".venv", "venv",
    "node_modules", "site-packages", "scdrs", "external", "tmp", "temp",
}

LEGACY_TERMS = {
    "backup", "bak", "old", "legacy", "deprecated", "archive",
    "obsolete", "superseded",
}

SOURCE_CATEGORIES = {
    "raw_counts_container", "raw_counts_matrix", "barcodes", "features",
    "spatial_coordinates", "spatial_scalefactors", "tissue_hires_image",
    "tissue_lowres_image",
}

RESULT_TABLE_CATEGORIES = {
    "cell_abundance_table", "niche_result", "pathology_annotation",
    "tumor_nat_distance", "spatial_association_result", "misty_result",
    "gsva_result",
}

FINAL_TABLE_TERMS = {
    "final", "summary", "stats", "statistic", "correlation", "colocal",
    "association", "relation_test", "spatial_lag", "score", "annotation",
    "assignment", "distance", "result",
}

INTERMEDIATE_TERMS = {
    "debug", "checkpoint", "cache", "temporary", "intermediate", "permutation_null",
}


def safe_relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(FORMAL_ROOT))
    except Exception:
        return ""


def under_formal(path: Path) -> bool:
    try:
        path.resolve().relative_to(FORMAL_ROOT)
        return True
    except Exception:
        return False


def sample_from_text(text: str) -> str:
    upper = text.upper()
    found = [sample for sample in ("GLMF1", "GLMF2") if sample in upper]
    if len(found) == 1:
        return found[0]
    return "SHARED"


def is_legacy(path: Path) -> bool:
    components = [x.lower() for x in path.parts]
    return any(any(term == comp or term in comp for term in LEGACY_TERMS) for comp in components)


def is_excluded_path(path: Path) -> bool:
    lower_parts = [x.lower() for x in path.parts]
    for comp in lower_parts:
        if comp in PRUNE_DIRS:
            return True
        if re.fullmatch(r"part[_-]?[123].*", comp):
            return True
    return False


def category_for(path: Path) -> str:
    name = path.name.lower()
    text = str(path).lower().replace("\\", "/")

    if path.resolve() in {x.resolve() for x in EXPECTED_H5AD.values()}:
        return "analyzed_h5ad"
    if name.endswith(".h5ad"):
        return "other_h5ad"

    if name.endswith(("filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5")):
        return "raw_counts_container"
    if name in {"matrix.mtx", "matrix.mtx.gz"}:
        return "raw_counts_matrix"
    if name in {"barcodes.tsv", "barcodes.tsv.gz"}:
        return "barcodes"
    if name in {"features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"}:
        return "features"
    if "tissue_positions" in name and name.endswith((".csv", ".csv.gz", ".parquet")):
        return "spatial_coordinates"
    if "scalefactors" in name and name.endswith(".json"):
        return "spatial_scalefactors"
    if "tissue_hires_image" in name:
        return "tissue_hires_image"
    if "tissue_lowres_image" in name:
        return "tissue_lowres_image"

    if any(term in text for term in ("patholog", "histolog", "tumor_nat", "tumor-nat", "tumor_nat")):
        if "distance" in name or "signed_dist" in name or "boundary" in name:
            return "tumor_nat_distance"
        return "pathology_annotation"
    if ("distance" in name or "signed_dist" in name) and any(
        term in text for term in ("tumor", "nat", "boundary")
    ):
        return "tumor_nat_distance"
    if "misty" in text:
        return "misty_result"
    if "gsva" in text:
        return "gsva_result"
    if any(term in name for term in ("colocal", "correlation", "association", "relation_test", "spatial_lag")):
        return "spatial_association_result"
    if "niche" in text:
        return "niche_result"
    if any(term in text for term in ("cell_abundance", "abundance_clean", "cell2location_abundance")):
        return "cell_abundance_table"
    if path.suffix in SCRIPT_SUFFIXES or path.suffix.lower() in {x.lower() for x in SCRIPT_SUFFIXES}:
        return "analysis_script"
    return "other"


def eligible_part4_file(path: Path) -> bool:
    if not path.is_file() or is_excluded_path(path) or is_legacy(path):
        return False
    if path.suffix.lower() not in {x.lower() for x in DATA_SUFFIXES} and not path.name.lower().endswith(
        (".mtx.gz", ".tsv.gz", ".csv.gz")
    ):
        return False
    category = category_for(path)
    return category in ({"analyzed_h5ad", "analysis_script"} | SOURCE_CATEGORIES | RESULT_TABLE_CATEGORIES)


def file_snapshot(root: Path):
    records = {}
    errors = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in {".git", "__pycache__"}]
        for filename in filenames:
            path = Path(dirpath) / filename
            try:
                st = path.stat()
                rel = str(path.relative_to(root))
                records[rel] = {
                    "size": int(st.st_size),
                    "mtime_ns": int(st.st_mtime_ns),
                    "dev": int(st.st_dev),
                    "inode": int(st.st_ino),
                }
            except Exception as exc:
                errors.append({"path": str(path), "error": f"{type(exc).__name__}: {exc}"})
    return records, errors


def write_snapshot(path: Path, records, errors):
    with path.open("w", encoding="utf-8") as out:
        json.dump({"files": records, "errors": errors}, out, ensure_ascii=False, indent=2)


def compare_snapshots(before, after):
    before_paths = set(before)
    after_paths = set(after)
    deleted = set(before_paths - after_paths)
    added = set(after_paths - before_paths)
    modified = {
        p for p in before_paths & after_paths
        if before[p] != after[p]
    }

    deleted_by_id = defaultdict(list)
    added_by_id = defaultdict(list)
    for p in deleted:
        r = before[p]
        deleted_by_id[(r["dev"], r["inode"], r["size"], r["mtime_ns"])].append(p)
    for p in added:
        r = after[p]
        added_by_id[(r["dev"], r["inode"], r["size"], r["mtime_ns"])].append(p)

    moved = []
    for key in set(deleted_by_id) & set(added_by_id):
        old_paths = sorted(deleted_by_id[key])
        new_paths = sorted(added_by_id[key])
        for old, new in zip(old_paths, new_paths):
            moved.append((old, new))
            deleted.discard(old)
            added.discard(new)

    return {
        "added": sorted(added),
        "deleted": sorted(deleted),
        "modified": sorted(modified),
        "moved": sorted(moved),
    }


def _safe_candidate_exists(candidate):
    """Return False for malformed or OS-invalid candidate paths."""
    try:
        return candidate.exists()
    except (OSError, ValueError):
        return False


def extract_existing_paths_from_controls():
    referenced = []
    excerpts = []
    patterns = [
        re.compile(r"`([^`]+)`"),
        re.compile(r"['\"]([^'\"]+)['\"]"),
        re.compile(r"(/mnt/BioDisk/[^\t\r\n|;,]+)"),
        re.compile(r"(results/part4_spatial/[^\t\r\n|;,]+)"),
    ]

    for control in CONTROL_FILES:
        if not control.is_file():
            excerpts.append(f"NOT_FOUND\t{control}")
            continue
        text = control.read_text(encoding="utf-8", errors="replace")
        for lineno, line in enumerate(text.splitlines(), start=1):
            lower = line.lower()
            if not any(term in lower for term in ("part4", "spatial", "glm", "cell2loc", "niche", "misty", "gsva", "tymp")):
                continue
            excerpts.append(f"{control}:{lineno}\t{line}")
            hint = sample_from_text(line)
            for pattern in patterns:
                for raw in pattern.findall(line):
                    raw = raw.strip().strip(".()[]{}<>").rstrip(".:,")
                    if not raw:
                        continue
                    candidate = Path(raw)
                    if not candidate.is_absolute():
                        candidate = FORMAL_ROOT / candidate
                    try:
                        candidate = candidate.resolve()
                    except Exception:
                        continue
                    if not under_formal(candidate) or not _safe_candidate_exists(candidate) or is_excluded_path(candidate):
                        continue
                    if PART4_ROOT not in candidate.parents and candidate != PART4_ROOT:
                        text2 = str(candidate).lower()
                        if not (hint in {"GLMF1", "GLMF2"} and any(x in text2 for x in ("spatial", "visium", "outs"))):
                            continue
                    referenced.append((candidate, hint, f"control:{control.name}:{lineno}"))

    return referenced, excerpts


def inventory_candidates():
    rows_by_path = {}
    scan_roots = []

    def add_file(path: Path, evidence: str, sample_hint: str = "SHARED"):
        try:
            path = path.resolve()
            if not eligible_part4_file(path) and category_for(path) not in SOURCE_CATEGORIES:
                return
            if not under_formal(path) or is_excluded_path(path) or is_legacy(path):
                return
            st = path.stat()
        except Exception:
            return
        sample = sample_from_text(str(path))
        if sample == "SHARED" and sample_hint in {"GLMF1", "GLMF2"}:
            sample = sample_hint
        category = category_for(path)
        current = rows_by_path.get(str(path))
        if current:
            evidence_set = set(current["evidence"].split(";"))
            evidence_set.add(evidence)
            current["evidence"] = ";".join(sorted(x for x in evidence_set if x))
            if current["sample"] == "SHARED" and sample != "SHARED":
                current["sample"] = sample
            return
        rows_by_path[str(path)] = {
            "absolute_path": str(path),
            "relative_to_project": safe_relative(path),
            "sample": sample,
            "data_category": category,
            "size_bytes": int(st.st_size),
            "size_mib": round(st.st_size / (1024 ** 2), 4),
            "modified_time": datetime.fromtimestamp(st.st_mtime).isoformat(),
            "evidence": evidence,
        }

    # 1. Strictly scan the formal Part 4 result tree.
    scan_roots.append((str(PART4_ROOT), "part4_whitelist"))
    for dirpath, dirnames, filenames in os.walk(PART4_ROOT):
        dirnames[:] = [d for d in dirnames if d.lower() not in PRUNE_DIRS]
        for filename in filenames:
            path = Path(dirpath) / filename
            if eligible_part4_file(path):
                add_file(path, "part4_whitelist")

    # 2. Add only existing paths explicitly referenced by the three control files.
    referenced, excerpts = extract_existing_paths_from_controls()
    for ref, hint, evidence in referenced:
        scan_roots.append((str(ref), evidence))
        if ref.is_file():
            add_file(ref, evidence, hint)
        elif ref.is_dir():
            for dirpath, dirnames, filenames in os.walk(ref):
                dirnames[:] = [d for d in dirnames if d.lower() not in PRUNE_DIRS]
                for filename in filenames:
                    path = Path(dirpath) / filename
                    if eligible_part4_file(path) or category_for(path) in SOURCE_CATEGORIES:
                        add_file(path, evidence, hint)

    # 3. Search the formal project only for GLMF1/GLMF2-specific Visium source components.
    # This walk does not admit unrelated files; only an exact source component is retained.
    for dirpath, dirnames, filenames in os.walk(FORMAL_ROOT):
        current = Path(dirpath)
        dirnames[:] = [
            d for d in dirnames
            if d.lower() not in PRUNE_DIRS
            and not re.fullmatch(r"part[_-]?[123].*", d.lower())
        ]
        current_sample = sample_from_text(str(current))
        if current_sample not in {"GLMF1", "GLMF2"}:
            continue
        for filename in filenames:
            path = current / filename
            if category_for(path) in SOURCE_CATEGORIES:
                add_file(path, "glm_sample_visium_source_signature", current_sample)
                scan_roots.append((str(current), f"{current_sample}_source_signature"))

    rows = sorted(rows_by_path.values(), key=lambda x: (x["sample"], x["data_category"], x["absolute_path"]))
    return rows, sorted(set(scan_roots)), excerpts


def matrix_sample_stats(matrix, n_rows=300, n_cols=2000):
    result = {
        "minimum": "", "maximum": "", "nonzero": "", "noninteger": "",
        "negative": "", "nonfinite": "", "error": "",
    }
    try:
        subset = matrix[:min(matrix.shape[0], n_rows), :min(matrix.shape[1], n_cols)]
        values = subset.data if sp.issparse(subset) else np.asarray(subset).ravel()
        values = np.asarray(values, dtype=float)
        if values.size:
            result.update({
                "minimum": float(np.nanmin(values)),
                "maximum": float(np.nanmax(values)),
                "nonzero": int(np.sum(values != 0)),
                "noninteger": int(np.sum(np.isfinite(values) & (np.abs(values - np.round(values)) > 1e-8))),
                "negative": int(np.sum(np.isfinite(values) & (values < 0))),
                "nonfinite": int(np.sum(~np.isfinite(values))),
            })
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
    return result


def inspect_h5ad(sample: str, path: Path):
    row = {
        "sample": sample,
        "absolute_path": str(path),
        "exists": path.is_file(),
        "size_bytes": path.stat().st_size if path.is_file() else "",
        "read_status": "FAIL",
        "error": "",
    }
    if not path.is_file():
        row["error"] = "FILE_NOT_FOUND"
        return row
    try:
        obj = ad.read_h5ad(path, backed="r")
        obs_fields = [str(x) for x in obj.obs.columns]
        layers = [str(x) for x in obj.layers.keys()]
        obsm = [str(x) for x in obj.obsm.keys()]
        uns = [str(x) for x in obj.uns.keys()]
        stats = matrix_sample_stats(obj.X)

        count_layers = [x for x in layers if "count" in x.lower()]
        if count_layers:
            raw_status = "CONFIRMED_COUNTS_LAYER_PRESENT"
        elif stats["error"]:
            raw_status = "UNDETERMINED_X_NOT_AUDITED"
        elif stats["negative"] == 0 and stats["nonfinite"] == 0 and stats["noninteger"] == 0:
            raw_status = "POSSIBLE_INTEGER_X_SAMPLED_NOT_CONFIRMED"
        else:
            raw_status = "X_SAMPLE_NOT_RAW_INTEGER_COUNTS"

        abundance_keys = [
            key for key in obsm
            if any(term in key.lower() for term in ("abundance", "means", "q05", "q95", "stds", "cell2loc"))
        ]
        abundance_details = {}
        for key in abundance_keys:
            try:
                value = obj.obsm[key]
                abundance_details[key] = {
                    "shape": list(value.shape),
                    "columns": [str(x) for x in value.columns] if isinstance(value, pd.DataFrame) else [],
                }
            except Exception as exc:
                abundance_details[key] = {"error": f"{type(exc).__name__}: {exc}"}

        images = []
        try:
            spatial_uns = obj.uns.get("spatial", {})
            if hasattr(spatial_uns, "items"):
                for library_id, value in spatial_uns.items():
                    if isinstance(value, dict) and isinstance(value.get("images"), dict):
                        for image_key in value["images"]:
                            images.append(f"{library_id}:{image_key}")
        except Exception:
            pass

        def match_fields(terms):
            return [field for field in obs_fields if any(term in field.lower() for term in terms)]

        row.update({
            "read_status": "PASS",
            "n_obs": int(obj.n_obs),
            "n_vars": int(obj.n_vars),
            "obs_fields": ";".join(obs_fields),
            "layers": ";".join(layers),
            "obsm_keys": ";".join(obsm),
            "uns_keys": ";".join(uns),
            "raw_slot_present": bool(obj.raw is not None),
            "raw_counts_status": raw_status,
            "x_sample_stats": json.dumps(stats, ensure_ascii=False, separators=(",", ":")),
            "spatial_coordinates": "spatial" in obsm,
            "embedded_images": ";".join(images),
            "abundance_keys": ";".join(abundance_keys),
            "abundance_details": json.dumps(abundance_details, ensure_ascii=False, separators=(",", ":")),
            "niche_fields": ";".join(match_fields(("niche",))),
            "pathology_fields": ";".join(match_fields(("patholog", "histolog", "region"))),
            "tumor_nat_fields": ";".join(match_fields(("tumor_nat", "tumor.nat", "tumor-nat", "tissue_status"))),
            "distance_fields": ";".join(match_fields(("distance", "signed_dist", "dist_"))),
        })
        if getattr(obj, "file", None) is not None:
            obj.file.close()
    except Exception as exc:
        row["error"] = f"{type(exc).__name__}: {exc}"
        row["traceback_tail"] = traceback.format_exc()[-2000:]
    return row


def is_final_small_table(row) -> bool:
    path = Path(row["absolute_path"])
    if path.suffix.lower() not in TABULAR_SUFFIXES:
        return False
    text = str(path).lower()
    if any(term in text for term in INTERMEDIATE_TERMS):
        return False
    return any(term in path.name.lower() for term in FINAL_TABLE_TERMS)


def classify_current(row):
    path = Path(row["absolute_path"])
    category = row["data_category"]
    if path.resolve() in {x.resolve() for x in EXPECTED_H5AD.values()}:
        return "P0_REQUIRED", "core analyzed H5AD for existing cell2location/spatial validation"
    if category in {
        "spatial_association_result", "pathology_annotation", "tumor_nat_distance",
        "misty_result", "gsva_result",
    } and is_final_small_table(row):
        return "P1_RECOMMENDED", "final small result/annotation table not safely replaced by H5AD alone"
    if category in {"niche_result", "cell_abundance_table"} and is_final_small_table(row):
        return "P2_OPTIONAL", "usually embedded in analyzed H5AD; retain only for traceability"
    if category == "analysis_script":
        return "P2_OPTIONAL", "analysis provenance; not required for immediate validation"
    return "EXCLUDE", "not required for current-result validation"


def classify_remap(row):
    category = row["data_category"]
    if category in SOURCE_CATEGORIES:
        return "P0_REQUIRED", "Visium counts/coordinates/scalefactors/image component for future remapping"
    if category in {"analyzed_h5ad", "pathology_annotation", "tumor_nat_distance"}:
        return "P1_RECOMMENDED", "current reference result or pathology context for remapping comparison"
    if category == "analysis_script":
        return "P2_OPTIONAL", "method provenance"
    return "EXCLUDE", "not required as a remapping input"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_recommendations(candidate_rows):
    recs = []
    for scenario, classifier in (
        ("CURRENT_RESULT_VALIDATION", classify_current),
        ("NEW_REFERENCE_REMAPPING", classify_remap),
    ):
        for row in candidate_rows:
            priority, reason = classifier(row)
            recs.append({
                "scenario": scenario,
                "priority": priority,
                "sample": row["sample"],
                "data_category": row["data_category"],
                "absolute_path": row["absolute_path"],
                "relative_to_project": row["relative_to_project"],
                "size_bytes": row["size_bytes"],
                "sha256": "",
                "duplicate_of": "",
                "reason": reason,
                "evidence": row["evidence"],
            })

    selected_paths = sorted({
        rec["absolute_path"] for rec in recs
        if rec["priority"] in {"P0_REQUIRED", "P1_RECOMMENDED"}
    })
    hashes = {}
    for index, path_text in enumerate(selected_paths, start=1):
        print(f"[SHA256 {index}/{len(selected_paths)}] {path_text}", flush=True)
        hashes[path_text] = sha256_file(Path(path_text))

    grouped = defaultdict(list)
    for path_text, digest in hashes.items():
        grouped[(digest, Path(path_text).stat().st_size)].append(path_text)

    duplicate_of = {}
    for _, paths in grouped.items():
        if len(paths) < 2:
            continue
        canonical = sorted(
            paths,
            key=lambda x: (
                0 if Path(x).resolve() in {p.resolve() for p in EXPECTED_H5AD.values()} else 1,
                0 if str(PART4_ROOT) in x else 1,
                len(x),
                x,
            ),
        )[0]
        for path_text in paths:
            if path_text != canonical:
                duplicate_of[path_text] = canonical

    for rec in recs:
        rec["sha256"] = hashes.get(rec["absolute_path"], "NOT_CALCULATED_NOT_SELECTED")
        if rec["absolute_path"] in duplicate_of and rec["priority"] in {"P0_REQUIRED", "P1_RECOMMENDED"}:
            rec["duplicate_of"] = duplicate_of[rec["absolute_path"]]
            rec["priority"] = "EXCLUDE_DUPLICATE"
            rec["reason"] = "byte-identical duplicate of canonical selected file"
    return recs


def build_component_matrix(candidate_rows, audits):
    audit_map = {row["sample"]: row for row in audits}
    result = []
    for sample in ("GLMF1", "GLMF2"):
        by_category = defaultdict(list)
        for row in candidate_rows:
            if row["sample"] == sample:
                by_category[row["data_category"]].append(row["absolute_path"])
        h5_container = bool(by_category["raw_counts_container"])
        mtx_triplet = bool(by_category["raw_counts_matrix"] and by_category["barcodes"] and by_category["features"])
        counts_complete = h5_container or mtx_triplet
        spatial_complete = bool(
            by_category["spatial_coordinates"]
            and by_category["spatial_scalefactors"]
            and (by_category["tissue_hires_image"] or by_category["tissue_lowres_image"])
        )
        audit = audit_map.get(sample, {})
        result.append({
            "sample": sample,
            "analyzed_h5ad": str(EXPECTED_H5AD[sample]) if EXPECTED_H5AD[sample].is_file() else "NOT_FOUND",
            "h5ad_dimensions": f"{audit.get('n_obs', '')}x{audit.get('n_vars', '')}" if audit.get("read_status") == "PASS" else "UNAVAILABLE",
            "raw_counts_container": ";".join(by_category["raw_counts_container"]) or "NOT_FOUND",
            "matrix_mtx": ";".join(by_category["raw_counts_matrix"]) or "NOT_FOUND",
            "barcodes": ";".join(by_category["barcodes"]) or "NOT_FOUND",
            "features": ";".join(by_category["features"]) or "NOT_FOUND",
            "spatial_coordinates": ";".join(by_category["spatial_coordinates"]) or "NOT_FOUND",
            "spatial_scalefactors": ";".join(by_category["spatial_scalefactors"]) or "NOT_FOUND",
            "hires_image": ";".join(by_category["tissue_hires_image"]) or "NOT_FOUND",
            "lowres_image": ";".join(by_category["tissue_lowres_image"]) or "NOT_FOUND",
            "counts_bundle_complete": counts_complete,
            "spatial_support_complete": spatial_complete,
            "future_remap_bundle_complete": counts_complete and spatial_complete,
            "h5ad_spatial_embedded": audit.get("spatial_coordinates", ""),
            "h5ad_images_embedded": audit.get("embedded_images", ""),
            "h5ad_abundance_keys": audit.get("abundance_keys", ""),
            "h5ad_niche_fields": audit.get("niche_fields", ""),
            "h5ad_pathology_fields": audit.get("pathology_fields", ""),
            "h5ad_tumor_nat_fields": audit.get("tumor_nat_fields", ""),
            "h5ad_distance_fields": audit.get("distance_fields", ""),
        })
    return result


def write_tsv(path: Path, rows):
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def ps_quote(value: str) -> str:
    return value.replace("'", "''")


def write_powershell(path: Path, records, title: str, local_folder: str):
    lines = [
        "$ErrorActionPreference = 'Stop'",
        f"$DestinationRoot = Join-Path 'C:\\Users\\31937\\Downloads' '{ps_quote(local_folder)}'",
        "New-Item -ItemType Directory -Force -Path $DestinationRoot | Out-Null",
        "$Files = @(",
    ]
    for rec in records:
        rel = rec["relative_to_project"] or Path(rec["absolute_path"]).name
        lines.extend([
            "    [pscustomobject]@{",
            f"        Remote = '{ps_quote(rec['absolute_path'])}'",
            f"        Relative = '{ps_quote(rel.replace('/', chr(92)))}'",
            f"        Bytes = [int64]{rec['size_bytes']}",
            f"        SHA256 = '{ps_quote(rec['sha256'])}'",
            "    }",
        ])
    lines.extend([
        ")",
        "",
        f"Write-Host '{ps_quote(title)}'",
        "Write-Host ('Files: {0}' -f $Files.Count)",
        "Write-Host ('Bytes: {0}' -f (($Files | Measure-Object -Property Bytes -Sum).Sum))",
        "",
        "foreach ($File in $Files) {",
        "    $LocalPath = Join-Path $DestinationRoot $File.Relative",
        "    $LocalDir = Split-Path -Parent $LocalPath",
        "    New-Item -ItemType Directory -Force -Path $LocalDir | Out-Null",
        "    Write-Host ('Downloading: {0}' -f $File.Remote)",
        "    & scp -- ('sgs10:' + $File.Remote) $LocalPath",
        "    if ($LASTEXITCODE -ne 0) { throw ('scp failed: {0}' -f $File.Remote) }",
        "    $ActualBytes = (Get-Item -LiteralPath $LocalPath).Length",
        "    if ($ActualBytes -ne $File.Bytes) { throw ('size mismatch: {0}' -f $LocalPath) }",
        "    $ActualHash = (Get-FileHash -Algorithm SHA256 -LiteralPath $LocalPath).Hash.ToLowerInvariant()",
        "    if ($ActualHash -ne $File.SHA256.ToLowerInvariant()) { throw ('SHA256 mismatch: {0}' -f $LocalPath) }",
        "}",
        "Write-Host 'DOWNLOAD_AND_SHA256_CHECK=PASS'",
        "",
    ])
    path.write_text("\r\n".join(lines), encoding="utf-8-sig")


def human_size(size: int) -> str:
    value = float(size)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024 or unit == "TiB":
            return f"{value:.2f} {unit}"
        value /= 1024
    return str(size)


def main():
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    before, before_errors = file_snapshot(FORMAL_ROOT)
    write_snapshot(OUT_ROOT / "formal_project_snapshot_before.json", before, before_errors)

    candidates, scan_roots, excerpts = inventory_candidates()
    write_tsv(OUT_ROOT / "all_candidate_files.tsv", candidates)
    write_tsv(
        OUT_ROOT / "scan_roots.tsv",
        [{"path": path, "reason": reason} for path, reason in scan_roots],
    )
    (OUT_ROOT / "project_control_references.txt").write_text("\n".join(excerpts) + "\n", encoding="utf-8")

    audits = []
    for sample, path in EXPECTED_H5AD.items():
        print(f"[H5AD] {sample} {path}", flush=True)
        audits.append(inspect_h5ad(sample, path))
    write_tsv(OUT_ROOT / "h5ad_content_audit.tsv", audits)

    components = build_component_matrix(candidates, audits)
    write_tsv(OUT_ROOT / "spatial_component_matrix.tsv", components)

    recommendations = build_recommendations(candidates)
    write_tsv(OUT_ROOT / "transfer_recommendations.tsv", recommendations)

    current_selected = [
        rec for rec in recommendations
        if rec["scenario"] == "CURRENT_RESULT_VALIDATION"
        and rec["priority"] in {"P0_REQUIRED", "P1_RECOMMENDED"}
    ]
    minimal_tables = [
        rec for rec in current_selected
        if rec["priority"] == "P1_RECOMMENDED"
        and Path(rec["absolute_path"]).suffix.lower() in TABULAR_SUFFIXES
        and int(rec["size_bytes"]) <= 50 * 1024 * 1024
    ]

    ps_full = OUT_ROOT / "download_spatial_p0_p1.ps1"
    ps_tables = OUT_ROOT / "download_spatial_minimal_tables.ps1"
    write_powershell(
        ps_full,
        current_selected,
        "project26210 current spatial validation: P0 plus P1",
        "project26210_spatial_p0_p1",
    )
    write_powershell(
        ps_tables,
        minimal_tables,
        "project26210 supplemental minimal spatial tables only",
        "project26210_spatial_minimal_tables",
    )

    after, after_errors = file_snapshot(FORMAL_ROOT)
    write_snapshot(OUT_ROOT / "formal_project_snapshot_after.json", after, after_errors)
    changes = compare_snapshots(before, after)
    change_rows = []
    for p in changes["added"]:
        change_rows.append({"change": "ADDED", "old_path": "", "new_path": p})
    for p in changes["deleted"]:
        change_rows.append({"change": "DELETED", "old_path": p, "new_path": ""})
    for p in changes["modified"]:
        change_rows.append({"change": "MODIFIED_METADATA", "old_path": p, "new_path": p})
    for old, new in changes["moved"]:
        change_rows.append({"change": "MOVED_SAME_INODE", "old_path": old, "new_path": new})
    write_tsv(OUT_ROOT / "formal_project_changes.tsv", change_rows)
    project_unchanged = not any(changes.values())

    priority_summary = defaultdict(lambda: {"count": 0, "bytes": 0})
    for rec in recommendations:
        key = (rec["scenario"], rec["priority"])
        priority_summary[key]["count"] += 1
        priority_summary[key]["bytes"] += int(rec["size_bytes"])

    duplicates = sorted({
        (rec["absolute_path"], rec["duplicate_of"])
        for rec in recommendations if rec["duplicate_of"]
    })

    not_confirmed = []
    for component in components:
        sample = component["sample"]
        if not component["counts_bundle_complete"]:
            not_confirmed.append(f"{sample}: complete raw-count bundle not confirmed")
        if not component["spatial_support_complete"]:
            not_confirmed.append(f"{sample}: complete coordinates/scalefactors/image bundle not confirmed")
        if not component["future_remap_bundle_complete"]:
            not_confirmed.append(f"{sample}: future remapping bundle incomplete or unconfirmed")
    for audit in audits:
        if audit.get("raw_counts_status") != "CONFIRMED_COUNTS_LAYER_PRESENT":
            not_confirmed.append(
                f"{audit['sample']}: analyzed H5AD raw counts not confirmed; status={audit.get('raw_counts_status', 'UNKNOWN')}"
            )

    summary_path = OUT_ROOT / "transfer_summary.md"
    with summary_path.open("w", encoding="utf-8") as out:
        out.write("# project26210 空间数据精准传输扫描 v2\n\n")
        out.write(f"- 扫描时间：{datetime.now().isoformat()}\n")
        out.write(f"- 服务器：{socket.gethostname()}\n")
        out.write(f"- 正式项目：`{FORMAL_ROOT}`\n")
        out.write(f"- 白名单第四部分：`{PART4_ROOT}`\n")
        out.write(f"- 候选文件数：{len(candidates)}\n")
        out.write(f"- 正式项目元数据快照是否保持不变：{'YES' if project_unchanged else 'NO'}\n")
        out.write("- 说明：变更检查比较路径、大小、mtime_ns、设备号和inode；不等同于对正式项目全部文件计算内容哈希。\n\n")

        out.write("## analyzed H5AD\n\n")
        out.write("| 样本 | 大小 | 维度 | 坐标 | 图像 | 丰度 | niche | raw counts判断 |\n")
        out.write("|---|---:|---:|---|---|---|---|---|\n")
        for audit in audits:
            dims = f"{audit.get('n_obs', '')} × {audit.get('n_vars', '')}" if audit.get("read_status") == "PASS" else "READ_FAIL"
            out.write(
                f"| {audit['sample']} | {human_size(int(audit['size_bytes'])) if audit.get('size_bytes') != '' else 'NA'} | "
                f"{dims} | {audit.get('spatial_coordinates', '')} | {audit.get('embedded_images', '') or 'NONE'} | "
                f"{audit.get('abundance_keys', '') or 'NONE'} | {audit.get('niche_fields', '') or 'NONE'} | "
                f"{audit.get('raw_counts_status', '')} |\n"
            )

        out.write("\n## 下载分级\n\n")
        out.write("| 场景 | 优先级 | 文件数 | 总大小 |\n")
        out.write("|---|---|---:|---:|\n")
        for key in sorted(priority_summary):
            data = priority_summary[key]
            out.write(f"| {key[0]} | {key[1]} | {data['count']} | {human_size(data['bytes'])} |\n")

        out.write("\n## 当前联合验证最小集合\n\n")
        for rec in current_selected:
            out.write(
                f"- `{rec['absolute_path']}` — {rec['priority']}; {human_size(int(rec['size_bytes']))}; SHA256 `{rec['sha256']}`\n"
            )
        out.write("\n`download_spatial_p0_p1.ps1` 下载上述核心对象及推荐表；"
                  "`download_spatial_minimal_tables.ps1` 只下载其中不超过50 MiB的P1补充表，不包含H5AD，不能单独完成丰度验证。\n")

        out.write("\n## 将来重新映射完整性\n\n")
        out.write("| 样本 | counts完整 | 空间支持完整 | 完整remap集合 |\n")
        out.write("|---|---|---|---|\n")
        for component in components:
            out.write(
                f"| {component['sample']} | {component['counts_bundle_complete']} | "
                f"{component['spatial_support_complete']} | {component['future_remap_bundle_complete']} |\n"
            )

        out.write("\n## 未找到或无法确认\n\n")
        if not_confirmed:
            for item in not_confirmed:
                out.write(f"- {item}\n")
        else:
            out.write("- NONE\n")

        out.write("\n## 重复对象\n\n")
        if duplicates:
            for duplicate, canonical in duplicates:
                out.write(f"- `{duplicate}` → `{canonical}`\n")
        else:
            out.write("- 未在最终P0/P1集合中发现内容完全相同的重复文件。\n")

        out.write("\n## 正式项目保护检查\n\n")
        out.write(f"- added={len(changes['added'])}\n")
        out.write(f"- deleted={len(changes['deleted'])}\n")
        out.write(f"- modified_metadata={len(changes['modified'])}\n")
        out.write(f"- moved_same_inode={len(changes['moved'])}\n")
        out.write(f"- unchanged={'YES' if project_unchanged else 'NO'}\n")

        out.write("\n## 输出文件\n\n")
        for output in (
            OUT_ROOT / "all_candidate_files.tsv",
            OUT_ROOT / "scan_roots.tsv",
            OUT_ROOT / "h5ad_content_audit.tsv",
            OUT_ROOT / "spatial_component_matrix.tsv",
            OUT_ROOT / "transfer_recommendations.tsv",
            OUT_ROOT / "formal_project_changes.tsv",
            ps_full,
            ps_tables,
        ):
            out.write(f"- `{output}`\n")

    print("\n===== PRECISE SCAN SUMMARY =====")
    print(f"candidate_files={len(candidates)}")
    for audit in audits:
        print(
            f"{audit['sample']}|path={audit['absolute_path']}|size_bytes={audit.get('size_bytes', '')}|"
            f"read_status={audit.get('read_status', '')}|n_obs={audit.get('n_obs', '')}|n_vars={audit.get('n_vars', '')}|"
            f"raw_counts_status={audit.get('raw_counts_status', '')}|spatial={audit.get('spatial_coordinates', '')}|"
            f"images={audit.get('embedded_images', '') or 'NONE'}|abundance={audit.get('abundance_keys', '') or 'NONE'}|"
            f"niche={audit.get('niche_fields', '') or 'NONE'}|pathology={audit.get('pathology_fields', '') or 'NONE'}|"
            f"distance={audit.get('distance_fields', '') or 'NONE'}"
        )

    print("\n===== PRIORITY COUNTS =====")
    for key in sorted(priority_summary):
        data = priority_summary[key]
        print(f"{key[0]}|{key[1]}|files={data['count']}|bytes={data['bytes']}")

    print("\n===== CURRENT P0 P1 FILES =====")
    for rec in current_selected:
        print(
            f"{rec['priority']}|{rec['sample']}|{rec['data_category']}|{rec['size_bytes']}|"
            f"{rec['sha256']}|{rec['absolute_path']}"
        )

    print("\n===== REMAP COMPLETENESS =====")
    for component in components:
        print(
            f"{component['sample']}|counts_bundle_complete={component['counts_bundle_complete']}|"
            f"spatial_support_complete={component['spatial_support_complete']}|"
            f"future_remap_bundle_complete={component['future_remap_bundle_complete']}"
        )

    print("\n===== DUPLICATES =====")
    if duplicates:
        for duplicate, canonical in duplicates:
            print(f"DUPLICATE|{duplicate}|canonical={canonical}")
    else:
        print("duplicate_files=0")

    print("\n===== FORMAL PROJECT SNAPSHOT COMPARISON =====")
    print(f"added={len(changes['added'])}")
    print(f"deleted={len(changes['deleted'])}")
    print(f"modified_metadata={len(changes['modified'])}")
    print(f"moved_same_inode={len(changes['moved'])}")
    print(f"formal_project_unchanged={'YES' if project_unchanged else 'NO'}")

    print("\n===== OUTPUT PATHS =====")
    print(f"summary={summary_path}")
    print(f"download_p0_p1={ps_full}")
    print(f"download_minimal_tables={ps_tables}")
    print("download_started=NO")


if __name__ == "__main__":
    main()
