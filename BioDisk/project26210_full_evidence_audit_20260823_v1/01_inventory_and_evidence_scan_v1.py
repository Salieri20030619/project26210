#!/usr/bin/env python3
import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import re
import sys
import tarfile
import zipfile
from collections import Counter, defaultdict
from pathlib import Path


PATTERNS = {
    "MTA": re.compile(r"(?i)(?<![A-Za-z0-9])(?:5[\u2032'\- ]*)?MTA(?![A-Za-z0-9])|methylthioadenosine|5[\u2032'\- ]*methylthioadenosine|5[\u2032'\- ]*\u7532\u57fa\u786b\u4ee3\u817a\u82f7"),
    "MTA_44X": re.compile(r"(?i)44\s*(?:x|fold|\u500d)|(?:x|fold)\s*44|\u72c2\u98d844\s*\u500d|44\s*\u500d.*MTA|MTA.*44\s*\u500d"),
    "TYMP": re.compile(r"(?i)(?<![A-Za-z0-9])TYMP(?![A-Za-z0-9])|thymidine phosphorylase|\u80f8\u82f7\u78f7\u9178\u5316\u9176"),
    "SPP1": re.compile(r"(?i)(?<![A-Za-z0-9])SPP1(?![A-Za-z0-9])|osteopontin|\u9aa8\u6865\u86cb\u767d"),
    "METABOLOMICS": re.compile(r"(?i)metabolomics|metabolome|untargeted metabol|non[- ]targeted metabol|\u975e\u9776\u4ee3\u8c22\u7ec4|\u975e\u9776\u5411\u4ee3\u8c22\u7ec4|\u4ee3\u8c22\u7ec4|\u4ee3\u8c22\u7269"),
    "WET_LAB": re.compile(r"(?i)wet[- ]lab|LC[- ]?MS|MS/MS|mass spectrom|ELISA|western blot|\bWB\b|qPCR|flow cyt|\u6e7f\u5b9e\u9a8c|\u8d28\u8c31|\u6d41\u5f0f|\u514d\u75ab\u7ec4\u5316"),
    "MTAP_POLYAMINE": re.compile(r"(?i)(?<![A-Za-z0-9])MTAP(?![A-Za-z0-9])|polyamine|methionine salvage|\u591a\u80fa|\u86cb\u6c28\u9178\u8865\u6551"),
}

TEXT_EXTS = {
    ".txt", ".tsv", ".csv", ".md", ".json", ".yaml", ".yml", ".xml",
    ".html", ".htm", ".r", ".py", ".sh", ".log", ".ini", ".toml",
}
OFFICE_EXTS = {".docx", ".xlsx", ".pptx"}
ARCHIVE_EXTS = {".zip", ".7z", ".tar", ".tgz", ".gz", ".bz2", ".xz"}
MAX_TEXT_BYTES = 20 * 1024 * 1024
MAX_ZIP_MEMBER_BYTES = 10 * 1024 * 1024


def utc_iso(ts):
    return dt.datetime.fromtimestamp(ts, tz=dt.timezone.utc).isoformat()


def safe_text(raw):
    for enc in ("utf-8", "utf-16", "gb18030", "latin-1"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            pass
    return raw.decode("utf-8", errors="replace")


def classify_path(path_lower):
    if any(k in path_lower for k in ("raw", "original", "source_data", "input", "metabolomics", "metabolome")):
        return "POTENTIAL_SOURCE_OR_INPUT"
    if any(k in path_lower for k in ("result", "output", "figure", "report", "summary", "qc", "validation")):
        return "GENERATED_OR_PROCESSED"
    return "UNCLASSIFIED"


def pattern_hits(text):
    return [name for name, pat in PATTERNS.items() if pat.search(text)]


def matching_contexts(text, max_per_pattern=4, radius=160):
    out = []
    compact = re.sub(r"\s+", " ", text)
    for name, pat in PATTERNS.items():
        count = 0
        for m in pat.finditer(compact):
            lo = max(0, m.start() - radius)
            hi = min(len(compact), m.end() + radius)
            out.append((name, compact[lo:hi]))
            count += 1
            if count >= max_per_pattern:
                break
    return out


def scan_regular_text(path, size):
    if size > MAX_TEXT_BYTES:
        return [], "SKIP_SIZE"
    try:
        raw = path.read_bytes()
        text = safe_text(raw)
        return matching_contexts(text), "OK"
    except Exception as exc:
        return [], f"ERROR:{type(exc).__name__}:{exc}"


def scan_office_zip(path):
    contexts = []
    try:
        with zipfile.ZipFile(path) as zf:
            for info in zf.infolist():
                name = info.filename.lower()
                if not name.endswith((".xml", ".rels")) or info.file_size > MAX_ZIP_MEMBER_BYTES:
                    continue
                text = safe_text(zf.read(info))
                for pat, ctx in matching_contexts(text):
                    contexts.append((pat, f"{info.filename}: {ctx}"))
        return contexts, "OK"
    except Exception as exc:
        return [], f"ERROR:{type(exc).__name__}:{exc}"


def archive_members(path):
    ext = path.suffix.lower()
    members = []
    try:
        if ext == ".zip":
            with zipfile.ZipFile(path) as zf:
                members = [i.filename for i in zf.infolist()]
        elif ext in {".tar", ".tgz", ".gz", ".bz2", ".xz"} and tarfile.is_tarfile(path):
            with tarfile.open(path) as tf:
                members = [m.name for m in tf.getmembers()]
        else:
            return [], "UNSUPPORTED_LISTING"
        return members, "OK"
    except Exception as exc:
        return [], f"ERROR:{type(exc).__name__}:{exc}"


def write_tsv(path, fieldnames, rows):
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", action="append", required=True, help="LABEL=/absolute/path")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    roots = []
    for item in args.root:
        label, value = item.split("=", 1)
        roots.append((label, Path(value)))

    inventory = []
    candidates = []
    content_hits = []
    archive_hits = []
    root_counts = Counter()
    root_bytes = Counter()
    ext_counts = Counter()
    errors = []

    for root_label, root in roots:
        if not root.exists():
            errors.append({"root_label": root_label, "path": str(root), "error": "ROOT_MISSING"})
            continue
        for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
            dirnames[:] = [d for d in dirnames if d not in {".git", "__pycache__", ".cache"}]
            for filename in filenames:
                path = Path(dirpath) / filename
                try:
                    st = path.stat()
                except Exception as exc:
                    errors.append({"root_label": root_label, "path": str(path), "error": repr(exc)})
                    continue
                suffix = path.suffix.lower() or "[none]"
                path_text = str(path)
                name_hits = pattern_hits(path_text)
                row = {
                    "root_label": root_label,
                    "path": path_text,
                    "size_bytes": st.st_size,
                    "mtime_utc": utc_iso(st.st_mtime),
                    "extension": suffix,
                    "path_class": classify_path(path_text.lower()),
                    "name_pattern_hits": ";".join(name_hits),
                }
                inventory.append(row)
                root_counts[root_label] += 1
                root_bytes[root_label] += st.st_size
                ext_counts[(root_label, suffix)] += 1
                if name_hits:
                    candidates.append(dict(row, candidate_basis="FILENAME_OR_PATH"))

                contexts = []
                status = "NOT_SCANNED"
                if suffix in TEXT_EXTS:
                    contexts, status = scan_regular_text(path, st.st_size)
                elif suffix in OFFICE_EXTS:
                    contexts, status = scan_office_zip(path)
                if contexts:
                    for pattern, context in contexts:
                        content_hits.append({
                            "root_label": root_label,
                            "path": path_text,
                            "extension": suffix,
                            "pattern": pattern,
                            "context": context,
                            "scan_status": status,
                        })

                if suffix in ARCHIVE_EXTS and (name_hits or "project26210" in path_text.lower()):
                    members, astatus = archive_members(path)
                    for member in members:
                        mhits = pattern_hits(member)
                        if mhits:
                            archive_hits.append({
                                "root_label": root_label,
                                "archive_path": path_text,
                                "member": member,
                                "pattern_hits": ";".join(mhits),
                                "listing_status": astatus,
                            })
                    if astatus != "OK" and name_hits:
                        archive_hits.append({
                            "root_label": root_label,
                            "archive_path": path_text,
                            "member": "",
                            "pattern_hits": ";".join(name_hits),
                            "listing_status": astatus,
                        })

    inventory.sort(key=lambda r: (r["root_label"], r["path"]))
    candidates.sort(key=lambda r: (r["root_label"], r["path"]))
    content_hits.sort(key=lambda r: (r["root_label"], r["path"], r["pattern"], r["context"]))
    archive_hits.sort(key=lambda r: (r["root_label"], r["archive_path"], r["member"]))

    write_tsv(outdir / "CURRENT_FILE_INVENTORY_v1.tsv",
              ["root_label", "path", "size_bytes", "mtime_utc", "extension", "path_class", "name_pattern_hits"], inventory)
    write_tsv(outdir / "EVIDENCE_FILENAME_CANDIDATES_v1.tsv",
              ["root_label", "path", "size_bytes", "mtime_utc", "extension", "path_class", "name_pattern_hits", "candidate_basis"], candidates)
    write_tsv(outdir / "EVIDENCE_CONTENT_HITS_v1.tsv",
              ["root_label", "path", "extension", "pattern", "context", "scan_status"], content_hits)
    write_tsv(outdir / "ARCHIVE_EVIDENCE_MEMBER_HITS_v1.tsv",
              ["root_label", "archive_path", "member", "pattern_hits", "listing_status"], archive_hits)
    write_tsv(outdir / "SCAN_ERRORS_v1.tsv", ["root_label", "path", "error"], errors)

    root_summary = [{
        "root_label": label,
        "root_path": str(root),
        "file_count": root_counts[label],
        "total_bytes": root_bytes[label],
        "exists": root.exists(),
    } for label, root in roots]
    write_tsv(outdir / "ROOT_SUMMARY_v1.tsv",
              ["root_label", "root_path", "file_count", "total_bytes", "exists"], root_summary)
    ext_summary = [{"root_label": k[0], "extension": k[1], "file_count": v}
                   for k, v in sorted(ext_counts.items())]
    write_tsv(outdir / "EXTENSION_SUMMARY_v1.tsv",
              ["root_label", "extension", "file_count"], ext_summary)

    metadata = {
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "python": sys.version,
        "roots": [{"label": label, "path": str(root)} for label, root in roots],
        "counts": {
            "files": len(inventory),
            "filename_candidates": len(candidates),
            "content_hits": len(content_hits),
            "archive_member_hits": len(archive_hits),
            "errors": len(errors),
        },
        "limits": {"max_text_bytes": MAX_TEXT_BYTES, "max_zip_member_bytes": MAX_ZIP_MEMBER_BYTES},
    }
    (outdir / "SCAN_METADATA_v1.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    outputs = sorted(p for p in outdir.iterdir() if p.is_file() and p.name != "OUTPUT_SHA256_v1.tsv")
    write_tsv(outdir / "OUTPUT_SHA256_v1.tsv", ["file", "sha256", "size_bytes"],
              [{"file": p.name, "sha256": sha256(p), "size_bytes": p.stat().st_size} for p in outputs])
    print(json.dumps(metadata, ensure_ascii=False))


if __name__ == "__main__":
    main()
