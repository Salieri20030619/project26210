#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from datetime import datetime
from pathlib import Path

import pandas as pd


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def parse_gse264586(metadata: Path, files_dir: Path) -> list[dict]:
    table = pd.read_csv(metadata, sep="\t", dtype=str).fillna("")
    rows: list[dict] = []
    for _, patient in table.iterrows():
        gsms = patient["geo_samples"].split(";")
        times = patient["timepoints"].split(";")
        sites = patient["biopsy_site"].split(";")
        organs = patient["organ"].split(";")
        if len(gsms) != len(times):
            raise RuntimeError(f"GSE264586 sample/timepoint mismatch: {patient['patient_id']}")
        for index, (gsm, timepoint) in enumerate(zip(gsms, times)):
            matches = sorted(files_dir.glob(f"{gsm}_*.matrix.mtx.gz"))
            if len(matches) != 1:
                raise RuntimeError(f"Expected one matrix for {gsm}, found {len(matches)}")
            stem = matches[0].name.removesuffix(".matrix.mtx.gz")
            rows.append({
                "cohort": "GSE264586", "gsm": gsm, "library_id": stem,
                "patient_id": patient["patient_id"], "timepoint": timepoint,
                "region": sites[min(index, len(sites) - 1)],
                "organ": organs[min(index, len(organs) - 1)],
                "treatment": "immunotherapy", "response_binary": patient["response_binary"],
                "response_detail": patient["best_objective_response"],
                "primary_response_eligible": timepoint == "Baseline",
                "matrix": str(matches[0]), "barcodes": str(files_dir / f"{stem}.barcodes.tsv.gz"),
                "features": str(files_dir / f"{stem}.features.tsv.gz"),
            })
    return rows


def parse_gse294109(metadata: Path, files_dir: Path) -> list[dict]:
    table = pd.read_csv(metadata, sep="\t", dtype=str).fillna("")
    rows: list[dict] = []
    for _, patient in table.iterrows():
        for gsm in patient["geo_samples"].split(";"):
            matches = sorted(files_dir.glob(f"{gsm}_*.matrix.mtx.gz"))
            if len(matches) != 1:
                raise RuntimeError(f"Expected one matrix for {gsm}, found {len(matches)}")
            stem = matches[0].name.removesuffix(".matrix.mtx.gz")
            match = re.match(r"GSM\d+_(.+)", stem)
            rows.append({
                "cohort": "GSE294109", "gsm": gsm, "library_id": stem,
                "patient_id": patient["patient_id"], "timepoint": "Pre_surgery",
                "region": match.group(1) if match else stem, "organ": patient["tissue"],
                "treatment": patient["pre_surgery_treatment"],
                "response_binary": patient["response_binary"],
                "response_detail": patient["pathological_response"],
                "primary_response_eligible": patient["response_binary"] in {"Responder", "Nonresponder"},
                "matrix": str(matches[0]), "barcodes": str(files_dir / f"{stem}.barcodes.tsv.gz"),
                "features": str(files_dir / f"{stem}.features.tsv.gz"),
            })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--download-root", type=Path, required=True)
    parser.add_argument("--metadata-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    rows = parse_gse264586(
        args.metadata_root / "GSE264586_patient_response_v1.tsv",
        args.download_root / "GSE264586" / "extracted",
    ) + parse_gse294109(
        args.metadata_root / "GSE294109_patient_metadata_v1.tsv",
        args.download_root / "GSE294109" / "files",
    )
    manifest = pd.DataFrame(rows).sort_values(["cohort", "patient_id", "gsm"]).reset_index(drop=True)
    if len(manifest) != 44:
        raise RuntimeError(f"Expected 44 libraries, found {len(manifest)}")
    if manifest[["cohort", "gsm"]].duplicated().any() or manifest["library_id"].duplicated().any():
        raise RuntimeError("Duplicate library identifier")
    for column in ["matrix", "barcodes", "features"]:
        missing = [path for path in manifest[column] if not Path(path).is_file()]
        if missing:
            raise FileNotFoundError(f"Missing {column}: {missing[:5]}")
    observed = manifest.groupby("cohort").size().to_dict()
    if observed != {"GSE264586": 19, "GSE294109": 25}:
        raise RuntimeError(f"Unexpected cohort library counts: {observed}")
    manifest.to_csv(args.output / "library_manifest_v1.tsv", sep="\t", index=False)
    summary = {
        "status": "PASS", "created_at": now(), "n_libraries": int(len(manifest)),
        "n_patients": int(manifest[["cohort", "patient_id"]].drop_duplicates().shape[0]),
        "cohort_libraries": observed,
        "cohort_patients": manifest.groupby("cohort")["patient_id"].nunique().to_dict(),
        "primary_response_patients": manifest.loc[manifest["primary_response_eligible"]]
        .groupby("cohort")["patient_id"].nunique().to_dict(),
    }
    (args.output / "library_manifest_summary_v1.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
