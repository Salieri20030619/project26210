#!/usr/bin/env python3
import csv
import hashlib
import json
import math
import sys
from pathlib import Path


MTA_NAME = "5-methylthioadenosine (MTA)"
DATASETS = {
    "RC18": ("matched_Harmonized_RawData_RC18.csv", "matched_tpm_RC18.csv"),
    "RC20": ("matched_Harmonized_RawData_RC20.csv", "matched_tpm_RC20.csv"),
    "CPTAC": ("matched_Harmonized_Met_CPTAC.csv", "matched_tpm_CPTAC.csv"),
    "CPTAC_val": ("matched_Harmonized_Met_CPTAC_val.csv", "matched_tpm_CPTAC_val.csv"),
}


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        rows = [row for row in reader if row]
    return header, [row[0] for row in rows], rows


def write_tsv(path: Path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: validate_unitedmet_download_v2.py PHASE98_ROOT")
    root = Path(sys.argv[1]).resolve()
    download_dir = root / "downloads_v1"
    record = json.loads((root / "zenodo_record_11286535_v1.json").read_text(encoding="utf-8"))
    validation_dir = root / "validation_v2"
    validation_dir.mkdir(parents=True, exist_ok=False)

    expected = {
        item["key"]: {
            "size": int(item["size"]),
            "md5": item["checksum"].split(":", 1)[1],
        }
        for item in record["files"]
    }
    file_rows = []
    all_file_pass = True
    for name, meta in sorted(expected.items()):
        path = download_dir / name
        exists = path.is_file()
        size = path.stat().st_size if exists else -1
        observed_md5 = md5_file(path) if exists else ""
        passed = exists and size == meta["size"] and observed_md5 == meta["md5"]
        all_file_pass &= passed
        file_rows.append({
            "file": name,
            "exists": str(exists).upper(),
            "expected_bytes": meta["size"],
            "observed_bytes": size,
            "expected_md5": meta["md5"],
            "observed_md5": observed_md5,
            "status": "PASS" if passed else "FAIL",
        })
    write_tsv(validation_dir / "download_file_validation_v2.tsv", file_rows)

    dataset_rows = []
    all_dataset_pass = True
    for dataset, (met_name, rna_name) in DATASETS.items():
        met_path = download_dir / met_name
        rna_path = download_dir / rna_name
        if not met_path.is_file() or not rna_path.is_file():
            dataset_rows.append({
                "dataset": dataset,
                "n_metabolomics_samples": 0,
                "n_rna_samples": 0,
                "sample_ids_ordered_exact": "FALSE",
                "mta_column_exact_count": 0,
                "mta_finite_n": 0,
                "mta_unique_n": 0,
                "TYMP_present": "FALSE",
                "MTAP_present": "FALSE",
                "SPP1_present": "FALSE",
                "status": "FAIL_MISSING_INPUT",
            })
            all_dataset_pass = False
            continue
        met_header, met_ids, met_rows = read_csv(met_path)
        rna_header, rna_ids, _ = read_csv(rna_path)
        mta_positions = [i for i, value in enumerate(met_header) if value == MTA_NAME]
        mta_values = []
        if len(mta_positions) == 1:
            position = mta_positions[0]
            for row in met_rows:
                try:
                    value = float(row[position])
                except (IndexError, TypeError, ValueError):
                    continue
                if math.isfinite(value):
                    mta_values.append(value)
        genes = {gene: gene in rna_header for gene in ("TYMP", "MTAP", "SPP1")}
        ids_match = met_ids == rna_ids
        passed = (
            len(mta_positions) == 1
            and len(mta_values) == len(met_ids)
            and len(set(mta_values)) > 1
            and all(genes.values())
            and ids_match
            and len(met_ids) == len(set(met_ids))
        )
        all_dataset_pass &= passed
        dataset_rows.append({
            "dataset": dataset,
            "n_metabolomics_samples": len(met_ids),
            "n_rna_samples": len(rna_ids),
            "sample_ids_ordered_exact": str(ids_match).upper(),
            "mta_column_exact_count": len(mta_positions),
            "mta_finite_n": len(mta_values),
            "mta_unique_n": len(set(mta_values)),
            "TYMP_present": str(genes["TYMP"]).upper(),
            "MTAP_present": str(genes["MTAP"]).upper(),
            "SPP1_present": str(genes["SPP1"]).upper(),
            "status": "PASS" if passed else "FAIL",
        })
    write_tsv(validation_dir / "paired_data_feasibility_v2.tsv", dataset_rows)

    part_files = sorted(path.name for path in download_dir.glob("*.part"))
    final_pass = all_file_pass and all_dataset_pass and not part_files
    status = {
        "status": "PASS" if final_pass else "FAIL",
        "validator_version": "v2",
        "zenodo_record": 11286535,
        "expected_files": len(expected),
        "file_validation_pass": all_file_pass,
        "paired_data_validation_pass": all_dataset_pass,
        "datasets": list(DATASETS),
        "total_paired_samples": sum(int(row["n_metabolomics_samples"]) for row in dataset_rows),
        "part_files_remaining": part_files,
        "scientific_scope": "Direct paired observational test of TYMP/MTAP RNA versus measured MTA; bulk SPP1 is not a fine-state macrophage measurement.",
    }
    (validation_dir / "validation_status_v2.json").write_text(
        json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(status, ensure_ascii=False))
    return 0 if final_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
