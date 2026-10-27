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


def read_csv_header_and_ids(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        ids = []
        rows = []
        for row in reader:
            if not row:
                continue
            ids.append(row[0])
            rows.append(row)
    return header, ids, rows


def write_tsv(path: Path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: validate_unitedmet_download_v1.py PHASE98_ROOT")
    root = Path(sys.argv[1]).resolve()
    download_dir = root / "downloads_v1"
    record_path = root / "zenodo_record_11286535_v1.json"
    validation_dir = root / "validation_v1"
    validation_dir.mkdir(parents=True, exist_ok=False)

    record = json.loads(record_path.read_text(encoding="utf-8"))
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
        file_rows.append(
            {
                "file": name,
                "exists": str(exists).upper(),
                "expected_bytes": meta["size"],
                "observed_bytes": size,
                "expected_md5": meta["md5"],
                "observed_md5": observed_md5,
                "status": "PASS" if passed else "FAIL",
            }
        )
    write_tsv(validation_dir / "download_file_validation_v1.tsv", file_rows)

    dataset_rows = []
    all_dataset_pass = True
    for dataset, (met_name, rna_name) in DATASETS.items():
        met_header, met_ids, met_rows = read_csv_header_and_ids(download_dir / met_name)
        rna_header, rna_ids, _ = read_csv_header_and_ids(download_dir / rna_name)
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
        genes_present = {gene: gene in rna_header for gene in ("TYMP", "MTAP", "SPP1")}
        ids_match = met_ids == rna_ids
        passed = (
            len(mta_positions) == 1
            and len(mta_values) == len(met_ids)
            and len(set(mta_values)) > 1
            and all(genes_present.values())
            and ids_match
            and len(met_ids) == len(set(met_ids))
        )
        all_dataset_pass &= passed
        dataset_rows.append(
            {
                "dataset": dataset,
                "n_metabolomics_samples": len(met_ids),
                "n_rna_samples": len(rna_ids),
                "sample_ids_ordered_exact": str(ids_match).upper(),
                "mta_column_exact_count": len(mta_positions),
                "mta_finite_n": len(mta_values),
                "mta_unique_n": len(set(mta_values)),
                "TYMP_present": str(genes_present["TYMP"]).upper(),
                "MTAP_present": str(genes_present["MTAP"]).upper(),
                "SPP1_present": str(genes_present["SPP1"]).upper(),
                "status": "PASS" if passed else "FAIL",
            }
        )
    write_tsv(validation_dir / "paired_data_feasibility_v1.tsv", dataset_rows)

    part_files = sorted(path.name for path in download_dir.glob("*.part"))
    status = {
        "status": "PASS" if all_file_pass and all_dataset_pass and not part_files else "FAIL",
        "zenodo_record": 11286535,
        "expected_files": len(expected),
        "file_validation_pass": all_file_pass,
        "paired_data_validation_pass": all_dataset_pass,
        "datasets": list(DATASETS),
        "total_paired_samples": sum(int(row["n_metabolomics_samples"]) for row in dataset_rows),
        "part_files_remaining": part_files,
        "scientific_scope": "Direct paired observational test of TYMP/MTAP RNA versus measured MTA; bulk SPP1 is not a fine-state macrophage measurement.",
    }
    (validation_dir / "validation_status_v1.json").write_text(
        json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(status, ensure_ascii=False))
    return 0 if status["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
