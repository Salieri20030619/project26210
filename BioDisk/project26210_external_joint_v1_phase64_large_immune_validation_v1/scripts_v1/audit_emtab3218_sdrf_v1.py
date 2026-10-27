#!/usr/bin/env python3
import csv
import os
import sys
from collections import Counter, defaultdict


NA_VALUES = {"", "NA", "N/A", "na", "null", "None"}


def clean(value):
    value = (value or "").strip()
    return "" if value in NA_VALUES else value


def write_tsv(path, fieldnames, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main(sdrf_path, output_dir):
    with open(sdrf_path, encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    required = {
        "Source Name",
        "Characteristics[individual]",
        "Characteristics[biopsy timepoint]",
        "Characteristics[BOR]",
        "Characteristics[BOR3]",
        "Characteristics[Response20pct]",
        "Characteristics[OS]",
        "Characteristics[PFS]",
        "Characteristics[BiopsySite]",
        "Characteristics[clinical trial arm]",
        "Array Data File",
    }
    missing_required = sorted(required - set(rows[0])) if rows else sorted(required)
    if missing_required:
        raise SystemExit(f"Missing required SDRF columns: {missing_required}")

    sample_rows = []
    by_patient = defaultdict(list)
    for raw in rows:
        patient = clean(raw["Characteristics[individual]"])
        timepoint_raw = clean(raw["Characteristics[biopsy timepoint]"])
        timepoint = "baseline" if timepoint_raw == "Screen" else ("week4" if timepoint_raw == "Week4" else timepoint_raw.lower())
        bor = clean(raw["Characteristics[BOR]"]).upper()
        bor3 = clean(raw["Characteristics[BOR3]"]).upper()
        response20 = clean(raw["Characteristics[Response20pct]"])
        out = {
            "sample_id": clean(raw["Source Name"]),
            "patient_id": patient,
            "timepoint": timepoint,
            "timepoint_raw": timepoint_raw,
            "bor": bor,
            "bor3": bor3,
            "orr_response": "1" if bor in {"CR", "PR"} else ("0" if bor in {"SD", "PD"} else ""),
            "response20pct": response20,
            "os_days": clean(raw["Characteristics[OS]"]),
            "pfs_days": clean(raw["Characteristics[PFS]"]),
            "os_event": clean(raw["Characteristics[OS_SOR]"]),
            "pfs_event": clean(raw["Characteristics[PFS_SOR]"]),
            "trial_arm": clean(raw["Characteristics[clinical trial arm]"]),
            "dose": clean(raw.get("Factor Value[dose]", "")),
            "biopsy_site": clean(raw["Characteristics[BiopsySite]"]),
            "matched_lesion": clean(raw.get("Characteristics[MatchedLesionBiopsy]", "")),
            "array_file": clean(raw["Array Data File"]),
        }
        sample_rows.append(out)
        by_patient[patient].append(out)

    patient_rows = []
    inconsistent = []
    patient_fields = ["bor", "bor3", "response20pct", "os_days", "pfs_days", "os_event", "pfs_event", "trial_arm", "dose"]
    for patient, patient_samples in sorted(by_patient.items(), key=lambda item: (int(item[0]) if item[0].isdigit() else 10**9, item[0])):
        baseline = [row for row in patient_samples if row["timepoint"] == "baseline"]
        week4 = [row for row in patient_samples if row["timepoint"] == "week4"]
        representative = baseline[0] if baseline else patient_samples[0]
        for field in patient_fields:
            values = {row[field] for row in patient_samples if row[field] != ""}
            if len(values) > 1:
                inconsistent.append((patient, field, ";".join(sorted(values))))
        patient_rows.append({
            "patient_id": patient,
            "n_samples": str(len(patient_samples)),
            "n_baseline": str(len(baseline)),
            "n_week4": str(len(week4)),
            "paired_baseline_week4": "1" if baseline and week4 else "0",
            **{field: representative[field] for field in patient_fields},
            "orr_response": representative["orr_response"],
            "baseline_sample_id": baseline[0]["sample_id"] if baseline else "",
            "week4_sample_id": week4[0]["sample_id"] if week4 else "",
            "baseline_biopsy_site": baseline[0]["biopsy_site"] if baseline else "",
            "week4_biopsy_site": week4[0]["biopsy_site"] if week4 else "",
        })

    def count_values(records, field):
        return ";".join(f"{key or 'MISSING'}:{count}" for key, count in sorted(Counter(row[field] for row in records).items()))

    valid_orr = [row for row in patient_rows if row["orr_response"] in {"0", "1"}]
    gates = [
        ("sdrf_rows", len(sample_rows), 114, "PASS" if len(sample_rows) == 114 else "FAIL"),
        ("unique_patients", len(patient_rows), 59, "PASS" if len(patient_rows) >= 55 else "FAIL"),
        ("baseline_samples", sum(row["timepoint"] == "baseline" for row in sample_rows), 55, "PASS" if sum(row["timepoint"] == "baseline" for row in sample_rows) >= 55 else "FAIL"),
        ("week4_samples", sum(row["timepoint"] == "week4" for row in sample_rows), 40, "PASS" if sum(row["timepoint"] == "week4" for row in sample_rows) >= 40 else "FAIL"),
        ("paired_patients", sum(row["paired_baseline_week4"] == "1" for row in patient_rows), 40, "PASS" if sum(row["paired_baseline_week4"] == "1" for row in patient_rows) >= 40 else "FAIL"),
        ("patients_with_binary_orr", len(valid_orr), 50, "PASS" if len(valid_orr) >= 50 else "FAIL"),
        ("clinical_field_consistency", len(inconsistent), 0, "PASS" if not inconsistent else "FAIL"),
        ("unique_array_files", len({row["array_file"] for row in sample_rows}), 114, "PASS" if len({row["array_file"] for row in sample_rows}) == len(sample_rows) else "FAIL"),
    ]
    gate_rows = [
        {"gate": gate, "observed": str(observed), "threshold_or_expected": str(expected), "status": status}
        for gate, observed, expected, status in gates
    ]
    gate_rows.extend([
        {"gate": "patient_BOR_distribution", "observed": count_values(patient_rows, "bor"), "threshold_or_expected": "descriptive", "status": "INFO"},
        {"gate": "patient_BOR3_distribution", "observed": count_values(patient_rows, "bor3"), "threshold_or_expected": "descriptive", "status": "INFO"},
        {"gate": "sample_timepoint_distribution", "observed": count_values(sample_rows, "timepoint"), "threshold_or_expected": "descriptive", "status": "INFO"},
        {"gate": "patient_arm_distribution", "observed": count_values(patient_rows, "trial_arm"), "threshold_or_expected": "descriptive", "status": "INFO"},
    ])

    sample_fields = list(sample_rows[0])
    patient_out_fields = list(patient_rows[0])
    write_tsv(os.path.join(output_dir, "sample_manifest_v1.tsv"), sample_fields, sample_rows)
    write_tsv(os.path.join(output_dir, "patient_clinical_v1.tsv"), patient_out_fields, patient_rows)
    write_tsv(os.path.join(output_dir, "clinical_gate_v1.tsv"), ["gate", "observed", "threshold_or_expected", "status"], gate_rows)
    write_tsv(
        os.path.join(output_dir, "clinical_inconsistencies_v1.tsv"),
        ["patient_id", "field", "values"],
        ({"patient_id": patient, "field": field, "values": values} for patient, field, values in inconsistent),
    )

    for row in gate_rows:
        print("\t".join(row[field] for field in ("gate", "observed", "threshold_or_expected", "status")))
    if any(row["status"] == "FAIL" for row in gate_rows):
        raise SystemExit(2)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
