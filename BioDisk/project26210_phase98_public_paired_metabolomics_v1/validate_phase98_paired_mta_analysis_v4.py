#!/usr/bin/env python3
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


MTA = "5-methylthioadenosine (MTA)"
PROGRAM = ["GPNMB", "LGALS3", "FABP5", "CD9", "LPL"]
DATASETS = {
    "RC18": ("matched_Harmonized_RawData_RC18.csv", "matched_tpm_RC18.csv"),
    "RC20": ("matched_Harmonized_RawData_RC20.csv", "matched_tpm_RC20.csv"),
    "CPTAC": ("matched_Harmonized_Met_CPTAC.csv", "matched_tpm_CPTAC.csv"),
    "CPTAC_val": ("matched_Harmonized_Met_CPTAC_val.csv", "matched_tpm_CPTAC_val.csv"),
}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def zrank(values):
    values = stats.rankdata(np.asarray(values, dtype=float), method="average")
    sd = values.std(ddof=0)
    return (values - values.mean()) / sd if sd > 0 else np.zeros_like(values)


def partial(x, y, covariates):
    x, y = zrank(x), zrank(y)
    design = np.column_stack([np.ones(len(x))] + [zrank(v) for v in covariates])
    x = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    y = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    return stats.pearsonr(x, y)


def bh(values):
    p = np.asarray(values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def add(checks, name, passed, detail):
    checks.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": str(detail)})


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: validate_phase98_paired_mta_analysis_v1.py PHASE98_ROOT")
    root = Path(sys.argv[1]).resolve()
    inp = root / "downloads_v1"
    result = root / "analysis_v2"
    out = root / "analysis_validation_v4"
    out.mkdir(parents=True, exist_ok=False)
    assoc = pd.read_csv(result / "cohort_associations_v2.tsv", sep="\t")
    sample_counts = pd.read_csv(result / "paired_sample_counts_v2.tsv", sep="\t")
    meta = pd.read_csv(result / "cross_cohort_meta_analysis_v2.tsv", sep="\t")
    med = pd.read_csv(result / "mediation_bootstrap_v2.tsv", sep="\t")
    status = json.loads((result / "analysis_status_v2.json").read_text())
    checks = []
    mapping = pd.read_csv(root / "MasterMapping_ccRCC.csv", dtype=str).fillna("")
    add(checks, "sample_count_rows", len(sample_counts) == 4, len(sample_counts))
    add(checks, "tumor_total", int(sample_counts["tumor"].sum()) == 255, int(sample_counts["tumor"].sum()))
    add(checks, "association_rows", len(assoc) == 40, len(assoc))
    add(checks, "meta_rows", len(meta) == 10, len(meta))
    add(checks, "mediation_rows", len(med) == 8, len(med))
    add(checks, "datasets_exact", set(assoc.dataset) == set(DATASETS), sorted(set(assoc.dataset)))
    add(checks, "mediation_bootstrap_count", set(med.bootstrap_replicates) == {10000}, sorted(set(med.bootstrap_replicates)))
    add(checks, "mediation_seed", set(med.seed) == {2621098}, sorted(set(med.seed)))

    recomputed = []
    sample_total = 0
    for dataset, (met_name, rna_name) in DATASETS.items():
        met = pd.read_csv(inp / met_name, index_col=0)
        rna = pd.read_csv(inp / rna_name, index_col=0)
        mapd = mapping[mapping["Dataset"] == dataset].copy()
        mapping_ok = len(mapd) == len(met) == len(rna) and set(mapd["MetabID"]) == set(met.index.astype(str)) and set(mapd["RNAID"]) == set(rna.index.astype(str))
        add(checks, f"{dataset}_mapping", mapping_ok, len(mapd))
        add(checks, f"{dataset}_mta_exact", MTA in met.columns, MTA)
        add(checks, f"{dataset}_genes", all(g in rna.columns for g in ["TYMP", "MTAP", "SPP1"] + PROGRAM), "required genes")
        met = met.loc[mapd["MetabID"]].copy()
        rna = rna.loc[mapd["RNAID"]].copy()
        paired_ids = mapd["CommonID"].where(mapd["CommonID"] != "", mapd["MetabID"]).tolist()
        met.index = paired_ids
        rna.index = paired_ids
        frame = pd.DataFrame(index=met.index)
        frame["MTA"] = pd.to_numeric(met[MTA], errors="coerce")
        frame["TN"] = mapd["TN"].to_numpy()
        for gene in ["TYMP", "MTAP", "SPP1"] + PROGRAM:
            frame[gene] = np.log2(pd.to_numeric(rna[gene], errors="coerce") + 1.0)
        frame = frame.replace([np.inf, -np.inf], np.nan).dropna()
        all_frame = frame.copy()
        frame = frame[frame["TN"] == "Tumor"].copy()
        z = (frame[PROGRAM] - frame[PROGRAM].mean()) / frame[PROGRAM].std(ddof=0)
        frame["SPP1_PROGRAM"] = z.mean(axis=1)
        sample_total += len(frame)
        specs = [
            ("A_B_TYMP_MTA_MARGINAL", "TYMP", "MTA", []),
            ("A_B_TYMP_MTA_PARTIAL_MTAP", "TYMP", "MTA", ["MTAP"]),
            ("B_C_MTA_SPP1_MARGINAL", "MTA", "SPP1", []),
            ("B_C_MTA_SPP1_PARTIAL_TYMP_MTAP", "MTA", "SPP1", ["TYMP", "MTAP"]),
            ("B_C_MTA_SPP1_PROGRAM_MARGINAL", "MTA", "SPP1_PROGRAM", []),
            ("B_C_MTA_SPP1_PROGRAM_PARTIAL_TYMP_MTAP", "MTA", "SPP1_PROGRAM", ["TYMP", "MTAP"]),
            ("A_C_TYMP_SPP1_DESCRIPTIVE", "TYMP", "SPP1", []),
            ("A_C_TYMP_SPP1_PROGRAM_DESCRIPTIVE", "TYMP", "SPP1_PROGRAM", []),
            ("MTAP_B_MTAP_MTA_CONTEXT", "MTAP", "MTA", []),
        ]
        sensitivity = all_frame[["TYMP", "MTA", "MTAP", "TN"]].dropna().copy()
        sensitivity["TN_BINARY"] = (sensitivity["TN"] == "Tumor").astype(float)
        srho, sp = partial(sensitivity["TYMP"], sensitivity["MTA"], [sensitivity["MTAP"], sensitivity["TN_BINARY"]])
        specs.append(("A_B_TYMP_MTA_ALL_PARTIAL_MTAP_TN_SENSITIVITY", "TYMP", "MTA", ["MTAP", "TN_BINARY"]))
        frame["TN_BINARY"] = 1.0
        validation_frames = {"A_B_TYMP_MTA_ALL_PARTIAL_MTAP_TN_SENSITIVITY": sensitivity}
        for endpoint, x, y, cov in specs:
            source_frame = validation_frames.get(endpoint, frame)
            values = source_frame[[x, y] + cov].dropna()
            rho, p = partial(values[x], values[y], [values[v] for v in cov]) if cov else stats.spearmanr(values[x], values[y])
            recomputed.append({"dataset": dataset, "endpoint": endpoint, "rho": float(rho), "p": float(p), "n": len(values)})
    rec = pd.DataFrame(recomputed)
    joined = assoc.merge(rec, on=["dataset", "endpoint"], suffixes=("_reported", "_recomputed"), validate="one_to_one")
    max_rho = float(np.max(np.abs(joined.rho_reported - joined.rho_recomputed)))
    max_p = float(np.max(np.abs(joined.p_reported - joined.p_recomputed)))
    add(checks, "rho_recomputed", max_rho < 1e-10, max_rho)
    add(checks, "p_recomputed", max_p < 1e-10, max_p)
    add(checks, "n_recomputed", bool((joined.n_reported == joined.n_recomputed).all()), sample_total)

    max_q = 0.0
    for family, group in assoc.groupby("family"):
        expected = bh(group.p.to_numpy())
        max_q = max(max_q, float(np.max(np.abs(expected - group.q_within_family.to_numpy()))))
    add(checks, "bh_recomputed", max_q < 1e-10, max_q)

    meta_marg = meta[meta.endpoint == "A_B_TYMP_MTA_MARGINAL"].iloc[0]
    primary = assoc[assoc.endpoint == "A_B_TYMP_MTA_MARGINAL"]
    rho = np.clip(primary.rho.to_numpy(), -0.999999, 0.999999)
    weights = primary.n.to_numpy() - 3
    zbar = float(np.sum(weights * np.arctanh(rho)) / np.sum(weights))
    fixed_rho = math.tanh(zbar)
    fixed_p = float(2 * stats.norm.sf(abs(zbar * math.sqrt(float(np.sum(weights))))))
    add(checks, "primary_meta_rho", abs(fixed_rho - meta_marg.fixed_rho) < 1e-10, fixed_rho)
    add(checks, "primary_meta_p", abs(fixed_p - meta_marg.fixed_p) < 1e-10, fixed_p)
    add(checks, "status_sample_total", int(status["total_samples"]) == int(primary.n.sum()), status["total_samples"])
    add(checks, "scope_boundary", "not macrophage fine-state" in status["scope"], status["scope"])

    files = [
        root / "PHASE98_PREREGISTRATION_v1.md",
        root / "PHASE98_PREREGISTRATION_AMENDMENT_A_v1.md",
        root / "MasterMapping_ccRCC.csv",
        root / "run_phase98_paired_mta_analysis_v2.py",
        result / "cohort_associations_v2.tsv",
        result / "paired_sample_counts_v2.tsv",
        result / "cross_cohort_meta_analysis_v2.tsv",
        result / "mediation_bootstrap_v2.tsv",
        result / "analysis_status_v2.json",
        result / "PHASE98_PAIRED_MTA_SUMMARY_v2.md",
    ]
    manifest = []
    for path in files:
        manifest.append({"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest).to_csv(out / "analysis_output_manifest_sha256_v4.tsv", sep="\t", index=False)
    pd.DataFrame(checks).to_csv(out / "validation_checks_v4.tsv", sep="\t", index=False)
    failed = [row for row in checks if row["status"] != "PASS"]
    final = {"status": "PASS" if not failed else "FAIL", "checks": len(checks), "passed": len(checks) - len(failed), "failed": len(failed)}
    (out / "validation_status_v4.json").write_text(json.dumps(final, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(final))
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
