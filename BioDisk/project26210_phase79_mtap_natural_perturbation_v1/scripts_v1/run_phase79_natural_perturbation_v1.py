#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

SPP1_STATE = ["GPNMB", "LGALS3", "FABP5", "CD9", "LPL"]
MTA_RESPONSE = ["ARG1", "VEGFA", "IL10", "TIMP1"]
IMMUNOSUPPRESSION = ["TGFB1", "VSIG4", "MARCO", "CD274", "IDO1"]
C1QC_STATE = ["C1QA", "C1QB", "C1QC"]
MACROPHAGE = ["CSF1R", "FCER1G", "TYROBP", "AIF1", "LST1"]
PROGRAMS = {
    "spp1_state": SPP1_STATE,
    "mta_response": MTA_RESPONSE,
    "immunosuppression": IMMUNOSUPPRESSION,
    "c1qc_specificity_control": C1QC_STATE,
    "macrophage_covariate": MACROPHAGE,
}
PRIMARY_OUTCOMES = ["spp1_state", "mta_response", "immunosuppression"]


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def scan_rows(path: Path, genes: set[str]) -> tuple[list[str], dict[str, np.ndarray]]:
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        samples = header[2:]
        rows = {}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            gene = fields[0].upper()
            if gene in genes and gene not in rows:
                rows[gene] = np.asarray(
                    [float(x) if x not in {"", "NA"} else np.nan for x in fields[2:]],
                    dtype=float,
                )
    return samples, rows


def zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    sd = np.nanstd(values, ddof=1)
    if np.isfinite(sd) and sd > 0:
        return (values - np.nanmean(values)) / sd
    return np.zeros_like(values)


def program_score(expression: dict[str, np.ndarray], genes: list[str]) -> np.ndarray:
    missing = sorted(set(genes).difference(expression))
    if missing:
        raise RuntimeError(f"Missing programme genes: {missing}")
    return np.nanmean(
        np.vstack([zscore(np.log2(expression[gene] + 1.0)) for gene in genes]),
        axis=0,
    )


def hc3_fit(y: np.ndarray, frame: pd.DataFrame, terms: list[str]) -> dict[str, dict[str, float]]:
    columns = [np.ones(len(frame))] + [frame[term].to_numpy(float) for term in terms]
    x = np.column_stack(columns)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(y) & np.isfinite(x).all(axis=1)
    x = x[ok]
    y = y[ok]
    names = ["intercept"] + terms
    if len(y) <= x.shape[1] + 2 or np.linalg.matrix_rank(x) < x.shape[1]:
        return {
            name: {"beta": np.nan, "se": np.nan, "p": np.nan, "n": int(len(y))}
            for name in names
        }
    xtx_inv = np.linalg.inv(x.T @ x)
    beta = xtx_inv @ x.T @ y
    residual = y - x @ beta
    leverage = np.sum((x @ xtx_inv) * x, axis=1)
    scaled = residual / np.maximum(1.0 - leverage, 1e-8)
    meat = x.T @ ((scaled ** 2)[:, None] * x)
    covariance = xtx_inv @ meat @ xtx_inv
    se = np.sqrt(np.maximum(np.diag(covariance), 0))
    tvalue = beta / se
    dof = len(y) - x.shape[1]
    p = 2.0 * stats.t.sf(np.abs(tvalue), dof)
    return {
        name: {
            "beta": float(beta[index]),
            "se": float(se[index]),
            "p": float(p[index]),
            "n": int(len(y)),
        }
        for index, name in enumerate(names)
    }


def bh(values: pd.Series) -> pd.Series:
    p = values.to_numpy(float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * len(pv) / np.arange(1, len(pv) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(pv)
    adjusted[order] = np.minimum(ranked, 1.0)
    out[np.flatnonzero(ok)] = adjusted
    return pd.Series(out, index=values.index)


def build_tcga(expression_path: Path, cna_path: Path, log2_cna_path: Path) -> pd.DataFrame:
    needed = {"TYMP", "MTAP"}
    for genes in PROGRAMS.values():
        needed.update(genes)
    samples, expression = scan_rows(expression_path, needed)
    cna_samples, cna = scan_rows(cna_path, {"MTAP"})
    log_samples, log_cna = scan_rows(log2_cna_path, {"MTAP"})
    common_samples = set(samples).intersection(cna_samples).intersection(log_samples)
    if len(common_samples) < 400:
        raise RuntimeError(f"Too few shared TCGA expression/CNA samples: {len(common_samples)}")
    cna_index = {sample: index for index, sample in enumerate(cna_samples)}
    log_index = {sample: index for index, sample in enumerate(log_samples)}
    selected = []
    seen = set()
    for index, sample in enumerate(samples):
        if sample not in common_samples:
            continue
        patient = sample[:12]
        if sample[13:15] != "01" or patient in seen:
            continue
        seen.add(patient)
        selected.append((index, sample, patient))
    scores = {
        name: program_score(expression, genes)
        for name, genes in PROGRAMS.items()
    }
    tymp = zscore(np.log2(expression["TYMP"] + 1.0))
    mtap_expr = zscore(np.log2(expression["MTAP"] + 1.0))
    rows = []
    for index, sample, patient in selected:
        cna_index_value = cna_index[sample]
        log_index_value = log_index[sample]
        row = {
            "patient": patient,
            "sample": sample,
            "tymp_z": tymp[index],
            "mtap_expr_z": mtap_expr[index],
            "mtap_cna": cna["MTAP"][cna_index_value],
            "mtap_log2_cna": log_cna["MTAP"][log_index_value],
        }
        for name, values in scores.items():
            row[name] = values[index]
        rows.append(row)
    frame = pd.DataFrame(rows)
    frame["mtap_loss"] = (frame["mtap_cna"] <= -1).astype(float)
    frame["neg_mtap_log2_cna_z"] = zscore(-frame["mtap_log2_cna"].to_numpy())
    frame["interaction"] = frame["tymp_z"] * frame["mtap_loss"]
    frame["interaction_continuous"] = frame["tymp_z"] * frame["neg_mtap_log2_cna_z"]
    return frame


def tcga_models(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    primary_frame = frame.loc[frame["mtap_cna"].isin([-2, -1, 0])].copy()
    rows = []
    for outcome in PRIMARY_OUTCOMES + ["c1qc_specificity_control"]:
        fit = hc3_fit(
            primary_frame[outcome].to_numpy(),
            primary_frame,
            ["tymp_z", "mtap_loss", "interaction", "macrophage_covariate"],
        )
        rows.append({
            "layer": "TCGA_KIRC_CNA",
            "outcome": outcome,
            "term": "TYMP_x_MTAP_loss",
            **fit["interaction"],
            "loss_n": int(primary_frame["mtap_loss"].sum()),
            "neutral_n": int((primary_frame["mtap_loss"] == 0).sum()),
        })
    result = pd.DataFrame(rows)
    mask = result["outcome"].isin(PRIMARY_OUTCOMES)
    result.loc[mask, "q"] = bh(result.loc[mask, "p"])

    sensitivity = []
    for outcome in PRIMARY_OUTCOMES:
        fit = hc3_fit(
            frame[outcome].to_numpy(),
            frame,
            ["tymp_z", "neg_mtap_log2_cna_z", "interaction_continuous", "macrophage_covariate"],
        )
        sensitivity.append({
            "layer": "TCGA_KIRC_CONTINUOUS_CNA",
            "outcome": outcome,
            "term": "TYMP_x_negative_MTAP_log2CNA",
            **fit["interaction_continuous"],
        })
        fit_mrna = hc3_fit(
            primary_frame[outcome].to_numpy(),
            primary_frame,
            ["tymp_z", "mtap_loss", "interaction", "macrophage_covariate", "mtap_expr_z"],
        )
        sensitivity.append({
            "layer": "TCGA_KIRC_CNA_ADJUST_MTAP_MRNA",
            "outcome": outcome,
            "term": "TYMP_x_MTAP_loss",
            **fit_mrna["interaction"],
        })
    return result, pd.DataFrame(sensitivity)


def molecular_gate(frame: pd.DataFrame) -> dict[str, float]:
    subset = frame.loc[frame["mtap_cna"].isin([-2, -1, 0])].copy()
    loss = subset.loc[subset["mtap_loss"] == 1, "mtap_expr_z"].to_numpy()
    neutral = subset.loc[subset["mtap_loss"] == 0, "mtap_expr_z"].to_numpy()
    test = stats.mannwhitneyu(loss, neutral, alternative="less")
    return {
        "loss_n": int(len(loss)),
        "neutral_n": int(len(neutral)),
        "loss_median_mtap_expr_z": float(np.median(loss)),
        "neutral_median_mtap_expr_z": float(np.median(neutral)),
        "one_sided_p": float(test.pvalue),
        "pass": bool(np.median(loss) < np.median(neutral) and test.pvalue < 0.05),
    }


def external_models(paths: list[Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    frames = []
    for path in paths:
        frame = pd.read_csv(path, sep="\t")
        frame["dataset"] = frame["dataset"].astype(str)
        frames.append(frame)
    all_data = pd.concat(frames, ignore_index=True)
    formal = all_data.groupby("dataset").filter(lambda x: len(x) >= 8).copy()
    for column, target in [
        ("tumor_tymp_mean", "tymp_z"),
        ("tumor_mtap_mean", "mtap_z"),
    ]:
        formal[target] = formal.groupby("dataset")[column].transform(
            lambda x: zscore(x.to_numpy())
        )
    formal["low_mtap_z"] = -formal["mtap_z"]
    formal["interaction"] = formal["tymp_z"] * formal["low_mtap_z"]
    formal["log_n_tumor"] = np.log1p(formal["n_tumor_cells"].to_numpy())
    formal["log_n_myeloid"] = np.log1p(formal["n_myeloid_cells"].to_numpy())
    datasets = sorted(formal["dataset"].unique())
    for dataset in datasets[1:]:
        formal[f"cohort_{dataset}"] = (formal["dataset"] == dataset).astype(float)
    dummies = [f"cohort_{dataset}" for dataset in datasets[1:]]
    outcomes = {
        "spp1_state": "myeloid_spp1_state_mean",
        "mta_response": "myeloid_mta_response_mean",
        "immunosuppression": "myeloid_immunosuppression_mean",
    }
    rows = []
    loo_rows = []
    terms = dummies + [
        "tymp_z", "low_mtap_z", "interaction", "log_n_tumor", "log_n_myeloid"
    ]
    for outcome, column in outcomes.items():
        fit = hc3_fit(formal[column].to_numpy(), formal, terms)
        rows.append({
            "layer": "EXTERNAL_SCRNA_POOLED",
            "outcome": outcome,
            "term": "TYMP_x_low_MTAP_expression",
            **fit["interaction"],
            "formal_cohorts": ";".join(datasets),
        })
        for omitted in datasets:
            subset = formal.loc[formal["dataset"] != omitted].copy()
            reduced_dummies = [
                term for term in dummies
                if term in subset and subset[term].nunique() > 1
            ]
            reduced_terms = reduced_dummies + [
                "tymp_z", "low_mtap_z", "interaction",
                "log_n_tumor", "log_n_myeloid",
            ]
            loo_fit = hc3_fit(subset[column].to_numpy(), subset, reduced_terms)
            loo_rows.append({
                "outcome": outcome,
                "omitted_cohort": omitted,
                **loo_fit["interaction"],
            })
    result = pd.DataFrame(rows)
    result["q"] = bh(result["p"])
    return result, pd.DataFrame(loo_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression", required=True, type=Path)
    parser.add_argument("--cna", required=True, type=Path)
    parser.add_argument("--log2-cna", required=True, type=Path)
    parser.add_argument("--patient-tables", nargs="+", required=True, type=Path)
    parser.add_argument("--preregistration", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.mkdir(parents=True)

    tcga = build_tcga(args.expression, args.cna, args.log2_cna)
    gate = molecular_gate(tcga)
    tcga_primary, tcga_sensitivity = tcga_models(tcga)
    external, external_loo = external_models(args.patient_tables)

    tcga.to_csv(
        args.output / "tcga_kirc_natural_perturbation_cohort_v1.tsv",
        sep="\t", index=False,
    )
    tcga_primary.to_csv(
        args.output / "tcga_kirc_interaction_primary_v1.tsv",
        sep="\t", index=False,
    )
    tcga_sensitivity.to_csv(
        args.output / "tcga_kirc_interaction_sensitivity_v1.tsv",
        sep="\t", index=False,
    )
    external.to_csv(
        args.output / "external_scrna_interaction_primary_v1.tsv",
        sep="\t", index=False,
    )
    external_loo.to_csv(
        args.output / "external_scrna_leave_one_cohort_out_v1.tsv",
        sep="\t", index=False,
    )
    (args.output / "mtap_molecular_validation_gate_v1.json").write_text(
        json.dumps(gate, indent=2) + "\n", encoding="utf-8"
    )

    tcga_spp1 = tcga_primary.loc[
        tcga_primary["outcome"] == "spp1_state"
    ].iloc[0]
    ext_spp1 = external.loc[external["outcome"] == "spp1_state"].iloc[0]
    loo_spp1 = external_loo.loc[external_loo["outcome"] == "spp1_state"]
    genomic_support = bool(
        gate["pass"] and tcga_spp1["beta"] > 0 and tcga_spp1["q"] < 0.05
    )
    external_support = bool(
        ext_spp1["beta"] > 0
        and ext_spp1["q"] < 0.05
        and len(loo_spp1) > 0
        and (loo_spp1["beta"] > 0).all()
    )
    decision = {
        "status": "PASS",
        "completed_at": now(),
        "genomic_support": genomic_support,
        "external_scrna_support": external_support,
        "full_natural_perturbation_support": bool(
            genomic_support and external_support
        ),
        "preregistration": str(args.preregistration),
        "preregistration_sha256": sha256(args.preregistration),
        "inputs": {
            "expression_sha256": sha256(args.expression),
            "cna_sha256": sha256(args.cna),
            "log2_cna_sha256": sha256(args.log2_cna),
            "patient_table_sha256": {
                str(path): sha256(path) for path in args.patient_tables
            },
        },
    }
    (args.output / "PHASE79_DECISION_v1.json").write_text(
        json.dumps(decision, indent=2) + "\n", encoding="utf-8"
    )
    (args.output / "PHASE79_PASS_v1").write_text(
        "status=PASS\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
