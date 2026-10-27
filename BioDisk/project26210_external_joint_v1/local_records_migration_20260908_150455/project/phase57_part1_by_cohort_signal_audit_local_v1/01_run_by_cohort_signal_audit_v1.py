from __future__ import annotations

import hashlib
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(r"D:\project26210_external_joint_v1")
PHASE = Path(__file__).resolve().parent
OUT = PHASE / "analysis_local_v1"

EXPECTED_STUDIES = [
    "GSE156632",
    "GSE159115",
    "GSE178481",
    "GSE207493",
    "GSE242299",
]
DOMINANT_STUDY = "GSE207493"
RNG_SEED = 2621057


def study_from_patient(patient_id: str) -> str:
    return str(patient_id).split("__", 1)[0]


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def make_frame(
    source: pd.DataFrame,
    signal_id: str,
    signal_label: str,
    evidence_role: str,
    effect_col: str,
    source_file: Path,
    study_col: str | None = None,
    patient_col: str = "curated_patient_id",
    n_cells_col: str = "n_cells",
) -> pd.DataFrame:
    out = pd.DataFrame(
        {
            "signal_id": signal_id,
            "signal_label": signal_label,
            "evidence_role": evidence_role,
            "study": (
                source[study_col].astype(str)
                if study_col
                else source[patient_col].map(study_from_patient)
            ),
            "curated_patient_id": source[patient_col].astype(str),
            "n_cells": pd.to_numeric(source[n_cells_col], errors="coerce"),
            "effect": pd.to_numeric(source[effect_col], errors="coerce"),
            "source_file": str(source_file),
        }
    )
    return out[np.isfinite(out["effect"])].copy()


def load_patient_effects() -> pd.DataFrame:
    frames: list[pd.DataFrame] = []

    p31 = (
        ROOT
        / "phase31_part1_m17_enzyme_matrix_local_v1"
        / "results_local_v1"
        / "m17_patient_gene_f1_effects_v1.tsv"
    )
    d31 = pd.read_csv(p31, sep="\t")
    d31 = d31[d31["endpoint"].eq("F1_continuous")]
    gene_specs = {
        "TYMP": ("F1__TYMP_RNA", "F1–TYMP RNA", "core"),
        "ODC1": ("F1__ODC1_RNA", "F1–ODC1 RNA", "secondary"),
        "TK1": ("F1__TK1_RNA", "F1–TK1 RNA", "exploratory"),
        "SRM": ("F1__SRM_RNA", "F1–SRM RNA", "exploratory"),
        "PCYT2": ("F1__PCYT2_RNA", "F1–PCYT2 RNA", "exploratory"),
    }
    for gene, spec in gene_specs.items():
        frames.append(make_frame(d31[d31["gene"].eq(gene)], *spec, "effect", p31))

    p30 = (
        ROOT
        / "phase30_part1_enzyme_function_axis_local_v1"
        / "results_local_v1"
        / "enzyme_axis_patient_effects_v1.tsv"
    )
    d30 = pd.read_csv(p30, sep="\t")
    d30 = d30[d30["endpoint"].eq("F1_continuous")]
    axis_specs = {
        "TYMP_axis_primary": (
            "F1__TYMP_axis_strict",
            "F1–strict TYMP functional axis",
            "core_sensitivity",
        ),
        "TYMP_axis_broad_sensitivity": (
            "F1__TYMP_axis_broad",
            "F1–broad TYMP/pyrimidine axis",
            "secondary",
        ),
        "RRM2B_axis_primary": (
            "F1__RRM2B_axis_strict",
            "F1–strict RRM2B functional axis",
            "core_sensitivity",
        ),
        "RRM2B_axis_broad_sensitivity": (
            "F1__RRM_axis_broad",
            "F1–broad RNR-linked axis",
            "opposing_sensitivity",
        ),
    }
    for feature, spec in axis_specs.items():
        frames.append(make_frame(d30[d30["feature"].eq(feature)], *spec, "effect", p30))

    p28 = (
        ROOT
        / "phase28_part1_hif1a_sensitivity_local_v1"
        / "server_results_sync_v1"
        / "analysis_v2"
        / "patient_hif1a_sensitivity_effects_v1.tsv"
    )
    d28 = pd.read_csv(p28, sep="\t")
    d28 = d28[d28["endpoint"].eq("continuous_F1_spearman")]
    modality_specs = {
        "frozen_downstream_5gene": (
            "F1__hypoxia_downstream_5gene",
            "F1–hypoxia/angiogenesis downstream program",
            "core",
        ),
        "metaVIPER_HIF1A": (
            "F1__metaVIPER_HIF1A",
            "F1–metaVIPER HIF1A activity",
            "independent_contrast",
        ),
    }
    for modality, spec in modality_specs.items():
        frames.append(
            make_frame(
                d28[d28["modality"].eq(modality)],
                *spec,
                "effect",
                p28,
                study_col="study",
            )
        )

    p25 = (
        ROOT
        / "phase25_part1_scdrs_rebuild_local_v1"
        / "server_results_sync_v1"
        / "analysis_v1"
        / "patient_within_association_v1.tsv"
    )
    d25 = pd.read_csv(p25, sep="\t")
    d25 = d25[as_bool(d25["eligible_primary"])]
    frames.append(
        make_frame(
            d25,
            "F1__scDRS_RCC_score",
            "F1–scDRS RCC susceptibility score",
            "independent",
            "spearman_rho_f1_vs_scdrs",
            p25,
            study_col="study",
        )
    )

    p26 = (
        ROOT
        / "phase26_part1_cnmf_pyscenic_local_v1"
        / "server_results_sync_v2"
        / "analysis_v3"
        / "patient_continuous_regulon_effects_v1.tsv.gz"
    )
    d26 = pd.read_csv(p26, sep="\t")
    d26 = d26[as_bool(d26["eligible"])]
    for tf in ["ILF2", "HOXD3", "HNF4A"]:
        for method in ["GENIE3", "GRNBoost2"]:
            subset = d26[d26["tf"].eq(tf) & d26["method"].eq(method)]
            frames.append(
                make_frame(
                    subset,
                    f"F1__{tf}_{method}",
                    f"F1–{tf} regulon ({method})",
                    "regulator_sensitivity",
                    "spearman_rho",
                    p26,
                    study_col="study",
                )
            )

    out = pd.concat(frames, ignore_index=True)
    out = out[out["study"].isin(EXPECTED_STUDIES)].copy()
    out = out.sort_values(["signal_id", "study", "curated_patient_id"])
    duplicated = out.duplicated(["signal_id", "curated_patient_id"], keep=False)
    if duplicated.any():
        raise RuntimeError("Duplicate patient effects detected:\n" + out[duplicated].to_string())
    return out


def wilcoxon_p(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    values = values[values != 0]
    if len(values) == 0:
        return 1.0
    try:
        return float(stats.wilcoxon(values, alternative="two-sided", method="auto").pvalue)
    except ValueError:
        return float("nan")


def bh_adjust(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.notna() & np.isfinite(pd.to_numeric(values, errors="coerce"))
    p = values[valid].astype(float)
    if p.empty:
        return result
    order = np.argsort(p.to_numpy())
    ranked = p.to_numpy()[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    ordered_index = p.index.to_numpy()[order]
    result.loc[ordered_index] = adjusted
    return result


def bootstrap_median_ci(values: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    if len(values) < 3:
        return float("nan"), float("nan")
    draw = rng.choice(values, size=(10000, len(values)), replace=True)
    medians = np.median(draw, axis=1)
    return tuple(float(x) for x in np.quantile(medians, [0.025, 0.975]))


def summarize_group(group: pd.DataFrame, rng: np.random.Generator) -> dict[str, float | int]:
    values = group["effect"].to_numpy(float)
    low, high = bootstrap_median_ci(values, rng)
    return {
        "n_patients": int(len(values)),
        "n_positive": int(np.sum(values > 0)),
        "n_negative": int(np.sum(values < 0)),
        "n_zero": int(np.sum(values == 0)),
        "median_effect": float(np.median(values)),
        "mean_effect": float(np.mean(values)),
        "bootstrap_median_ci_low": low,
        "bootstrap_median_ci_high": high,
        "wilcoxon_p_two_sided": wilcoxon_p(values),
        "total_cells": int(group["n_cells"].fillna(0).sum()),
    }


def random_effects_meta(group: pd.DataFrame) -> dict[str, float | int]:
    rows: list[tuple[float, float]] = []
    for _, cohort in group.groupby("study", sort=False):
        vals = np.clip(cohort["effect"].to_numpy(float), -0.999999, 0.999999)
        if len(vals) < 2:
            continue
        zvals = np.arctanh(vals)
        variance = float(np.var(zvals, ddof=1) / len(zvals))
        if not np.isfinite(variance):
            continue
        rows.append((float(np.mean(zvals)), max(variance, 1e-8)))
    if len(rows) < 2:
        return {
            "meta_k_estimable_studies": len(rows),
            "meta_random_rho": np.nan,
            "meta_ci_low": np.nan,
            "meta_ci_high": np.nan,
            "meta_p_two_sided": np.nan,
            "tau2": np.nan,
            "I2_percent": np.nan,
            "Q": np.nan,
            "Q_p": np.nan,
        }
    yi = np.array([x[0] for x in rows])
    vi = np.array([x[1] for x in rows])
    wi = 1.0 / vi
    fixed = float(np.sum(wi * yi) / np.sum(wi))
    q = float(np.sum(wi * (yi - fixed) ** 2))
    df = len(yi) - 1
    c = float(np.sum(wi) - np.sum(wi**2) / np.sum(wi))
    tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0
    wr = 1.0 / (vi + tau2)
    pooled = float(np.sum(wr * yi) / np.sum(wr))
    se = math.sqrt(1.0 / float(np.sum(wr)))
    p = float(2 * stats.norm.sf(abs(pooled / se))) if se > 0 else np.nan
    q_p = float(stats.chi2.sf(q, df))
    i2 = max(0.0, (q - df) / q) * 100 if q > 0 else 0.0
    return {
        "meta_k_estimable_studies": len(rows),
        "meta_random_rho": float(np.tanh(pooled)),
        "meta_ci_low": float(np.tanh(pooled - 1.96 * se)),
        "meta_ci_high": float(np.tanh(pooled + 1.96 * se)),
        "meta_p_two_sided": p,
        "tau2": tau2,
        "I2_percent": i2,
        "Q": q,
        "Q_p": q_p,
    }


def classify_signal(
    overall: pd.Series,
    per_cohort: pd.DataFrame,
    loo: pd.DataFrame,
    dominant_vs_others: pd.Series,
) -> tuple[str, str]:
    overall_effect = float(overall["median_effect"])
    overall_sign = int(np.sign(overall_effect))
    cohort_effects = per_cohort["median_effect"].to_numpy(float)
    cohort_signs = np.sign(cohort_effects[cohort_effects != 0])
    same_study_fraction = (
        float(np.mean(cohort_signs == overall_sign)) if len(cohort_signs) else np.nan
    )
    without_dom = loo[loo["left_out_study"].eq(DOMINANT_STUDY)]
    without_dom_effect = (
        float(without_dom.iloc[0]["median_effect"]) if len(without_dom) else np.nan
    )
    dom_effect = float(dominant_vs_others["dominant_median_effect"])
    other_effect = float(dominant_vs_others["other_median_effect"])
    dom_share = float(dominant_vs_others["dominant_patient_fraction"])
    overall_p = float(overall["wilcoxon_fdr_bh_audit18"])

    if (
        np.isfinite(other_effect)
        and np.isfinite(dom_effect)
        and abs(other_effect) >= 0.01
        and abs(dom_effect) >= 0.01
        and np.sign(other_effect) != 0
        and np.sign(dom_effect) != np.sign(other_effect)
    ):
        return (
            "COHORT_HETEROGENEITY_OPPOSITE",
            "GSE207493 and the other cohorts have opposite median directions.",
        )
    if (
        np.isfinite(without_dom_effect)
        and overall_sign != 0
        and overall_p < 0.05
        and np.sign(without_dom_effect) == overall_sign
        and abs(without_dom_effect) >= 0.4 * abs(overall_effect)
        and same_study_fraction >= 0.75
    ):
        return (
            "CROSS_COHORT_DIRECTIONAL_SUPPORT",
            "The pooled patient-level test is supported; removing GSE207493 preserves direction and at least 40% of the pooled median, and at least 75% of cohort medians agree.",
        )
    if (
        np.isfinite(without_dom_effect)
        and overall_sign != 0
        and overall_p < 0.05
        and np.sign(without_dom_effect) == overall_sign
        and same_study_fraction >= 0.75
        and abs(without_dom_effect) < 0.4 * abs(overall_effect)
    ):
        return (
            "GSE207493_WEIGHTED_DIRECTION_REPLICATED",
            "Other cohorts preserve the direction, but the median magnitude falls by more than 60% after removing GSE207493.",
        )
    if (
        np.isfinite(without_dom_effect)
        and overall_sign != 0
        and overall_p < 0.05
        and (
            np.sign(without_dom_effect) != overall_sign
            or abs(without_dom_effect) < 0.25 * abs(overall_effect)
        )
        and dom_share >= 0.40
    ):
        return (
            "GSE207493_DEPENDENT",
            "The pooled direction flips or loses at least 75% of its median magnitude after removing GSE207493.",
        )
    if same_study_fraction >= 0.75 and overall_p >= 0.05:
        return (
            "DIRECTIONALLY_CONSISTENT_UNDERPOWERED",
            "Most cohort medians agree, but the pooled patient-level test is not supported.",
        )
    return (
        "MIXED_OR_UNDERPOWERED",
        "Cohort directions are mixed and available cohort-level sample sizes do not support a stronger classification.",
    )


def fmt(value: float, digits: int = 3) -> str:
    if value is None or not np.isfinite(value):
        return "NA"
    return f"{value:.{digits}f}"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)
    effects = load_patient_effects()
    effects.to_csv(OUT / "patient_effects_long_v1.tsv", sep="\t", index=False)

    cohort_rows: list[dict] = []
    for (signal_id, study), group in effects.groupby(["signal_id", "study"], sort=True):
        first = group.iloc[0]
        cohort_rows.append(
            {
                "signal_id": signal_id,
                "signal_label": first["signal_label"],
                "evidence_role": first["evidence_role"],
                "study": study,
                **summarize_group(group, rng),
            }
        )
    per_cohort = pd.DataFrame(cohort_rows)
    per_cohort["wilcoxon_fdr_bh_within_study"] = per_cohort.groupby(
        "study", group_keys=False
    )["wilcoxon_p_two_sided"].apply(bh_adjust)
    per_cohort.to_csv(OUT / "per_cohort_effects_v1.tsv", sep="\t", index=False)

    loo_rows: list[dict] = []
    synthesis_rows: list[dict] = []
    comparison_rows: list[dict] = []
    for signal_id, group in effects.groupby("signal_id", sort=True):
        first = group.iloc[0]
        overall_stats = summarize_group(group, rng)
        cohort_subset = per_cohort[per_cohort["signal_id"].eq(signal_id)]
        study_medians = cohort_subset["median_effect"].to_numpy(float)
        meta = random_effects_meta(group)
        synthesis_rows.append(
            {
                "signal_id": signal_id,
                "signal_label": first["signal_label"],
                "evidence_role": first["evidence_role"],
                "n_evaluable_studies": int(group["study"].nunique()),
                "n_positive_study_medians": int(np.sum(study_medians > 0)),
                "n_negative_study_medians": int(np.sum(study_medians < 0)),
                "n_zero_study_medians": int(np.sum(study_medians == 0)),
                "unweighted_median_of_study_medians": float(np.median(study_medians)),
                **overall_stats,
                **meta,
            }
        )
        for left_out in EXPECTED_STUDIES:
            kept = group[~group["study"].eq(left_out)]
            if kept.empty:
                continue
            loo_rows.append(
                {
                    "signal_id": signal_id,
                    "signal_label": first["signal_label"],
                    "left_out_study": left_out,
                    "n_remaining_studies": int(kept["study"].nunique()),
                    **summarize_group(kept, rng),
                }
            )
        dom = group[group["study"].eq(DOMINANT_STUDY)]
        other = group[~group["study"].eq(DOMINANT_STUDY)]
        dom_values = dom["effect"].to_numpy(float)
        other_values = other["effect"].to_numpy(float)
        mw_p = (
            float(stats.mannwhitneyu(dom_values, other_values, alternative="two-sided").pvalue)
            if len(dom_values) >= 2 and len(other_values) >= 2
            else np.nan
        )
        comparison_rows.append(
            {
                "signal_id": signal_id,
                "signal_label": first["signal_label"],
                "dominant_study": DOMINANT_STUDY,
                "dominant_n_patients": len(dom_values),
                "other_n_patients": len(other_values),
                "dominant_patient_fraction": len(dom_values) / len(group),
                "dominant_median_effect": float(np.median(dom_values)) if len(dom_values) else np.nan,
                "other_median_effect": float(np.median(other_values)) if len(other_values) else np.nan,
                "dominant_minus_other_median": (
                    float(np.median(dom_values) - np.median(other_values))
                    if len(dom_values) and len(other_values)
                    else np.nan
                ),
                "mannwhitney_p_two_sided": mw_p,
            }
        )

    synthesis = pd.DataFrame(synthesis_rows)
    synthesis["wilcoxon_fdr_bh_audit18"] = bh_adjust(synthesis["wilcoxon_p_two_sided"])
    synthesis["meta_fdr_bh_audit18"] = bh_adjust(synthesis["meta_p_two_sided"])
    loo = pd.DataFrame(loo_rows)
    comparison = pd.DataFrame(comparison_rows)
    loo.to_csv(OUT / "leave_one_dataset_out_v1.tsv", sep="\t", index=False)
    comparison.to_csv(OUT / "gse207493_vs_other_cohorts_v1.tsv", sep="\t", index=False)

    classifications: list[dict] = []
    for _, overall in synthesis.iterrows():
        signal_id = overall["signal_id"]
        category, rationale = classify_signal(
            overall,
            per_cohort[per_cohort["signal_id"].eq(signal_id)],
            loo[loo["signal_id"].eq(signal_id)],
            comparison[comparison["signal_id"].eq(signal_id)].iloc[0],
        )
        classifications.append(
            {
                "signal_id": signal_id,
                "signal_label": overall["signal_label"],
                "evidence_role": overall["evidence_role"],
                "audit_classification": category,
                "rationale": rationale,
            }
        )
    classification = pd.DataFrame(classifications)
    synthesis = synthesis.merge(
        classification[["signal_id", "audit_classification"]], on="signal_id", how="left"
    )
    synthesis.to_csv(OUT / "cross_cohort_synthesis_v1.tsv", sep="\t", index=False)
    classification.to_csv(OUT / "signal_classification_v1.tsv", sep="\t", index=False)

    key_order = [
        "F1__TYMP_RNA",
        "F1__hypoxia_downstream_5gene",
        "F1__metaVIPER_HIF1A",
        "F1__scDRS_RCC_score",
        "F1__TYMP_axis_strict",
        "F1__TYMP_axis_broad",
        "F1__RRM2B_axis_strict",
        "F1__RRM_axis_broad",
        "F1__ODC1_RNA",
        "F1__ILF2_GENIE3",
        "F1__ILF2_GRNBoost2",
        "F1__HOXD3_GENIE3",
        "F1__HOXD3_GRNBoost2",
        "F1__HNF4A_GENIE3",
        "F1__HNF4A_GRNBoost2",
    ]
    key = synthesis.set_index("signal_id").loc[key_order].reset_index()
    lines = [
        "# Part 1 by-cohort signal audit v1",
        "",
        "## Scope and method",
        "",
        "This audit reuses frozen patient-level effects from the validated Part 1 analyses. No cell-level matrix, embedding, CNV call, cNMF factor, or score was recomputed. Cohorts are split across GSE156632, GSE159115, GSE178481, GSE207493, and GSE242299. The audit reports patient-equal medians, per-cohort direction, leave-one-dataset-out results, GSE207493-versus-other-cohort contrasts, and an exploratory DerSimonian-Laird random-effects synthesis of cohort mean Fisher-z effects when at least two patients are available per contributing cohort.",
        "",
        "## Key results",
        "",
    ]
    for _, row in key.iterrows():
        comp = comparison[comparison["signal_id"].eq(row["signal_id"])].iloc[0]
        no_dom = loo[
            loo["signal_id"].eq(row["signal_id"])
            & loo["left_out_study"].eq(DOMINANT_STUDY)
        ].iloc[0]
        lines.append(
            f"- **{row['signal_label']}**: all-patient median {fmt(row['median_effect'])} "
            f"({int(row['n_positive'])}/{int(row['n_patients'])} positive; Wilcoxon P={fmt(row['wilcoxon_p_two_sided'], 4)}, audit BH q={fmt(row['wilcoxon_fdr_bh_audit18'], 4)}); "
            f"study medians {int(row['n_positive_study_medians'])} positive / {int(row['n_negative_study_medians'])} negative; "
            f"GSE207493 median {fmt(comp['dominant_median_effect'])}, other cohorts median {fmt(comp['other_median_effect'])}; "
            f"after removing GSE207493 median {fmt(no_dom['median_effect'])}. "
            f"Classification: `{row['audit_classification']}`."
        )
    lines.extend(
        [
            "",
            "## Interpretation guardrails",
            "",
            "- A positive all-patient test is not by itself evidence of cross-cohort replication.",
            "- `CROSS_COHORT_DIRECTIONAL_SUPPORT` means the pooled patient-level test passes the 18-signal audit BH correction, direction survives removal of GSE207493 with at least 40% of the pooled median magnitude, and at least 75% of evaluable cohort medians agree; it does not require every cohort to be individually significant.",
            "- `GSE207493_WEIGHTED_DIRECTION_REPLICATED` means the other cohorts preserve direction, but most of the effect magnitude is contributed by GSE207493.",
            "- `COHORT_HETEROGENEITY_OPPOSITE` means the dominant cohort and other cohorts point in opposite directions; pooling can cancel effects and should not be described as simple dilution.",
            "- Cohorts with one eligible strict-malignant patient contribute direction only and cannot estimate within-cohort uncertainty for the random-effects calculation.",
            "- F1-state representation remains highly imbalanced (previous audit: 99.66% of F1-dominant cells from GSE207493). Therefore cross-cohort directional support for an association does not prove that the F1 state itself was independently rediscovered in every cohort.",
            "",
        ]
    )
    (PHASE / "AUDIT_SUMMARY_v1.md").write_text("\n".join(lines), encoding="utf-8")

    checks = [
        ("five_expected_studies_only", set(effects["study"]) <= set(EXPECTED_STUDIES)),
        ("all_five_studies_represented_somewhere", set(effects["study"]) == set(EXPECTED_STUDIES)),
        ("no_duplicate_signal_patient", not effects.duplicated(["signal_id", "curated_patient_id"]).any()),
        ("finite_effects", bool(np.isfinite(effects["effect"]).all())),
        ("core_tymp_present", "F1__TYMP_RNA" in set(effects["signal_id"])),
        ("all_key_signals_present", set(key_order) <= set(effects["signal_id"])),
        ("loo_gse207493_present_every_signal", int((loo["left_out_study"] == DOMINANT_STUDY).sum()) == effects["signal_id"].nunique()),
        ("classification_complete", classification["audit_classification"].notna().all()),
    ]
    check_df = pd.DataFrame(checks, columns=["check", "pass"])
    check_df.to_csv(PHASE / "validation_checks_v1.tsv", sep="\t", index=False)
    if not check_df["pass"].all():
        raise RuntimeError("Audit validation failed")

    manifest_rows = []
    for path in sorted(
        [
            *OUT.glob("*"),
            PHASE / "AUDIT_SUMMARY_v1.md",
            PHASE / "AUDIT_SUMMARY_ZH_v1.md",
            PHASE / "validation_checks_v1.tsv",
        ]
    ):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        manifest_rows.append(
            {"file": str(path.relative_to(PHASE)), "bytes": path.stat().st_size, "sha256": digest}
        )
    pd.DataFrame(manifest_rows).to_csv(
        PHASE / "FILE_MANIFEST_SHA256_v1.tsv", sep="\t", index=False
    )
    print(
        f"PASS: {effects['signal_id'].nunique()} signals, {len(effects)} patient-effects, "
        f"{len(per_cohort)} cohort summaries"
    )


if __name__ == "__main__":
    main()
