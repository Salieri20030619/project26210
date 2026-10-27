#!/usr/bin/env python3
import argparse
import hashlib
import json
import math
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats


SCORES = {
    "F1": "score_primary__f1_top100_weighted",
    "TYMP": "score_primary__tymp_expression",
    "SPP1_state_no_SPP1": "score_primary__spp1_state_no_spp1",
    "C1QC_state": "score_primary__c1qc_state",
    "APOE_TAM": "score_primary__apoe_tam_expression_proxy",
    "monocyte_program": "score_primary__monocyte_program",
    "myeloid_suppression": "score_primary__myeloid_immunosuppression",
    "vascular_remodeling": "score_primary__vascular_remodeling",
    "T_cell_dysfunction": "score_primary__tcell_dysfunction",
    "cytotoxicity": "score_primary__cytotoxicity",
    "POSTN_ECM": "score_primary__postn_ecm_proxy",
    "HLA_II_support": "score_primary__hla_ii_supporting",
}

ABUND = {
    "Epithelial": ["Epithelial cell"],
    "Fibroblast": ["Fibroblast"],
    "Macrophage": ["Macrophage"],
    "Monocyte": ["Monocyte"],
    "Myeloid": ["Macrophage", "Monocyte"],
    "Vascular": ["Endothelial", "Pericyte/SMC"],
    "T_cell": ["T cell"],
    "Cytotoxic_lymphocyte": ["T cell", "NK cell"],
    "HLA_II_compartment": ["B cell", "Dendritic cell", "Macrophage"],
}

CONSTRUCTS = [
    ("C01_F1_to_Epithelial", "F1", "Epithelial", 1),
    ("C02_POSTN_to_Fibroblast", "POSTN_ECM", "Fibroblast", 1),
    ("C03_APOE_to_Macrophage", "APOE_TAM", "Macrophage", 1),
    ("C04_C1QC_to_Macrophage", "C1QC_state", "Macrophage", 1),
    ("C05_SPP1_to_Macrophage", "SPP1_state_no_SPP1", "Macrophage", 1),
    ("C06_Monocyte_program_to_Monocyte", "monocyte_program", "Monocyte", 1),
    ("C07_Myeloid_suppression_to_Myeloid", "myeloid_suppression", "Myeloid", 1),
    ("C08_Vascular_to_Vascular", "vascular_remodeling", "Vascular", 1),
    ("C09_Tdysfunction_to_Tcell", "T_cell_dysfunction", "T_cell", 1),
    ("C10_Cytotoxicity_to_T_or_NK", "cytotoxicity", "Cytotoxic_lymphocyte", 1),
    ("C11_HLAII_to_APC", "HLA_II_support", "HLA_II_compartment", 1),
    ("C12_TYMP_to_Myeloid", "TYMP", "Myeloid", 1),
]

BRIDGES = [
    ("B01_TYMP_to_SPP1", "TYMP", "SPP1_state_no_SPP1", ["Macrophage", "Monocyte"], 1),
    ("B02_TYMP_to_C1QC", "TYMP", "C1QC_state", ["Macrophage", "Monocyte"], 1),
    ("B03_TYMP_to_APOE", "TYMP", "APOE_TAM", ["Macrophage", "Monocyte"], 1),
    ("B04_F1_to_SPP1", "F1", "SPP1_state_no_SPP1", ["Epithelial", "Macrophage", "Monocyte"], 1),
    ("B05_F1_to_C1QC", "F1", "C1QC_state", ["Epithelial", "Macrophage", "Monocyte"], 1),
    ("B06_F1_to_APOE", "F1", "APOE_TAM", ["Epithelial", "Macrophage", "Monocyte"], 1),
    ("B07_POSTN_to_SPP1", "POSTN_ECM", "SPP1_state_no_SPP1", ["Fibroblast", "Macrophage"], 1),
    ("B08_SPP1_to_Tdysfunction", "SPP1_state_no_SPP1", "T_cell_dysfunction", ["Macrophage", "T_cell"], 1),
    ("B09_SPP1_to_Vascular", "SPP1_state_no_SPP1", "vascular_remodeling", ["Macrophage", "Vascular"], 1),
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bh_adjust(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def prep_for(sample_id):
    n = int(sample_id.replace("GSM", ""))
    return "FFPE" if n <= 5924041 else "FROZEN"


def abundance_column(prefix, label):
    return f"{prefix}__{label}"


def build_abundance(df, prefix):
    out = {}
    for name, members in ABUND.items():
        cols = [abundance_column(prefix, member) for member in members]
        missing = [col for col in cols if col not in df.columns]
        if missing:
            raise RuntimeError(f"Missing cell2location abundance columns: {missing}")
        out[name] = df[cols].sum(axis=1).to_numpy(dtype=float)
    return out


def partial_spearman(x, y, covariates):
    matrix = np.column_stack([x, y] + covariates)
    keep = np.all(np.isfinite(matrix), axis=1)
    matrix = matrix[keep]
    n = matrix.shape[0]
    k = len(covariates)
    if n <= k + 4:
        return math.nan, math.nan, n
    ranked = np.column_stack([stats.rankdata(matrix[:, i], method="average") for i in range(matrix.shape[1])])
    design = np.column_stack([np.ones(n), ranked[:, 2:]])
    rx = ranked[:, 0] - design @ np.linalg.lstsq(design, ranked[:, 0], rcond=None)[0]
    ry = ranked[:, 1] - design @ np.linalg.lstsq(design, ranked[:, 1], rcond=None)[0]
    if np.std(rx) == 0 or np.std(ry) == 0:
        return math.nan, math.nan, n
    rho = float(np.corrcoef(rx, ry)[0, 1])
    dfree = n - k - 2
    if dfree <= 0 or abs(rho) >= 1:
        p = 0.0 if abs(rho) >= 1 else math.nan
    else:
        tval = rho * math.sqrt(dfree / max(1e-15, 1.0 - rho * rho))
        p = float(2 * stats.t.sf(abs(tval), dfree))
    return rho, p, n


def meta_analyze(per_section, group_cols):
    rows = []
    for keys, sub in per_section.groupby(group_cols, sort=True, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        values = sub["rho"].dropna().to_numpy(dtype=float)
        expected = int(sub["expected_direction"].iloc[0])
        if len(values) == 0:
            p = math.nan
        elif np.allclose(values, 0):
            p = 1.0
        else:
            p = float(stats.wilcoxon(values, alternative="two-sided", zero_method="wilcox").pvalue)
        row = dict(zip(group_cols, keys))
        row.update(
            n_sections=len(values),
            median_rho=float(np.median(values)) if len(values) else math.nan,
            min_rho=float(np.min(values)) if len(values) else math.nan,
            max_rho=float(np.max(values)) if len(values) else math.nan,
            expected_direction=expected,
            expected_direction_fraction=float(np.mean(values * expected > 0)) if len(values) else math.nan,
            wilcoxon_signed_rank_p=p,
        )
        rows.append(row)
    meta = pd.DataFrame(rows)
    meta["wilcoxon_bh_q"] = np.nan
    for _, idx in meta.groupby(["preparation", "analysis_family"], sort=True).groups.items():
        idx = list(idx)
        valid = [i for i in idx if np.isfinite(meta.loc[i, "wilcoxon_signed_rank_p"])]
        if valid:
            meta.loc[valid, "wilcoxon_bh_q"] = bh_adjust(meta.loc[valid, "wilcoxon_signed_rank_p"].to_numpy())
    meta["support_status"] = np.where(
        (meta["wilcoxon_bh_q"] < 0.05)
        & (meta["median_rho"] * meta["expected_direction"] > 0)
        & (meta["expected_direction_fraction"] >= (8 / 12)),
        "SUPPORTED",
        "NOT_SUPPORTED",
    )
    return meta


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase54", required=True)
    parser.add_argument("--phase55", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    phase54 = Path(args.phase54).resolve()
    phase55 = Path(args.phase55).resolve()
    out = Path(args.out).resolve()
    if out.exists():
        raise RuntimeError(f"Refusing to overwrite existing output: {out}")
    out_part = Path(str(out) + ".part")
    if out_part.exists():
        raise RuntimeError(f"Refusing to overwrite partial output: {out_part}")
    out_part.mkdir(parents=True)

    pass54 = phase54 / "postrun_validation_v1" / "POSTRUN_VALIDATION_PASS_v1"
    pass55 = phase55 / "postrun_validation_v2" / "POSTRUN_VALIDATION_PASS_v2"
    if not pass54.is_file() or not pass55.is_file():
        raise RuntimeError("Frozen Phase54/Phase55 PASS gates are missing")

    sec54 = phase54 / "run_v1" / "sections_v1"
    sec55 = phase55 / "run_v1" / "sections_v1"
    ids54 = {p.name for p in sec54.iterdir() if p.is_dir()}
    ids55 = {p.name for p in sec55.iterdir() if p.is_dir()}
    if ids54 != ids55 or len(ids54) != 24:
        raise RuntimeError(f"Section mismatch: phase54={len(ids54)}, phase55={len(ids55)}")

    section_qc = []
    construct_rows = []
    bridge_rows = []
    input_manifest = []

    for sample_id in sorted(ids54):
        score_path = sec54 / sample_id / "spot_metadata_scores_v1.tsv.gz"
        abundance_path = sec55 / sample_id / "spot_cell_abundance_v1.tsv.gz"
        input_manifest.extend(
            [
                {"sample_id": sample_id, "input_type": "phase54_scores", "path": str(score_path), "bytes": score_path.stat().st_size, "sha256": sha256(score_path)},
                {"sample_id": sample_id, "input_type": "phase55_abundance", "path": str(abundance_path), "bytes": abundance_path.stat().st_size, "sha256": sha256(abundance_path)},
            ]
        )
        scores = pd.read_csv(score_path, sep="\t")
        abundance = pd.read_csv(abundance_path, sep="\t")
        scores = scores.loc[scores["qc_pass"].astype(str).str.lower().isin(["true", "1"])].copy()
        required_score_cols = ["barcode", "total_counts"] + list(SCORES.values())
        missing_scores = [c for c in required_score_cols if c not in scores.columns]
        if missing_scores:
            raise RuntimeError(f"{sample_id}: missing Phase54 columns {missing_scores}")
        if scores["barcode"].duplicated().any() or abundance["barcode"].duplicated().any():
            raise RuntimeError(f"{sample_id}: duplicated barcodes")
        merged = scores.merge(abundance, on="barcode", how="inner", validate="one_to_one", suffixes=("", "__c2l"))
        exact = set(scores["barcode"]) == set(abundance["barcode"]) == set(merged["barcode"])
        if not exact:
            raise RuntimeError(f"{sample_id}: Phase54/55 barcode sets do not match exactly")
        prep = prep_for(sample_id)
        log_umi = np.log1p(merged["total_counts"].to_numpy(dtype=float))
        mean_abund = build_abundance(merged, "mean")
        q05_abund = build_abundance(merged, "q05")
        score_values = {name: merged[col].to_numpy(dtype=float) for name, col in SCORES.items()}
        section_qc.append(
            {
                "sample_id": sample_id,
                "preparation": prep,
                "spots": len(merged),
                "exact_barcode_match": "PASS",
                "finite_score_fraction": float(np.isfinite(merged[list(SCORES.values())].to_numpy(dtype=float)).mean()),
                "finite_mean_abundance_fraction": float(np.isfinite(merged.filter(regex=r"^mean__").to_numpy(dtype=float)).mean()),
            }
        )

        for relation_id, score_name, abundance_name, direction in CONSTRUCTS:
            for abundance_metric, abundance_values in [("mean", mean_abund), ("q05", q05_abund)]:
                rho, p = stats.spearmanr(score_values[score_name], abundance_values[abundance_name], nan_policy="omit")
                construct_rows.append(
                    {
                        "sample_id": sample_id,
                        "preparation": prep,
                        "analysis_family": "supporting_construct_mean" if abundance_metric == "mean" else "sensitivity_construct_q05",
                        "relation_id": relation_id,
                        "score": score_name,
                        "abundance": abundance_name,
                        "abundance_metric": abundance_metric,
                        "expected_direction": direction,
                        "n_spots": len(merged),
                        "rho": float(rho),
                        "spot_level_p_descriptive": float(p),
                    }
                )

        for relation_id, left, right, cov_names, direction in BRIDGES:
            covs = [mean_abund[name] for name in cov_names] + [log_umi]
            rho, p, n = partial_spearman(score_values[left], score_values[right], covs)
            simple_rho, simple_p = stats.spearmanr(score_values[left], score_values[right], nan_policy="omit")
            bridge_rows.append(
                {
                    "sample_id": sample_id,
                    "preparation": prep,
                    "analysis_family": "primary_composition_adjusted_bridge",
                    "relation_id": relation_id,
                    "left_score": left,
                    "right_score": right,
                    "covariates": ";".join(cov_names + ["log1p_total_counts"]),
                    "expected_direction": direction,
                    "n_spots": n,
                    "rho": rho,
                    "spot_level_p_descriptive": p,
                    "unadjusted_spearman_rho": float(simple_rho),
                    "unadjusted_spot_level_p_descriptive": float(simple_p),
                }
            )

    qc = pd.DataFrame(section_qc)
    constructs = pd.DataFrame(construct_rows)
    bridges = pd.DataFrame(bridge_rows)
    construct_meta = meta_analyze(constructs, ["preparation", "analysis_family", "relation_id", "score", "abundance", "abundance_metric"])
    bridge_meta = meta_analyze(bridges, ["preparation", "analysis_family", "relation_id", "left_score", "right_score", "covariates"])

    decision_rows = []
    for relation_id, sub in bridge_meta.groupby("relation_id", sort=True):
        statuses = dict(zip(sub["preparation"], sub["support_status"]))
        n_supported = sum(v == "SUPPORTED" for v in statuses.values())
        decision_rows.append(
            {
                "relation_id": relation_id,
                "ffpe_status": statuses.get("FFPE", "MISSING"),
                "frozen_status": statuses.get("FROZEN", "MISSING"),
                "preparations_supported": n_supported,
                "decision": "TWO_PREPARATION_SUPPORTED" if n_supported == 2 else ("ONE_PREPARATION_SUPPORTED" if n_supported == 1 else "NOT_SUPPORTED"),
            }
        )
    decisions = pd.DataFrame(decision_rows)

    gates = [
        ("phase54_pass_gate", int(pass54.is_file()), 1),
        ("phase55_pass_gate", int(pass55.is_file()), 1),
        ("sections_total", len(qc), 24),
        ("ffpe_sections", int((qc["preparation"] == "FFPE").sum()), 12),
        ("frozen_sections", int((qc["preparation"] == "FROZEN").sum()), 12),
        ("exact_barcode_matches", int((qc["exact_barcode_match"] == "PASS").sum()), 24),
        ("construct_rows", len(constructs), 24 * len(CONSTRUCTS) * 2),
        ("bridge_rows", len(bridges), 24 * len(BRIDGES)),
        ("construct_meta_rows", len(construct_meta), 2 * len(CONSTRUCTS) * 2),
        ("bridge_meta_rows", len(bridge_meta), 2 * len(BRIDGES)),
        ("decision_rows", len(decisions), len(BRIDGES)),
        ("finite_bridge_rho", int(np.isfinite(bridges["rho"]).sum()), len(bridges)),
    ]
    status = pd.DataFrame(gates, columns=["gate", "observed", "expected"])
    status["status"] = np.where(status["observed"] == status["expected"], "PASS", "FAIL")
    if (status["status"] != "PASS").any():
        raise RuntimeError("Phase70 validation gates failed before finalization")

    qc.to_csv(out_part / "section_join_qc_v1.tsv", sep="\t", index=False)
    constructs.to_csv(out_part / "construct_correlations_by_section_v1.tsv", sep="\t", index=False)
    bridges.to_csv(out_part / "composition_adjusted_bridges_by_section_v1.tsv", sep="\t", index=False)
    construct_meta.to_csv(out_part / "construct_cross_section_meta_v1.tsv", sep="\t", index=False)
    bridge_meta.to_csv(out_part / "composition_adjusted_bridge_meta_v1.tsv", sep="\t", index=False)
    decisions.to_csv(out_part / "primary_decisions_v1.tsv", sep="\t", index=False)
    status.to_csv(out_part / "PHASE70_RUN_STATUS_v1.tsv", sep="\t", index=False)
    pd.DataFrame(input_manifest).to_csv(out_part / "INPUT_MANIFEST_SHA256_v1.tsv", sep="\t", index=False)

    supported2 = decisions.loc[decisions["decision"] == "TWO_PREPARATION_SUPPORTED", "relation_id"].tolist()
    supported1 = decisions.loc[decisions["decision"] == "ONE_PREPARATION_SUPPORTED", "relation_id"].tolist()
    summary = [
        "# Phase70 composition-aware spatial validation",
        "",
        f"- Phase54 and Phase55 spot barcodes matched exactly in {len(qc)}/24 sections ({int(qc['spots'].sum()):,} spots).",
        f"- Two-preparation supported primary bridges: {len(supported2)}/{len(BRIDGES)}.",
        f"- One-preparation-only primary bridges: {len(supported1)}/{len(BRIDGES)}.",
        "",
        "## Two-preparation supported",
        "",
    ]
    summary += [f"- {x}" for x in supported2] if supported2 else ["- None."]
    summary += ["", "## One-preparation-only", ""]
    summary += [f"- {x}" for x in supported1] if supported1 else ["- None."]
    summary += [
        "",
        "## Interpretation boundary",
        "",
        "Section-level coefficients, not spots, are the inferential units. Composition adjustment uses coarse cell2location classes and cannot prove cell of origin or causality. This stage does not alter the frozen response/PFS evidence grade.",
    ]
    (out_part / "PHASE70_COMPOSITION_AWARE_SPATIAL_SUMMARY_v1.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    metadata = {
        "status": "PASS",
        "phase54": str(phase54),
        "phase55": str(phase55),
        "sections": len(qc),
        "spots": int(qc["spots"].sum()),
        "two_preparation_supported": len(supported2),
        "one_preparation_supported": len(supported1),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
    }
    (out_part / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    output_manifest = []
    for path in sorted(out_part.iterdir()):
        if path.is_file() and path.name != "OUTPUT_MANIFEST_SHA256_v1.tsv":
            output_manifest.append({"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(output_manifest).to_csv(out_part / "OUTPUT_MANIFEST_SHA256_v1.tsv", sep="\t", index=False)
    (out_part / "PHASE70_PASS_v1").write_text("status=PASS\n", encoding="utf-8")
    out_part.rename(out)


if __name__ == "__main__":
    main()
