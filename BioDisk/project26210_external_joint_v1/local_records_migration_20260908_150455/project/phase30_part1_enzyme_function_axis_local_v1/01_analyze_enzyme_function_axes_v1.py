#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats


EXPECTED_HASHES = {
    "Human_M168_information.symbols.csv": "97d00e32b40dd9cb27f5af705c73c7e21f36fcdfb2d39c7d6eb0447637af4a92",
    "module_gene_m168.csv": "0e96d1c82e51038123d6a7c3972f91263397a6018ec6973356ba757b74ae6ca9",
    "m13_tymp_rrm2b_cell_expression_v1.tsv.gz": "3b5e9f34b3ca87b101acc7ddc436496f8e9323516491d1643a6c14fcff338806",
    "scFEA_four_axis_cell_scores_v1.tsv.gz": "a2162caacc31f1c82b2b94f1b235d56c69bb3dbb7922d7580682d0044f7355c4",
    "scFEA_M168_cell_flux.tsv.gz": "35d2ca0a431a1bfd67c9ea474601639cde2faa4ae2222acfd7b66bd98924a8f4",
    "trajectory_all_methods_cell_scores_v1.tsv.gz": "8162d22930858f7aa83975de232dd2111cedab38376fe7b3426aeafad695796a",
}

REACTIONS = {
    "M_45": ("Pyrimidine", "B-Alanine", "TYMP"),
    "M_164": ("dTMP", "Succinyl-CoA", "TYMP"),
    "M_157": ("CDP", "dCDP", "RRM2B"),
    "M_138": ("AMP", "Deoxyadenosine", "RRM2B"),
    "M_140": ("ADP", "Deoxyadenosine", "RRM2B"),
}

METHODS = {
    "Monocle3": "Monocle3_pseudotime_patient_v1",
    "Slingshot": "Slingshot_pseudotime_patient_v1",
    "scTour": "scTour_pseudotime_patient_v1",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def safe_patient_z(values: pd.Series) -> pd.Series:
    x = values.astype(float)
    sd = float(x.std(ddof=1))
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.zeros(len(x)), index=x.index, dtype=float)
    return (x - float(x.mean())) / sd


def signed_rank_p(values: pd.Series) -> float:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 3:
        return np.nan
    if np.allclose(x, 0):
        return 1.0
    return float(stats.wilcoxon(x, zero_method="wilcox", correction=False,
                                alternative="two-sided", method="auto").pvalue)


def bh(values: pd.Series) -> pd.Series:
    x = values.to_numpy(dtype=float)
    out = np.full(len(x), np.nan)
    ok = np.isfinite(x)
    if not ok.any():
        return pd.Series(out, index=values.index)
    idx = np.flatnonzero(ok)
    order = idx[np.argsort(x[idx])]
    ranked = x[order] * len(order) / np.arange(1, len(order) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out[order] = np.minimum(ranked, 1.0)
    return pd.Series(out, index=values.index)


def read_module_genes(path: Path) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.reader(handle):
            if row and row[0] in REACTIONS:
                result[row[0]] = [x for x in row[1:] if x]
    return result


def summarize_effects(detail: pd.DataFrame) -> pd.DataFrame:
    keys = ["feature", "feature_class", "endpoint", "method"]
    rows = []
    for key, d in detail.groupby(keys, sort=False, dropna=False):
        x = d["effect"].to_numpy(dtype=float)
        x = x[np.isfinite(x)]
        rows.append({
            **dict(zip(keys, key)),
            "n_eligible_patients": int(len(x)),
            "n_positive": int(np.sum(x > 0)),
            "n_negative": int(np.sum(x < 0)),
            "n_zero": int(np.sum(x == 0)),
            "median_patient_effect": float(np.median(x)) if len(x) else np.nan,
            "wilcoxon_p_two_sided": signed_rank_p(pd.Series(x)),
        })
    out = pd.DataFrame(rows)
    out["BH_FDR"] = out.groupby(["endpoint", "method", "feature_class"], dropna=False)[
        "wilcoxon_p_two_sided"
    ].transform(bh)
    n = out["n_eligible_patients"]
    positive = out["median_patient_effect"] > 0
    negative = out["median_patient_effect"] < 0
    inferential = n >= 10
    small = n.between(3, 9)
    out["positive_support"] = (
        (inferential & positive & (out["BH_FDR"] < 0.05))
        | (small & positive & (out["wilcoxon_p_two_sided"] < 0.05))
    )
    out["negative_support"] = (
        (inferential & negative & (out["BH_FDR"] < 0.05))
        | (small & negative & (out["wilcoxon_p_two_sided"] < 0.05))
    )
    out["reporting_tier"] = np.select(
        [inferential, small, n < 3],
        ["FDR_AT_LEAST_10", "NOMINAL_P_3_TO_9", "NOT_TESTED_LT3"],
        default="NOT_TESTED",
    )
    out["direction"] = np.select(
        [positive, negative], ["positive", "negative"], default="zero"
    )
    return out


def main() -> None:
    root = Path(__file__).resolve().parent
    input_dir = root / "input_snapshot_local_v1"
    output = root / "results_local_v1"
    part = root / "results_local_v1.part"
    if output.exists() or part.exists():
        raise FileExistsError("refusing to overwrite an existing Phase 30 result directory")
    part.mkdir(parents=True)

    manifest_rows = []
    for name, expected in EXPECTED_HASHES.items():
        path = input_dir / name
        observed = sha256(path)
        if observed != expected:
            raise ValueError(f"input hash mismatch: {name}: {observed}")
        manifest_rows.append({"file": name, "size_bytes": path.stat().st_size,
                              "sha256": observed})
    pd.DataFrame(manifest_rows).to_csv(part / "input_manifest_v1.tsv", sep="\t", index=False)

    info = pd.read_csv(input_dir / "Human_M168_information.symbols.csv")
    info = info.rename(columns={info.columns[0]: "row_id"})
    genes = read_module_genes(input_dir / "module_gene_m168.csv")
    reaction_rows = []
    for module, (source, target, required_gene) in REACTIONS.items():
        row = info.loc[info["row_id"].astype(str) == module]
        if len(row) != 1:
            raise ValueError(f"reaction row gate failed for {module}")
        observed_source = str(row.iloc[0]["Compound_IN_name"])
        observed_target = str(row.iloc[0]["Compound_OUT_name"])
        if (observed_source, observed_target) != (source, target):
            raise ValueError(f"reaction definition mismatch for {module}")
        if module not in genes or required_gene not in genes[module]:
            raise ValueError(f"required gene missing from {module}")
        reaction_rows.append({
            "module": module,
            "compound_in": source,
            "compound_out": target,
            "required_gene": required_gene,
            "n_module_genes": len(genes[module]),
            "module_genes": ";".join(genes[module]),
            "role": {
                "M_45": "TYMP_primary_broad_catabolic",
                "M_164": "TYMP_primary_broad_catabolic",
                "M_157": "RRM2B_primary_RNR_linked",
                "M_138": "RRM2B_broad_sensitivity",
                "M_140": "RRM2B_broad_sensitivity",
            }[module],
        })
    pd.DataFrame(reaction_rows).to_csv(part / "enzyme_axis_reaction_registry_v1.tsv",
                                       sep="\t", index=False)

    meta = pd.read_csv(input_dir / "scFEA_four_axis_cell_scores_v1.tsv.gz", sep="\t")
    flux = pd.read_csv(input_dir / "scFEA_M168_cell_flux.tsv.gz", sep="\t",
                       usecols=["cell", *REACTIONS.keys()]).rename(columns={"cell": "exact_cell_id"})
    traj_cols = ["exact_cell_id", "curated_patient_id", "study", *METHODS.values()]
    traj = pd.read_csv(input_dir / "trajectory_all_methods_cell_scores_v1.tsv.gz", sep="\t",
                       usecols=traj_cols)
    expr = pd.read_csv(input_dir / "m13_tymp_rrm2b_cell_expression_v1.tsv.gz", sep="\t")
    frames = {"metadata": meta, "flux": flux, "trajectory": traj, "expression": expr}
    for label, frame in frames.items():
        if len(frame) != 11950 or frame["exact_cell_id"].nunique() != 11950:
            raise ValueError(f"{label} 11,950-cell identity gate failed")
    cell_set = set(meta["exact_cell_id"])
    if any(set(frame["exact_cell_id"]) != cell_set for frame in frames.values()):
        raise ValueError("cross-input exact-cell set mismatch")

    data = meta.merge(flux, on="exact_cell_id", how="left", validate="one_to_one")
    data = data.merge(traj, on="exact_cell_id", how="left", validate="one_to_one",
                      suffixes=("", "_trajectory"))
    data = data.merge(expr, on="exact_cell_id", how="left", validate="one_to_one")
    if not (data["curated_patient_id"] == data["curated_patient_id_trajectory"]).all():
        raise ValueError("patient key mismatch")
    if not (data["study"] == data["study_trajectory"]).all():
        raise ValueError("study key mismatch")
    if not np.allclose(data["TYMP_log1p_CP10K_v1"], data["TYMP_log1p_cp10k"], atol=1e-12):
        raise ValueError("independent TYMP expression reconstruction mismatch")
    numeric_required = ["F1", "Pyrimidine_flux", "TYMP_log1p_cp10k",
                        "RRM2B_log1p_cp10k", *REACTIONS.keys(), *METHODS.values()]
    if not np.isfinite(data[numeric_required].to_numpy(dtype=float)).all():
        raise ValueError("non-finite joined numeric input")

    component_source = {
        "z_TYMP_RNA": "TYMP_log1p_cp10k",
        "z_TYMP_M45_flux": "M_45",
        "z_TYMP_M164_flux": "M_164",
        "z_TYMP_Pyrimidine_flux": "Pyrimidine_flux",
        "z_RRM2B_RNA": "RRM2B_log1p_cp10k",
        "z_RRM_M157_flux": "M_157",
        "z_RRM_M138_flux": "M_138",
        "z_RRM_M140_flux": "M_140",
    }
    constant_rows = []
    for target, source in component_source.items():
        stds = data.groupby("curated_patient_id", sort=False)[source].std(ddof=1)
        for patient, sd in stds.items():
            constant_rows.append({"component": target, "source_column": source,
                                  "curated_patient_id": patient,
                                  "source_sd": float(sd) if np.isfinite(sd) else np.nan,
                                  "constant_component": bool(not np.isfinite(sd) or sd == 0)})
        data[target] = data.groupby("curated_patient_id", sort=False)[source].transform(safe_patient_z)
    pd.DataFrame(constant_rows).to_csv(part / "enzyme_axis_component_qc_v1.tsv",
                                       sep="\t", index=False)

    data["TYMP_axis_primary"] = data[["z_TYMP_RNA", "z_TYMP_M45_flux",
                                      "z_TYMP_M164_flux"]].mean(axis=1)
    data["TYMP_axis_broad_sensitivity"] = data[["z_TYMP_RNA",
                                                "z_TYMP_Pyrimidine_flux"]].mean(axis=1)
    data["RRM2B_axis_primary"] = data[["z_RRM2B_RNA", "z_RRM_M157_flux"]].mean(axis=1)
    data["RRM2B_axis_broad_sensitivity"] = data[["z_RRM2B_RNA", "z_RRM_M138_flux",
                                                  "z_RRM_M140_flux",
                                                  "z_RRM_M157_flux"]].mean(axis=1)

    features = {
        "TYMP_axis_primary": "axis_primary",
        "RRM2B_axis_primary": "axis_primary",
        "TYMP_axis_broad_sensitivity": "axis_sensitivity",
        "RRM2B_axis_broad_sensitivity": "axis_sensitivity",
        **{name: "component" for name in component_source},
    }
    feature_registry = pd.DataFrame([
        {"feature": name, "feature_class": klass,
         "interpretation": {
             "TYMP_axis_primary": "TYMP-associated pyrimidine-catabolism function proxy",
             "RRM2B_axis_primary": "RRM2B-linked ribonucleotide-reduction function proxy",
             "TYMP_axis_broad_sensitivity": "TYMP RNA plus broad pyrimidine-flux sensitivity proxy",
             "RRM2B_axis_broad_sensitivity": "RRM2B RNA plus broad RNR-linked sensitivity proxy",
         }.get(name, "axis component; explanatory only")}
        for name, klass in features.items()
    ])
    feature_registry["is_VIPER_activity"] = False
    feature_registry["is_direct_biochemical_activity"] = False
    feature_registry.to_csv(part / "enzyme_axis_feature_registry_v1.tsv", sep="\t", index=False)

    score_columns = ["exact_cell_id", "study", "curated_patient_id", "F1", "clone_group",
                     *METHODS.values(), *component_source.keys(), *features.keys()]
    score_columns = list(dict.fromkeys(score_columns))
    data[score_columns].to_csv(part / "enzyme_axis_cell_scores_v1.tsv.gz", sep="\t",
                               index=False, compression="gzip")

    detail_rows = []
    patient_col = "curated_patient_id"
    for feature, feature_class in features.items():
        for patient, d in data.groupby(patient_col, sort=False):
            if len(d) >= 20 and d["F1"].nunique() > 1 and d[feature].nunique() > 1:
                effect = float(stats.spearmanr(d["F1"], d[feature]).statistic)
                detail_rows.append({"feature": feature, "feature_class": feature_class,
                                    "endpoint": "F1_continuous", "method": "Spearman",
                                    patient_col: patient, "n_cells": len(d), "n_low": np.nan,
                                    "n_high": np.nan, "effect": effect})

            low = d[d["clone_group"] == "other_primary_clones"]
            high = d[d["clone_group"] == "F1_high_primary_clone"]
            if len(low) >= 10 and len(high) >= 10:
                detail_rows.append({"feature": feature, "feature_class": feature_class,
                                    "endpoint": "F1_dominant_group", "method": "median_delta",
                                    patient_col: patient, "n_cells": len(low) + len(high),
                                    "n_low": len(low), "n_high": len(high),
                                    "effect": float(high[feature].median() - low[feature].median())})

            if len(d) >= 20 and d["F1"].nunique() > 1:
                q25, q75 = d["F1"].quantile([0.25, 0.75]).to_numpy()
                low_q = d[d["F1"] <= q25]
                high_q = d[d["F1"] >= q75]
                if q75 > q25 and len(low_q) >= 5 and len(high_q) >= 5:
                    detail_rows.append({"feature": feature, "feature_class": feature_class,
                                        "endpoint": "F1_patient_quartile", "method": "median_delta",
                                        patient_col: patient, "n_cells": len(d),
                                        "n_low": len(low_q), "n_high": len(high_q),
                                        "effect": float(high_q[feature].median() - low_q[feature].median())})

            for method, trajectory_col in METHODS.items():
                if len(d) < 20 or d[trajectory_col].nunique() <= 1:
                    continue
                q20, q80 = d[trajectory_col].quantile([0.20, 0.80]).to_numpy()
                early = d[d[trajectory_col] <= q20]
                late = d[d[trajectory_col] >= q80]
                if q80 > q20 and len(early) >= 5 and len(late) >= 5:
                    detail_rows.append({"feature": feature, "feature_class": feature_class,
                                        "endpoint": "trajectory_terminal_shift", "method": method,
                                        patient_col: patient, "n_cells": len(d),
                                        "n_low": len(early), "n_high": len(late),
                                        "effect": float(late[feature].median() - early[feature].median())})

    detail = pd.DataFrame(detail_rows)
    detail.to_csv(part / "enzyme_axis_patient_effects_v1.tsv", sep="\t", index=False,
                  na_rep="NA")
    summary = summarize_effects(detail)
    summary.to_csv(part / "enzyme_axis_endpoint_summary_v1.tsv", sep="\t", index=False,
                   na_rep="NA")

    evidence_rows = []
    for axis in ["TYMP_axis_primary", "RRM2B_axis_primary"]:
        s = summary[summary["feature"] == axis]
        support = {
            "F1_continuous": bool(s.loc[s["endpoint"] == "F1_continuous", "positive_support"].iloc[0]),
            "F1_dominant_group": bool(s.loc[s["endpoint"] == "F1_dominant_group", "positive_support"].iloc[0]),
            "F1_patient_quartile": bool(s.loc[s["endpoint"] == "F1_patient_quartile", "positive_support"].iloc[0]),
        }
        traj_support = int(s.loc[s["endpoint"] == "trajectory_terminal_shift",
                                 "positive_support"].sum())
        support["trajectory_convergence_2_of_3"] = traj_support >= 2
        count = int(sum(support.values()))
        tier = "HIGH_SUPPORT" if count == 4 else "SECONDARY_SUPPORT" if count >= 3 else "LIMITED_OR_UNSUPPORTED"
        evidence_rows.append({"axis": axis, **support,
                              "n_positive_trajectory_methods": traj_support,
                              "n_supported_endpoints_of_4": count, "evidence_tier": tier})
    evidence = pd.DataFrame(evidence_rows)
    evidence.to_csv(part / "enzyme_axis_integrated_evidence_v1.tsv", sep="\t", index=False)

    terminal = detail[(detail["endpoint"] == "trajectory_terminal_shift") &
                      (detail["feature"].isin(["TYMP_axis_primary", "RRM2B_axis_primary"]))]
    joint_rows = []
    for method in METHODS:
        a = terminal[(terminal["method"] == method) &
                     (terminal["feature"] == "TYMP_axis_primary")][[patient_col, "effect"]].rename(columns={"effect": "TYMP_effect"})
        b = terminal[(terminal["method"] == method) &
                     (terminal["feature"] == "RRM2B_axis_primary")][[patient_col, "effect"]].rename(columns={"effect": "RRM2B_effect"})
        d = a.merge(b, on=patient_col, validate="one_to_one")
        rho = float(stats.spearmanr(d["TYMP_effect"], d["RRM2B_effect"]).statistic) if len(d) >= 3 else np.nan
        joint_rows.append({"method": method, "n_patients": len(d),
                           "n_both_positive": int(((d["TYMP_effect"] > 0) & (d["RRM2B_effect"] > 0)).sum()),
                           "fraction_both_positive": float(((d["TYMP_effect"] > 0) & (d["RRM2B_effect"] > 0)).mean()),
                           "median_joint_min_shift": float(np.median(np.minimum(d["TYMP_effect"], d["RRM2B_effect"]))),
                           "spearman_patient_shift": rho})
    pd.DataFrame(joint_rows).to_csv(part / "enzyme_axis_joint_terminal_summary_v1.tsv",
                                    sep="\t", index=False)

    status = {
        "status": "PASS_ANALYSIS_COMPLETE_AWAITING_INDEPENDENT_VALIDATION",
        "n_cells": int(len(data)),
        "n_patients_total": int(data[patient_col].nunique()),
        "n_primary_axes": 2,
        "n_sensitivity_axes": 2,
        "n_primary_pseudotime_methods": 3,
        "viper_activity_inferred": False,
        "direct_biochemical_activity_measured": False,
        "figures_created": False,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
    }
    (part / "status.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    os.replace(part, output)


if __name__ == "__main__":
    main()
