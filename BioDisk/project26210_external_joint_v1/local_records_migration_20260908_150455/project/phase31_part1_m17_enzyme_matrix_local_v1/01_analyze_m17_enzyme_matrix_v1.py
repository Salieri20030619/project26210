#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
from itertools import combinations
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy
from scipy import sparse, stats


PANELS = {
    "thymidine_transport_salvage": ["TYMP", "TK1", "SLC29A1", "SLC29A2"],
    "polyamine_MTA": ["ODC1", "SRM", "SMS", "MTAP"],
    "kennedy_phospholipid": ["CHKA", "PCYT1A", "PCYT2", "CEPT1"],
}
GENES = [gene for panel in PANELS.values() for gene in panel]
EXPECTED_HASHES = {
    "malignant_eligible_raw_counts_phase25_v1.h5ad": "e9309d1cbccb2fb0621f38d89d36407866d51e94ba26e9b03c812dd6a3f5e35b",
    "scdrs_factor_tymp_cell_source_v1.tsv.gz": "3eeef447714e7a172d22317c9d9f64ac6fb9afe082fe0d125ddb55e2c9a5a247",
    "INPUT_FILE_MANIFEST_v1.tsv": "a0276e5332cf19ae59f9f98acc153ce8b3e81ad16070f883e3eb8713443c0203",
    "phase25_input_status.json": "d4a35291d3cac7e32d78887d10f80816d8c3c635d8c072f97a0e4104abd6c0a7",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def decode(values: np.ndarray) -> np.ndarray:
    return np.asarray([x.decode("utf-8") if isinstance(x, (bytes, np.bytes_)) else str(x)
                       for x in values])


def read_obs(handle: h5py.File, name: str) -> np.ndarray:
    obj = handle["obs"][name]
    if isinstance(obj, h5py.Group):
        categories = decode(obj["categories"][:])
        codes = obj["codes"][:]
        if (codes < 0).any():
            raise ValueError(f"missing categorical value in obs/{name}")
        return categories[codes]
    values = obj[:]
    return decode(values) if values.dtype.kind in "OSU" else values


def signed_rank_p(values: np.ndarray) -> float:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 3:
        return np.nan
    if np.allclose(x, 0):
        return 1.0
    return float(stats.wilcoxon(x, zero_method="wilcox", correction=False,
                                alternative="two-sided", method="auto").pvalue)


def bh(values: pd.Series) -> pd.Series:
    x = values.to_numpy(dtype=float)
    out = np.full(len(x), np.nan)
    good = np.isfinite(x)
    idx = np.flatnonzero(good)
    if not len(idx):
        return pd.Series(out, index=values.index)
    order = idx[np.argsort(x[idx])]
    adjusted = x[order] * len(order) / np.arange(1, len(order) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out[order] = np.minimum(adjusted, 1.0)
    return pd.Series(out, index=values.index)


def classify_summary(frame: pd.DataFrame, family: list[str]) -> pd.DataFrame:
    out = frame.copy()
    out["BH_FDR"] = out.groupby(family, dropna=False)["wilcoxon_p_two_sided"].transform(bh)
    n = out["n_eligible_patients"]
    out["statistical_support"] = ((n >= 10) & (out["BH_FDR"] < 0.05)) | (
        n.between(3, 9) & (out["wilcoxon_p_two_sided"] < 0.05)
    )
    out["reporting_tier"] = np.select(
        [n >= 10, n.between(3, 9), n < 3],
        ["FDR_AT_LEAST_10", "NOMINAL_P_3_TO_9", "NOT_TESTED_LT3"],
        default="NOT_TESTED",
    )
    out["direction"] = np.select(
        [out["median_patient_effect"] > 0, out["median_patient_effect"] < 0],
        ["positive", "negative"], default="zero"
    )
    out["supported_direction"] = np.where(out["statistical_support"], out["direction"], "unsupported")
    return out


def summarize(detail: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    rows = []
    for key, d in detail.groupby(keys, sort=False, dropna=False):
        key = key if isinstance(key, tuple) else (key,)
        x = d["effect"].to_numpy(dtype=float)
        x = x[np.isfinite(x)]
        rows.append({**dict(zip(keys, key)),
                     "n_eligible_patients": len(x),
                     "n_positive": int(np.sum(x > 0)),
                     "n_negative": int(np.sum(x < 0)),
                     "n_zero": int(np.sum(x == 0)),
                     "median_patient_effect": float(np.median(x)) if len(x) else np.nan,
                     "wilcoxon_p_two_sided": signed_rank_p(x)})
    return pd.DataFrame(rows)


def main() -> None:
    root = Path(__file__).resolve().parent
    inputs = root / "input_snapshot_local_v1"
    output = root / "results_local_v1"
    part = root / "results_local_v1.part"
    if output.exists() or part.exists():
        raise FileExistsError("refusing to overwrite Phase 31 results")
    part.mkdir(parents=True)

    manifest_rows = []
    for name, expected in EXPECTED_HASHES.items():
        path = inputs / name
        observed = sha256(path)
        if observed != expected:
            raise ValueError(f"input hash mismatch: {name}")
        manifest_rows.append({"file": name, "size_bytes": path.stat().st_size, "sha256": observed})
    pd.DataFrame(manifest_rows).to_csv(part / "input_manifest_v1.tsv", sep="\t", index=False)

    h5_path = inputs / "malignant_eligible_raw_counts_phase25_v1.h5ad"
    with h5py.File(h5_path, "r") as handle:
        shape = tuple(int(x) for x in handle["X"].attrs["shape"])
        if shape != (11950, 37795) or handle["X"].attrs.get("encoding-type") != "csr_matrix":
            raise ValueError("H5AD matrix contract failed")
        matrix_data = handle["X/data"][:]
        matrix_indices = handle["X/indices"][:]
        matrix_indptr = handle["X/indptr"][:]
        if (not np.isfinite(matrix_data).all() or (matrix_data < 0).any()
                or not np.equal(matrix_data, np.floor(matrix_data)).all()):
            raise ValueError("raw-count safety gate failed")
        matrix = sparse.csr_matrix((matrix_data, matrix_indices, matrix_indptr), shape=shape)
        genes = decode(handle["var/gene"][:])
        exact_cell_id = read_obs(handle, "exact_cell_id")
        index_id = decode(handle["obs/_index"][:])
        study = read_obs(handle, "study")
        patient = read_obs(handle, "curated_patient_id")
        f1 = read_obs(handle, "F1").astype(float)
        dominant = read_obs(handle, "dominant_factor_rank7_v1")
        library_recorded = read_obs(handle, "raw_library_size_v1").astype(float)
    del matrix_data, matrix_indices, matrix_indptr
    if len(set(GENES)) != 12 or not set(GENES).issubset(set(genes)):
        raise ValueError("target gene registry gate failed")
    if len(np.unique(exact_cell_id)) != 11950 or not np.array_equal(exact_cell_id, index_id):
        raise ValueError("cell identity gate failed")
    if len(np.unique(patient)) != 26 or len(np.unique(study)) != 5:
        raise ValueError("cohort gate failed")
    if not np.isfinite(f1).all():
        raise ValueError("non-finite F1")
    library = np.asarray(matrix.sum(axis=1)).ravel().astype(float)
    if (library <= 0).any() or not np.allclose(library, library_recorded, atol=0, rtol=0):
        raise ValueError("library-size conservation failed")

    cross = pd.read_csv(inputs / "scdrs_factor_tymp_cell_source_v1.tsv.gz", sep="\t",
                        usecols=["exact_cell_id", "study", "curated_patient_id", "F1"])
    cross = cross.set_index("exact_cell_id").loc[exact_cell_id].reset_index()
    if not (np.array_equal(cross["study"].astype(str), study)
            and np.array_equal(cross["curated_patient_id"].astype(str), patient)
            and np.allclose(cross["F1"], f1, atol=1e-15)):
        raise ValueError("independent cell-metadata cross-check failed")

    target_indices = [int(np.flatnonzero(genes == gene)[0]) for gene in GENES]
    raw = matrix[:, target_indices].toarray().astype(float)
    log_cp10k = np.log1p(raw / library[:, None] * 10000.0)
    data = pd.DataFrame({"exact_cell_id": exact_cell_id, "study": study,
                         "curated_patient_id": patient, "F1": f1,
                         "dominant_factor_rank7_v1": dominant,
                         "raw_library_size": library})
    qc_rows = []
    for i, gene in enumerate(GENES):
        data[f"{gene}_raw_count"] = raw[:, i]
        data[f"{gene}_log1p_cp10k"] = log_cp10k[:, i]
        panel = next(k for k, v in PANELS.items() if gene in v)
        qc_rows.append({"gene": gene, "panel": panel, "h5ad_gene_index": target_indices[i],
                        "detected_cells": int(np.sum(raw[:, i] > 0)),
                        "detection_fraction": float(np.mean(raw[:, i] > 0)),
                        "total_counts": float(raw[:, i].sum()),
                        "n_patients_with_detection": int(pd.DataFrame({"patient": patient, "x": raw[:, i]})
                                                         .groupby("patient")["x"].apply(lambda x: (x > 0).any()).sum())})
    pd.DataFrame(qc_rows).to_csv(part / "m17_gene_detection_qc_v1.tsv", sep="\t", index=False)
    registry = pd.DataFrame([{"gene": gene, "panel": panel, "prespecified": True}
                             for panel, genes_in_panel in PANELS.items() for gene in genes_in_panel])
    registry.to_csv(part / "m17_gene_registry_v1.tsv", sep="\t", index=False)
    data.to_csv(part / "m17_cell_expression_v1.tsv.gz", sep="\t", index=False, compression="gzip")

    effect_rows = []
    for gene in GENES:
        value_col = f"{gene}_log1p_cp10k"
        panel = next(k for k, v in PANELS.items() if gene in v)
        for pid, d in data.groupby("curated_patient_id", sort=False):
            if len(d) >= 20 and d["F1"].nunique() > 1 and d[value_col].nunique() > 1:
                effect_rows.append({"gene": gene, "panel": panel, "endpoint": "F1_continuous",
                                    "method": "Spearman", "curated_patient_id": pid,
                                    "n_cells": len(d), "n_low": np.nan, "n_high": np.nan,
                                    "effect": float(stats.spearmanr(d["F1"], d[value_col]).statistic)})
            low = d[d["dominant_factor_rank7_v1"] != "F1"]
            high = d[d["dominant_factor_rank7_v1"] == "F1"]
            if len(low) >= 10 and len(high) >= 10:
                effect_rows.append({"gene": gene, "panel": panel, "endpoint": "F1_dominant_group",
                                    "method": "median_delta", "curated_patient_id": pid,
                                    "n_cells": len(low) + len(high), "n_low": len(low),
                                    "n_high": len(high),
                                    "effect": float(high[value_col].median() - low[value_col].median())})
            if len(d) >= 20 and d["F1"].nunique() > 1:
                q25, q75 = d["F1"].quantile([0.25, 0.75]).to_numpy()
                low_q = d[d["F1"] <= q25]
                high_q = d[d["F1"] >= q75]
                if q75 > q25 and len(low_q) >= 5 and len(high_q) >= 5:
                    effect_rows.append({"gene": gene, "panel": panel, "endpoint": "F1_patient_quartile",
                                        "method": "median_delta", "curated_patient_id": pid,
                                        "n_cells": len(d), "n_low": len(low_q), "n_high": len(high_q),
                                        "effect": float(high_q[value_col].median() - low_q[value_col].median())})
    effects = pd.DataFrame(effect_rows)
    effects.to_csv(part / "m17_patient_gene_f1_effects_v1.tsv", sep="\t", index=False, na_rep="NA")
    gene_summary = summarize(effects, ["gene", "panel", "endpoint", "method"])
    gene_summary = classify_summary(gene_summary, ["endpoint", "method"])
    gene_summary.to_csv(part / "m17_gene_f1_summary_v1.tsv", sep="\t", index=False, na_rep="NA")

    integrated_rows = []
    for gene in GENES:
        s = gene_summary[gene_summary["gene"] == gene]
        supported = s[s["statistical_support"]]
        n_support = len(supported)
        n_pos = int((supported["direction"] == "positive").sum())
        n_neg = int((supported["direction"] == "negative").sum())
        direction = "positive" if n_pos and not n_neg else "negative" if n_neg and not n_pos else "mixed" if n_pos and n_neg else "unsupported"
        tier = ("STRONG_ASSOCIATION" if n_support == 3 else
                "SECONDARY_ASSOCIATION" if n_support == 2 else
                "EXPLORATORY_ASSOCIATION" if n_support == 1 else "NO_SUPPORTED_ASSOCIATION")
        integrated_rows.append({"gene": gene, "panel": s["panel"].iloc[0],
                                "n_supported_endpoints_of_3": n_support,
                                "n_positive_supported_endpoints": n_pos,
                                "n_negative_supported_endpoints": n_neg,
                                "association_direction": direction,
                                "evidence_tier": tier})
    integrated = pd.DataFrame(integrated_rows)
    integrated.to_csv(part / "m17_gene_integrated_evidence_v1.tsv", sep="\t", index=False)

    pair_rows = []
    for gene_a, gene_b in combinations(GENES, 2):
        panel_a = next(k for k, v in PANELS.items() if gene_a in v)
        panel_b = next(k for k, v in PANELS.items() if gene_b in v)
        for pid, d in data.groupby("curated_patient_id", sort=False):
            a = d[f"{gene_a}_log1p_cp10k"]
            b = d[f"{gene_b}_log1p_cp10k"]
            if len(d) >= 20 and a.nunique() > 1 and b.nunique() > 1:
                pair_rows.append({"gene_a": gene_a, "gene_b": gene_b,
                                  "panel_a": panel_a, "panel_b": panel_b,
                                  "curated_patient_id": pid, "n_cells": len(d),
                                  "effect": float(stats.spearmanr(a, b).statistic)})
    pair_detail = pd.DataFrame(pair_rows)
    pair_detail.to_csv(part / "m17_patient_gene_pair_correlations_v1.tsv", sep="\t", index=False)
    pair_summary = summarize(pair_detail, ["gene_a", "gene_b", "panel_a", "panel_b"])
    pair_summary = classify_summary(pair_summary, []) if False else pair_summary
    pair_summary["BH_FDR"] = bh(pair_summary["wilcoxon_p_two_sided"])
    n = pair_summary["n_eligible_patients"]
    pair_summary["statistical_support"] = ((n >= 10) & (pair_summary["BH_FDR"] < 0.05)) | (
        n.between(3, 9) & (pair_summary["wilcoxon_p_two_sided"] < 0.05))
    pair_summary["reporting_tier"] = np.select([n >= 10, n.between(3, 9), n < 3],
                                                ["FDR_AT_LEAST_10", "NOMINAL_P_3_TO_9", "NOT_TESTED_LT3"],
                                                default="NOT_TESTED")
    pair_summary["direction"] = np.select([pair_summary["median_patient_effect"] > 0,
                                            pair_summary["median_patient_effect"] < 0],
                                           ["positive", "negative"], default="zero")
    pair_summary["supported_direction"] = np.where(pair_summary["statistical_support"],
                                                    pair_summary["direction"], "unsupported")
    pair_summary.to_csv(part / "m17_gene_pair_correlation_summary_v1.tsv", sep="\t", index=False, na_rep="NA")

    status = {"status": "PASS_M17_ANALYSIS_AWAITING_INDEPENDENT_VALIDATION",
              "n_cells": len(data), "n_genes": len(GENES), "n_gene_pairs": len(pair_summary),
              "n_patients_total": int(data["curated_patient_id"].nunique()),
              "direction_prespecified": False, "two_sided_testing": True,
              "figures_created": False, "python": sys.version.split()[0],
              "platform": platform.platform(), "numpy": np.__version__,
              "pandas": pd.__version__, "scipy": scipy.__version__, "h5py": h5py.__version__}
    (part / "status.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    os.replace(part, output)


if __name__ == "__main__":
    main()

