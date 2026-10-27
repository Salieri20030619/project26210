#!/usr/bin/env python3
import argparse
import hashlib
import json
import math
import os
import platform
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "2")

import anndata as ad
import numpy as np
import pandas as pd
import scipy
from scipy import sparse, stats
from scipy.spatial import cKDTree

RESPONSE_GENES = ["ARG1", "VEGFA", "IL10", "TIMP1"]
VASCULAR_NO_VEGFA = ["MMP9", "MMP14", "ANGPTL4", "IL1B", "LGALS3"]

SCORE_COLUMNS = {
    "TYMP": "score_primary__tymp_expression",
    "SPP1_STATE": "score_primary__spp1_state_no_spp1",
    "T_DYSFUNCTION": "score_primary__tcell_dysfunction",
    "CYTOTOXICITY": "score_primary__cytotoxicity",
}

ABUNDANCE_GROUPS = {
    "Epithelial": ["Epithelial cell"],
    "Myeloid": ["Macrophage", "Monocyte"],
    "Vascular": ["Endothelial", "Pericyte/SMC"],
    "T_cell": ["T cell"],
    "Cytotoxic_lymphocyte": ["T cell", "NK cell"],
}

RELATIONS = [
    ("P01_TYMP_to_MTA_response", "TYMP", "MTA_RESPONSE", ["Epithelial", "Myeloid", "Vascular"], 1),
    ("P02_MTA_response_to_SPP1_state", "MTA_RESPONSE", "SPP1_STATE", ["Epithelial", "Myeloid", "Vascular"], 1),
    ("P03_TYMP_to_SPP1_state", "TYMP", "SPP1_STATE", ["Epithelial", "Myeloid"], 1),
    ("P04_MTA_response_to_vascular_no_VEGFA", "MTA_RESPONSE", "VASCULAR_NO_VEGFA", ["Epithelial", "Myeloid", "Vascular"], 1),
    ("P05_MTA_response_to_T_dysfunction", "MTA_RESPONSE", "T_DYSFUNCTION", ["Epithelial", "Myeloid", "T_cell"], 1),
    ("P06_MTA_response_to_cytotoxicity", "MTA_RESPONSE", "CYTOTOXICITY", ["Epithelial", "Myeloid", "Cytotoxic_lymphocyte"], -1),
]

CHAINS = [
    ("C01_TYMP_MTA_response_SPP1_state", "P01_TYMP_to_MTA_response", "P02_MTA_response_to_SPP1_state", 1),
    ("C02_TYMP_MTA_response_vascular_no_VEGFA", "P01_TYMP_to_MTA_response", "P04_MTA_response_to_vascular_no_VEGFA", 1),
    ("C03_TYMP_MTA_response_T_dysfunction", "P01_TYMP_to_MTA_response", "P05_MTA_response_to_T_dysfunction", 1),
    ("C04_TYMP_MTA_response_cytotoxicity", "P01_TYMP_to_MTA_response", "P06_MTA_response_to_cytotoxicity", -1),
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bh_adjust(values):
    values = np.asarray(values, dtype=float)
    out = np.full(len(values), np.nan)
    good = np.isfinite(values)
    if not good.any():
        return out
    vals = values[good]
    order = np.argsort(vals)
    ranked = vals[order] * len(vals) / np.arange(1, len(vals) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    restored = np.empty(len(vals))
    restored[order] = np.minimum(ranked, 1.0)
    out[np.where(good)[0]] = restored
    return out


def prep_for(sample_id):
    return "FFPE" if int(sample_id.replace("GSM", "")) <= 5924041 else "FROZEN"


def dense_vector(adata, gene):
    names = np.asarray([str(x).upper() for x in adata.var_names])
    idx = np.where(names == gene)[0]
    if len(idx) != 1:
        raise RuntimeError(f"Expected exactly one {gene}, found {len(idx)}")
    vec = adata.X[:, int(idx[0])]
    if sparse.issparse(vec):
        vec = vec.toarray()
    return np.asarray(vec, dtype=float).reshape(-1)


def zscore(values):
    sd = float(np.std(values, ddof=0))
    if not np.isfinite(sd) or sd == 0:
        return np.zeros(len(values), dtype=float), False
    return (values - float(np.mean(values))) / sd, True


def build_abundance(df):
    result = {}
    for name, members in ABUNDANCE_GROUPS.items():
        cols = [f"mean__{member}" for member in members]
        missing = [c for c in cols if c not in df.columns]
        if missing:
            raise RuntimeError(f"Missing abundance columns: {missing}")
        result[name] = df[cols].sum(axis=1).to_numpy(dtype=float)
    return result


def residual(values, covariates):
    ranked_y = stats.rankdata(values, method="average")
    ranked_cov = np.column_stack([stats.rankdata(x, method="average") for x in covariates])
    design = np.column_stack([np.ones(len(values)), ranked_cov])
    return ranked_y - design @ np.linalg.lstsq(design, ranked_y, rcond=None)[0]


def partial_spearman(x, y, covariates):
    matrix = np.column_stack([x, y] + covariates)
    keep = np.all(np.isfinite(matrix), axis=1)
    matrix = matrix[keep]
    if len(matrix) <= len(covariates) + 4:
        return math.nan, math.nan, int(len(matrix))
    rx = residual(matrix[:, 0], [matrix[:, i] for i in range(2, matrix.shape[1])])
    ry = residual(matrix[:, 1], [matrix[:, i] for i in range(2, matrix.shape[1])])
    if np.std(rx) == 0 or np.std(ry) == 0:
        return math.nan, math.nan, int(len(matrix))
    rho = float(np.corrcoef(rx, ry)[0, 1])
    dfree = len(matrix) - len(covariates) - 2
    tval = rho * math.sqrt(dfree / max(1e-15, 1.0 - rho * rho))
    return rho, float(2 * stats.t.sf(abs(tval), dfree)), int(len(matrix))


def neighbor_partial(x, y, current_covs, neighbors):
    neighbor_y = np.mean(y[neighbors], axis=1)
    neighbor_covs = [np.mean(cov[neighbors], axis=1) for cov in current_covs]
    matrix = np.column_stack([x, neighbor_y] + current_covs + neighbor_covs)
    keep = np.all(np.isfinite(matrix), axis=1)
    matrix = matrix[keep]
    k = len(current_covs)
    if len(matrix) <= 2 * k + 4:
        return math.nan, math.nan, int(len(matrix))
    rx = residual(matrix[:, 0], [matrix[:, i] for i in range(2, 2 + k)])
    ry = residual(matrix[:, 1], [matrix[:, i] for i in range(2 + k, 2 + 2 * k)])
    if np.std(rx) == 0 or np.std(ry) == 0:
        return math.nan, math.nan, int(len(matrix))
    rho = float(np.corrcoef(rx, ry)[0, 1])
    dfree = len(matrix) - 2 * k - 2
    tval = rho * math.sqrt(dfree / max(1e-15, 1.0 - rho * rho))
    return rho, float(2 * stats.t.sf(abs(tval), dfree)), int(len(matrix))


def meta_analyze(rows, id_col, family_col):
    output = []
    for keys, sub in rows.groupby(["preparation", family_col, id_col], sort=True):
        prep, family, relation = keys
        values = sub["rho"].dropna().to_numpy(dtype=float)
        expected = int(sub["expected_direction"].iloc[0])
        p = 1.0 if len(values) and np.allclose(values, 0) else float(stats.wilcoxon(values, alternative="two-sided", zero_method="wilcox").pvalue)
        output.append({
            "preparation": prep,
            family_col: family,
            id_col: relation,
            "n_sections": len(values),
            "median_rho": float(np.median(values)),
            "min_rho": float(np.min(values)),
            "max_rho": float(np.max(values)),
            "expected_direction": expected,
            "expected_direction_fraction": float(np.mean(values * expected > 0)),
            "wilcoxon_signed_rank_p": p,
        })
    meta = pd.DataFrame(output)
    meta["wilcoxon_bh_q"] = np.nan
    for _, idx in meta.groupby(["preparation", family_col], sort=True).groups.items():
        idx = list(idx)
        meta.loc[idx, "wilcoxon_bh_q"] = bh_adjust(meta.loc[idx, "wilcoxon_signed_rank_p"])
    meta["support_status"] = np.where(
        (meta["wilcoxon_bh_q"] < 0.05)
        & (meta["median_rho"] * meta["expected_direction"] > 0)
        & (meta["expected_direction_fraction"] >= 8 / 12),
        "SUPPORTED", "NOT_SUPPORTED"
    )
    return meta


def decisions_from_meta(meta, id_col, family_value, family_col):
    rows = []
    subset = meta[meta[family_col] == family_value]
    for relation, sub in subset.groupby(id_col, sort=True):
        status = dict(zip(sub["preparation"], sub["support_status"]))
        n = sum(status.get(x) == "SUPPORTED" for x in ["FFPE", "FROZEN"])
        rows.append({
            id_col: relation,
            "ffpe_status": status.get("FFPE", "MISSING"),
            "frozen_status": status.get("FROZEN", "MISSING"),
            "preparations_supported": n,
            "decision": "TWO_PREPARATION_SUPPORTED" if n == 2 else ("ONE_PREPARATION_SUPPORTED" if n == 1 else "NOT_SUPPORTED"),
        })
    return pd.DataFrame(rows)


def atomic_write(frame, path):
    tmp = Path(str(path) + ".part")
    frame.to_csv(tmp, sep="\t", index=False, na_rep="NA")
    os.replace(tmp, path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase54", required=True)
    parser.add_argument("--phase55", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    phase54 = Path(args.phase54).resolve()
    phase55 = Path(args.phase55).resolve()
    out = Path(args.out).resolve()
    part = Path(str(out) + ".part")
    if out.exists() or part.exists():
        raise RuntimeError("Refusing existing output or staging directory")
    part.mkdir(parents=True)

    pass54 = phase54 / "postrun_validation_v1" / "POSTRUN_VALIDATION_PASS_v1"
    pass55 = phase55 / "postrun_validation_v2" / "POSTRUN_VALIDATION_PASS_v2"
    if not pass54.is_file() or not pass55.is_file():
        raise RuntimeError("Frozen Phase54/55 PASS gate missing")

    sec54 = phase54 / "run_v1" / "sections_v1"
    sec55 = phase55 / "run_v1" / "sections_v1"
    ids54 = {p.name for p in sec54.iterdir() if p.is_dir()}
    ids55 = {p.name for p in sec55.iterdir() if p.is_dir()}
    if ids54 != ids55 or len(ids54) != 24:
        raise RuntimeError(f"Section mismatch Phase54={len(ids54)} Phase55={len(ids55)}")

    coverage_rows, relation_rows, manifest_rows, section_rows = [], [], [], []
    for sample_id in sorted(ids54):
        h5ad_path = sec54 / sample_id / "processed_section_v1.h5ad"
        score_path = sec54 / sample_id / "spot_metadata_scores_v1.tsv.gz"
        abundance_path = sec55 / sample_id / "spot_cell_abundance_v1.tsv.gz"
        for kind, path in [("phase54_h5ad", h5ad_path), ("phase54_scores", score_path), ("phase55_abundance", abundance_path)]:
            manifest_rows.append({"sample_id": sample_id, "input_type": kind, "path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)})

        adata = ad.read_h5ad(h5ad_path)
        scores = pd.read_csv(score_path, sep="\t")
        abundance = pd.read_csv(abundance_path, sep="\t")
        scores = scores[scores["qc_pass"].astype(str).str.lower().isin(["true", "1"])].copy()
        if scores["barcode"].duplicated().any() or abundance["barcode"].duplicated().any():
            raise RuntimeError(f"{sample_id}: duplicated barcode")
        merged = scores.merge(abundance, on="barcode", how="inner", validate="one_to_one")
        if set(adata.obs_names) != set(scores["barcode"]) or set(scores["barcode"]) != set(abundance["barcode"]):
            raise RuntimeError(f"{sample_id}: barcode mismatch")
        order = pd.Index(adata.obs_names).get_indexer(merged["barcode"])
        if (order < 0).any():
            raise RuntimeError(f"{sample_id}: H5AD join failure")

        response_z, response_raw, variable_response = [], [], []
        for gene in RESPONSE_GENES:
            raw = dense_vector(adata, gene)[order]
            z, variable = zscore(raw)
            response_raw.append(raw)
            response_z.append(z)
            variable_response.append(variable)
            coverage_rows.append({"sample_id": sample_id, "program": "MTA_RESPONSE", "gene": gene, "present": True, "nonzero_variance": variable})
        if sum(variable_response) < 2:
            raise RuntimeError(f"{sample_id}: fewer than two variable response genes")
        mta_response = np.mean(np.column_stack(response_z), axis=1)
        mta_response_raw = np.mean(np.column_stack(response_raw), axis=1)

        vascular_z, variable_vascular = [], []
        for gene in VASCULAR_NO_VEGFA:
            raw = dense_vector(adata, gene)[order]
            z, variable = zscore(raw)
            vascular_z.append(z)
            variable_vascular.append(variable)
            coverage_rows.append({"sample_id": sample_id, "program": "VASCULAR_NO_VEGFA", "gene": gene, "present": True, "nonzero_variance": variable})
        if sum(variable_vascular) < 3:
            raise RuntimeError(f"{sample_id}: fewer than three variable vascular genes")
        vascular_no_vegfa = np.mean(np.column_stack(vascular_z), axis=1)

        values = {name: merged[col].to_numpy(dtype=float) for name, col in SCORE_COLUMNS.items()}
        values["MTA_RESPONSE"] = mta_response
        values["MTA_RESPONSE_RAW"] = mta_response_raw
        values["VASCULAR_NO_VEGFA"] = vascular_no_vegfa
        abund = build_abundance(merged)
        log_umi = np.log1p(merged["total_counts"].to_numpy(dtype=float))
        coords = np.asarray(adata.obsm["spatial"])[order]
        _, neighbors = cKDTree(coords).query(coords, k=7)
        neighbors = np.asarray(neighbors[:, 1:], dtype=int)
        prep = prep_for(sample_id)

        section_rows.append({
            "sample_id": sample_id, "preparation": prep, "n_spots": len(merged),
            "exact_barcode_match": "PASS", "variable_response_genes": int(sum(variable_response)),
            "variable_vascular_no_vegfa_genes": int(sum(variable_vascular)),
        })

        for relation_id, left, right, cov_names, expected in RELATIONS:
            covs = [abund[x] for x in cov_names] + [log_umi]
            rho, p, n = partial_spearman(values[left], values[right], covs)
            relation_rows.append({
                "sample_id": sample_id, "preparation": prep, "analysis_family": "PRIMARY_SAME_SPOT",
                "relation_id": relation_id, "predictor": left, "outcome": right,
                "covariates": ";".join(cov_names + ["log1p_total_counts"]),
                "expected_direction": expected, "n_spots": n, "rho": rho, "spot_level_p_descriptive": p,
            })
            nrho, npval, nn = neighbor_partial(values[left], values[right], covs, neighbors)
            relation_rows.append({
                "sample_id": sample_id, "preparation": prep, "analysis_family": "SUPPORTING_K6_NEIGHBOR",
                "relation_id": relation_id, "predictor": left, "outcome": right,
                "covariates": ";".join(cov_names + ["log1p_total_counts"]),
                "expected_direction": expected, "n_spots": nn, "rho": nrho, "spot_level_p_descriptive": npval,
            })
        del adata

    sections = pd.DataFrame(section_rows)
    relations = pd.DataFrame(relation_rows)
    coverage = pd.DataFrame(coverage_rows)
    relation_meta = meta_analyze(relations, "relation_id", "analysis_family")
    primary_decisions = decisions_from_meta(relation_meta, "relation_id", "PRIMARY_SAME_SPOT", "analysis_family")
    neighbor_decisions = decisions_from_meta(relation_meta, "relation_id", "SUPPORTING_K6_NEIGHBOR", "analysis_family")

    chain_rows = []
    primary = relations[relations["analysis_family"] == "PRIMARY_SAME_SPOT"]
    for sample_id, sub in primary.groupby("sample_id", sort=True):
        lookup = dict(zip(sub["relation_id"], sub["rho"]))
        prep = sub["preparation"].iloc[0]
        for chain_id, first, second, expected in CHAINS:
            chain_rows.append({
                "sample_id": sample_id, "preparation": prep, "analysis_family": "PRIMARY_CHAIN_PRODUCT",
                "chain_id": chain_id, "first_relation": first, "second_relation": second,
                "expected_direction": expected, "rho": float(lookup[first] * lookup[second]),
            })
    chains = pd.DataFrame(chain_rows)
    chain_meta = meta_analyze(chains, "chain_id", "analysis_family")
    chain_decisions = decisions_from_meta(chain_meta, "chain_id", "PRIMARY_CHAIN_PRODUCT", "analysis_family")

    gates = [
        ("sections", len(sections), 24),
        ("ffpe_sections", int((sections["preparation"] == "FFPE").sum()), 12),
        ("frozen_sections", int((sections["preparation"] == "FROZEN").sum()), 12),
        ("barcode_pass", int((sections["exact_barcode_match"] == "PASS").sum()), 24),
        ("coverage_rows", len(coverage), 24 * 9),
        ("relation_rows", len(relations), 24 * 6 * 2),
        ("finite_relation_rho", int(np.isfinite(relations["rho"]).sum()), 24 * 6 * 2),
        ("relation_meta_rows", len(relation_meta), 2 * 2 * 6),
        ("primary_decisions", len(primary_decisions), 6),
        ("neighbor_decisions", len(neighbor_decisions), 6),
        ("chain_rows", len(chains), 24 * 4),
        ("finite_chain_products", int(np.isfinite(chains["rho"]).sum()), 24 * 4),
        ("chain_meta_rows", len(chain_meta), 2 * 4),
        ("chain_decisions", len(chain_decisions), 4),
    ]
    status = pd.DataFrame(gates, columns=["gate", "observed", "expected"])
    status["status"] = np.where(status["observed"] == status["expected"], "PASS", "FAIL")
    if (status["status"] != "PASS").any():
        raise RuntimeError("Phase75 gates failed")

    atomic_write(sections, part / "section_qc_v1.tsv")
    atomic_write(coverage, part / "program_gene_coverage_v1.tsv")
    atomic_write(relations, part / "relations_by_section_v1.tsv")
    atomic_write(relation_meta, part / "relation_meta_v1.tsv")
    atomic_write(primary_decisions, part / "primary_relation_decisions_v1.tsv")
    atomic_write(neighbor_decisions, part / "neighbor_relation_decisions_v1.tsv")
    atomic_write(chains, part / "chain_products_by_section_v1.tsv")
    atomic_write(chain_meta, part / "chain_product_meta_v1.tsv")
    atomic_write(chain_decisions, part / "chain_product_decisions_v1.tsv")
    atomic_write(status, part / "PHASE75_RUN_STATUS_v1.tsv")
    atomic_write(pd.DataFrame(manifest_rows), part / "INPUT_MANIFEST_SHA256_v1.tsv")

    two_primary = primary_decisions.loc[primary_decisions["decision"] == "TWO_PREPARATION_SUPPORTED", "relation_id"].tolist()
    one_primary = primary_decisions.loc[primary_decisions["decision"] == "ONE_PREPARATION_SUPPORTED", "relation_id"].tolist()
    two_chain = chain_decisions.loc[chain_decisions["decision"] == "TWO_PREPARATION_SUPPORTED", "chain_id"].tolist()
    summary = [
        "# Phase75 MTA-response spatial transfer",
        "",
        f"- Sections: {len(sections)}/24; spots: {int(sections['n_spots'].sum()):,}.",
        f"- Two-preparation primary relations: {len(two_primary)}/6.",
        f"- One-preparation primary relations: {len(one_primary)}/6.",
        f"- Two-preparation chain-consistency products: {len(two_chain)}/4.",
        "",
        "## Two-preparation primary relations", "",
        *([f"- {x}" for x in two_primary] if two_primary else ["- None."]),
        "", "## One-preparation primary relations", "",
        *([f"- {x}" for x in one_primary] if one_primary else ["- None."]),
        "", "## Two-preparation chain products", "",
        *([f"- {x}" for x in two_chain] if two_chain else ["- None."]),
        "", "## Boundary", "",
        "The transferred score is a transcriptomic response proxy, not measured 5'-MTA. Chain products are directional consistency statistics, not causal indirect effects. Sections are the inferential units and FFPE/frozen are analyzed separately.",
    ]
    (part / "PHASE75_SPATIAL_TRANSFER_SUMMARY_v1.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    metadata = {
        "status": "PASS", "phase54": str(phase54), "phase55": str(phase55),
        "sections": len(sections), "spots": int(sections["n_spots"].sum()),
        "python": platform.python_version(), "anndata": ad.__version__, "numpy": np.__version__,
        "pandas": pd.__version__, "scipy": scipy.__version__, "inference_unit": "section",
        "preparations_separate": True, "mta_measured": False,
    }
    (part / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    os.replace(part, out)
    print(json.dumps({"status": "PASS_PHASE75", "sections": len(sections), "spots": int(sections["n_spots"].sum())}))


if __name__ == "__main__":
    main()

