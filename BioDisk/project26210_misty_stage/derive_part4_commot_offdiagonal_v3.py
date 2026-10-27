#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


ROOT = Path("/mnt/BioDisk/ludingyu/project26210/results/part4_spatial")
OUTPUT_ROOT = ROOT / "advanced_patch_commot_offdiagonal_v3"
EXPECTED_CELLTYPES = {
    "APOE+C1QC+ TAM", "B cell", "Dendritic cell", "Endothelial", "Epithelial-Normal",
    "Mast cell", "Monocyte", "NK cell", "Other Macrophage", "SMC/Pericyte",
    "SPP1+ TAM", "T cell", "TYMP-High Malignant", "TYMP-Low Malignant",
}


def main(sample: str) -> None:
    if sample not in {"GLMF1", "GLMF2"}:
        raise ValueError("sample must be GLMF1 or GLMF2")
    source_dir = ROOT / f"advanced_patch_commot_{sample.lower()}_v2"
    output_dir = OUTPUT_ROOT / sample
    output_dir.mkdir(parents=True, exist_ok=False)
    manifest = pd.read_csv(source_dir / f"{sample}_commot_candidate_manifest_v2.tsv", sep="\t")
    full_summary = pd.read_csv(source_dir / f"{sample}_commot_transport_summary_v2.tsv", sep="\t")
    model = ad.read_h5ad(source_dir / f"{sample}_commot_candidate_v2.h5ad")
    analyzed = ad.read_h5ad(ROOT / "analysis" / f"{sample}_analyzed.h5ad")
    analyzed = analyzed[model.obs_names].copy()
    abundance = analyzed.obsm["q05_cell_abundance_w_sf"]
    if not isinstance(abundance, pd.DataFrame):
        abundance = pd.DataFrame(abundance, index=analyzed.obs_names, columns=analyzed.uns["mod"]["factor_names"])
    abundance = abundance.rename(columns=lambda value: str(value).split("_fg_", 1)[-1])
    if len(abundance.columns) != len(set(abundance.columns)) or set(abundance.columns) != EXPECTED_CELLTYPES:
        raise RuntimeError("Unexpected cell-abundance columns")
    values = abundance.to_numpy(dtype=float)
    if not np.all(np.isfinite(values)) or np.any(values < 0) or np.any(values.sum(axis=1) <= 0):
        raise RuntimeError("Invalid cell-abundance matrix")
    fractions = values / values.sum(axis=1, keepdims=True)
    celltypes = abundance.columns.tolist()

    transport_rows = []
    celltype_rows = []
    hypothesis_rows = []
    spot = pd.DataFrame({"barcode": model.obs_names})
    specifications = [
        ("secreted_200um", "candidate200", manifest[manifest["signaling_type"] == "Secreted Signaling"]),
        ("secreted_500um", "candidate500", manifest[manifest["signaling_type"] == "Secreted Signaling"]),
        ("contact_110um", "contact110", manifest[manifest["signaling_type"] == "Cell-Cell Contact"]),
    ]
    for scale, database_name, candidates in specifications:
        for candidate in candidates.to_dict(orient="records"):
            key = f"commot-{database_name}-{candidate['ligand']}-{candidate['receptor']}"
            matrix = model.obsp[key].tocsr(copy=True)
            full_total = float(matrix.sum())
            diagonal_score = float(matrix.diagonal().sum())
            matrix.setdiag(0)
            matrix.eliminate_zeros()
            offdiagonal_score = float(matrix.sum())
            expected = full_summary[(full_summary["scale"] == scale) & (full_summary["axis_id"] == candidate["axis_id"])]
            if len(expected) != 1 or not np.isclose(offdiagonal_score, expected.iloc[0]["offdiagonal_score"], rtol=1e-10, atol=1e-10):
                raise RuntimeError(f"Off-diagonal score mismatch for {sample} {scale} {candidate['axis_id']}")
            transport_rows.append({
                "sample": sample, "scale": scale, "axis_id": candidate["axis_id"], "ligand": candidate["ligand"],
                "receptor": candidate["receptor"], "full_total_score": full_total, "diagonal_score": diagonal_score,
                "diagonal_share": diagonal_score / full_total, "offdiagonal_score": offdiagonal_score,
                "offdiagonal_share": offdiagonal_score / full_total, "offdiagonal_edges": int(matrix.nnz),
                "interpretation_limit": "spot-to-spot transport after diagonal removal; mixed spots and spatial proximity remain non-cell-resolved",
            })
            sender = np.asarray(matrix.sum(axis=1)).ravel()
            receiver = np.asarray(matrix.sum(axis=0)).ravel()
            spot[f"sender__{scale}__{candidate['axis_id']}"] = sender
            spot[f"receiver__{scale}__{candidate['axis_id']}"] = receiver

            mass = fractions.T @ matrix.dot(fractions)
            mass_total = float(mass.sum())
            if mass_total <= 0:
                raise RuntimeError(f"No off-diagonal mass for {sample} {scale} {candidate['axis_id']}")
            mass_fraction = mass / mass_total
            source_marginal = mass_fraction.sum(axis=1)
            target_marginal = mass_fraction.sum(axis=0)
            baseline = np.outer(source_marginal, target_marginal)
            enrichment = np.divide(mass_fraction, baseline, out=np.full_like(mass_fraction, np.nan), where=baseline > 0)
            ranks = pd.Series(mass_fraction.ravel()).rank(method="min", ascending=False).to_numpy().reshape(mass_fraction.shape)
            for source_index, source_celltype in enumerate(celltypes):
                for target_index, target_celltype in enumerate(celltypes):
                    enrich = float(enrichment[source_index, target_index])
                    celltype_rows.append({
                        "sample": sample, "scale": scale, "axis_id": candidate["axis_id"], "ligand": candidate["ligand"],
                        "receptor": candidate["receptor"], "source_celltype": source_celltype, "target_celltype": target_celltype,
                        "offdiagonal_weighted_mass": float(mass[source_index, target_index]),
                        "offdiagonal_weighted_mass_fraction": float(mass_fraction[source_index, target_index]),
                        "offdiagonal_expected_fraction": float(baseline[source_index, target_index]),
                        "offdiagonal_proxy_enrichment": enrich,
                        "offdiagonal_proxy_log2_enrichment": float(np.log2(enrich)) if enrich > 0 else math.nan,
                        "offdiagonal_mass_rank_among_196": int(ranks[source_index, target_index]),
                        "interpretation_limit": "composition-weighted off-diagonal mixed-spot proxy; not cell-resolved communication",
                    })
            source_index = celltypes.index(candidate["source_proxy"])
            target_index = celltypes.index(candidate["target_proxy"])
            enrich = float(enrichment[source_index, target_index])
            hypothesis_rows.append({
                "sample": sample, "scale": scale, "axis_id": candidate["axis_id"], "ligand": candidate["ligand"],
                "receptor": candidate["receptor"], "source_proxy": candidate["source_proxy"], "target_proxy": candidate["target_proxy"],
                "offdiagonal_proxy_weighted_mass_fraction": float(mass_fraction[source_index, target_index]),
                "offdiagonal_proxy_expected_fraction": float(baseline[source_index, target_index]),
                "offdiagonal_proxy_enrichment": enrich,
                "offdiagonal_proxy_log2_enrichment": float(np.log2(enrich)) if enrich > 0 else math.nan,
                "offdiagonal_proxy_mass_rank_among_196": int(ranks[source_index, target_index]),
                "interpretation_limit": "prespecified composition-weighted off-diagonal mixed-spot proxy; not cell-resolved attribution or causality",
            })

    transport = pd.DataFrame(transport_rows)
    celltype = pd.DataFrame(celltype_rows)
    hypothesis = pd.DataFrame(hypothesis_rows)
    transport.to_csv(output_dir / f"{sample}_commot_offdiagonal_transport_v3.tsv", sep="\t", index=False)
    celltype.to_csv(output_dir / f"{sample}_commot_offdiagonal_celltype_edge_mass_proxy_v3.tsv", sep="\t", index=False)
    hypothesis.to_csv(output_dir / f"{sample}_commot_offdiagonal_hypothesis_proxy_support_v3.tsv", sep="\t", index=False)
    spot.to_csv(output_dir / f"{sample}_commot_offdiagonal_spot_sender_receiver_v3.tsv", sep="\t", index=False)
    metadata = {
        "status": "complete", "sample": sample, "source_commot_version": "v2", "n_axes_scales": len(transport),
        "n_celltype_proxy_rows": len(celltype), "diagonal_operation": "set diagonal to zero, eliminate zeros, then recompute marginals",
        "reason": "primary spot-to-spot interpretation must not be dominated by within-spot self transport",
        "interpretation_limit": "off-diagonal targeted transport remains a mixed-spot proximity proxy, not direct binding or causal cell-cell communication",
    }
    (output_dir / "run_metadata_v3.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (output_dir / "status.json").write_text(json.dumps({"status": "complete", "sample": sample}, indent=2) + "\n")
    print(f"COMMOT_OFFDIAGONAL_V3_ALL_PASS sample={sample} rows={len(transport)}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: derive_part4_commot_offdiagonal_v3.py GLMF1|GLMF2")
    main(sys.argv[1])
