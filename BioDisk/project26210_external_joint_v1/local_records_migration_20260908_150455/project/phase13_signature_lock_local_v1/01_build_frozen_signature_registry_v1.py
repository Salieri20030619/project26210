#!/usr/bin/env python3
"""Build the immutable Phase 13 signature and analysis registry.

The program definitions are copied from pre-existing project scripts before any
external-cohort patient-level outcomes are calculated.  Outputs are written via
temporary ``.part`` files and atomically renamed.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
import sys
from datetime import datetime
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = HERE / "source_project_v1"
OUT = HERE / "registry_v1"
F1_SOURCE = (
    ROOT
    / "phase8_spatial_alignment_local_v1"
    / "frozen_reference_inputs_v1"
    / "cnmf_rank7"
    / "rank7_program_top100_genes_v1.tsv"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_text(path: Path, text: str) -> None:
    part = path.with_name(path.name + ".part")
    if path.exists() or part.exists():
        raise FileExistsError(f"Refusing existing output: {path}")
    with part.open("w", encoding="utf-8", newline="") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(part, path)


def atomic_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    part = path.with_name(path.name + ".part")
    if path.exists() or part.exists():
        raise FileExistsError(f"Refusing existing output: {path}")
    with part.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(part, path)


PROJECT_SIGNATURES = [
    {
        "signature_id": "tymp_expression",
        "display_name": "TYMP expression",
        "family": "malignant_program",
        "compartment": "malignant",
        "role": "primary_component",
        "genes": ["TYMP"],
        "source": "run_tymp_spp1_state_analysis_v1.py",
        "source_lines": "373-377",
    },
    {
        "signature_id": "spp1_state_no_spp1",
        "display_name": "SPP1 state without SPP1",
        "family": "myeloid_state",
        "compartment": "myeloid",
        "role": "primary",
        "genes": ["GPNMB", "LGALS3", "FABP5", "CD9", "LPL"],
        "source": "run_tymp_spp1_state_analysis_v1.py",
        "source_lines": "46-48",
    },
    {
        "signature_id": "spp1_state_full",
        "display_name": "SPP1 state including SPP1",
        "family": "myeloid_state",
        "compartment": "myeloid",
        "role": "sensitivity",
        "genes": ["SPP1", "GPNMB", "LGALS3", "FABP5", "CD9", "LPL"],
        "source": "run_tymp_spp1_state_analysis_v1.py",
        "source_lines": "46-47",
    },
    {
        "signature_id": "c1qc_state",
        "display_name": "C1QC state",
        "family": "myeloid_state",
        "compartment": "myeloid",
        "role": "primary",
        "genes": ["C1QA", "C1QB", "C1QC"],
        "source": "run_tymp_spp1_state_analysis_v1.py",
        "source_lines": "48-48",
    },
    {
        "signature_id": "apoe_tam_expression_proxy",
        "display_name": "APOE/C1QC TAM expression proxy",
        "family": "myeloid_state",
        "compartment": "myeloid",
        "role": "primary",
        "genes": ["APOE", "C1QA", "C1QB", "C1QC", "CD68", "CD163", "MRC1", "LGALS3", "CTSB", "LPL"],
        "source": "analyze_postn_apoe_singlecell_support_v1.py",
        "source_lines": "35-38",
    },
    {
        "signature_id": "monocyte_program",
        "display_name": "Monocyte program",
        "family": "myeloid_state",
        "compartment": "myeloid",
        "role": "supporting",
        "genes": ["FCN1", "VCAN", "S100A8", "S100A9", "CTSS", "LILRB1"],
        "source": "run_monocyte_spp1_trajectory_v1.py",
        "source_lines": "22-25",
    },
    {
        "signature_id": "myeloid_immunosuppression",
        "display_name": "Myeloid immunosuppression",
        "family": "myeloid_state",
        "compartment": "myeloid",
        "role": "primary",
        "genes": ["SPP1", "CD163", "MRC1", "IL10", "TGFB1", "TREM2", "APOE", "C1QC", "VSIG4", "MARCO"],
        "source": "analyze_tymp_immune_ecology_v2.py",
        "source_lines": "97-108",
    },
    {
        "signature_id": "immunosuppression_core",
        "display_name": "Immunosuppression core",
        "family": "myeloid_state",
        "compartment": "myeloid",
        "role": "sensitivity",
        "genes": ["IL10", "TGFB1", "VSIG4", "MARCO", "CD274", "IDO1"],
        "source": "run_monocyte_spp1_trajectory_v1.py",
        "source_lines": "28-28",
    },
    {
        "signature_id": "vascular_remodeling",
        "display_name": "Vascular remodeling",
        "family": "vascular_hypoxia",
        "compartment": "myeloid;endothelial",
        "role": "supporting",
        "genes": ["VEGFA", "MMP9", "MMP14", "ANGPTL4", "IL1B", "LGALS3"],
        "source": "run_monocyte_spp1_trajectory_v1.py",
        "source_lines": "27-27",
    },
    {
        "signature_id": "hypoxia_angiogenesis",
        "display_name": "Hypoxia-angiogenesis",
        "family": "vascular_hypoxia",
        "compartment": "malignant;endothelial",
        "role": "primary",
        "genes": ["CA9", "VEGFA", "HIF1A", "ADM", "ANGPTL4", "LDHA"],
        "source": "run_tymp_vascular_niche_analysis_v1.py",
        "source_lines": "40-44",
    },
    {
        "signature_id": "invasion_emt",
        "display_name": "Invasion-EMT",
        "family": "malignant_program",
        "compartment": "malignant",
        "role": "supporting",
        "genes": ["VIM", "FN1", "MMP2", "MMP9", "ITGA5", "AXL"],
        "source": "run_tymp_vascular_niche_analysis_v1.py",
        "source_lines": "40-44",
    },
    {
        "signature_id": "tcell_dysfunction",
        "display_name": "T-cell dysfunction",
        "family": "tcell_function",
        "compartment": "t_nk",
        "role": "primary",
        "genes": ["PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX", "CXCL13", "ENTPD1"],
        "source": "analyze_tymp_immune_ecology_v2.py",
        "source_lines": "85-96",
    },
    {
        "signature_id": "cytotoxicity",
        "display_name": "Cytotoxicity",
        "family": "tcell_function",
        "compartment": "t_nk",
        "role": "primary",
        "genes": ["NKG7", "GNLY", "GZMB", "PRF1", "CTSW"],
        "source": "analyze_tymp_immune_ecology_v2.py",
        "source_lines": "96-96",
    },
    {
        "signature_id": "ifn_response",
        "display_name": "IFN response",
        "family": "immune_context",
        "compartment": "myeloid;t_nk",
        "role": "supporting",
        "genes": ["IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "ISG15", "IDO1"],
        "source": "analyze_tymp_immune_ecology_v2.py",
        "source_lines": "109-109",
    },
    {
        "signature_id": "polyamine_mta_context",
        "display_name": "Polyamine-MTA context",
        "family": "metabolic_proxy",
        "compartment": "malignant;myeloid",
        "role": "primary",
        "genes": ["MAT2A", "MAT2B", "AMD1", "SRM", "SMS"],
        "source": "analyze_tymp_immune_ecology_v2.py",
        "source_lines": "111-123,317-325",
    },
    {
        "signature_id": "postn_ecm_proxy",
        "display_name": "POSTN ECM proxy",
        "family": "stromal_context",
        "compartment": "fibroblast",
        "role": "exploratory",
        "genes": ["POSTN", "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL5A2", "FN1", "FAP", "PDGFRA", "LOX", "LOXL2", "MMP2", "MMP14", "THBS1", "TGFBI"],
        "source": "analyze_postn_apoe_singlecell_support_v1.py",
        "source_lines": "27-32",
    },
]


# HLA-II was requested in the external-validation plan but no exact HLA-II gene
# set existed in the recovered project scripts.  It is therefore frozen here,
# before outcomes, as a supporting (not project-derived primary) module.
PRESPECIFIED_SUPPORTING = [
    {
        "signature_id": "hla_ii_supporting",
        "display_name": "HLA-II supporting module",
        "family": "myeloid_state",
        "compartment": "myeloid",
        "role": "supporting_prespecified_not_project_derived",
        "genes": ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA"],
        "source": "phase13_prespecification",
        "source_lines": "01_build_frozen_signature_registry_v1.py",
    }
]


def main() -> None:
    required = [SOURCE / item["source"] for item in PROJECT_SIGNATURES if item["source"] != "phase13_prespecification"]
    required.extend([F1_SOURCE])
    missing = sorted({str(path) for path in required if not path.is_file()})
    if missing:
        raise FileNotFoundError("Missing frozen sources:\n" + "\n".join(missing))

    OUT.mkdir(parents=False, exist_ok=False)

    with F1_SOURCE.open("r", encoding="utf-8", newline="") as handle:
        f1_all = list(csv.DictReader(handle, delimiter="\t"))
    f1 = [row for row in f1_all if row["factor"] == "F1" and int(row["gene_rank"]) <= 100]
    if len(f1) != 100 or len({row["gene"] for row in f1}) != 100:
        raise ValueError("F1 top-100 source is incomplete or duplicated")
    if any(float(row["loading"]) < 0 for row in f1):
        raise ValueError("F1 contains a negative loading")

    signatures = [
        {
            "signature_id": "f1_top100_weighted",
            "display_name": "F1 cNMF top-100 weighted program",
            "family": "malignant_program",
            "compartment": "malignant",
            "role": "primary",
            "genes": [row["gene"] for row in f1],
            "source": str(F1_SOURCE.relative_to(ROOT)).replace("\\", "/"),
            "source_lines": "factor=F1,gene_rank<=100",
            "weighted": True,
        },
        *PROJECT_SIGNATURES,
        *PRESPECIFIED_SUPPORTING,
    ]

    registry_rows = []
    gene_rows = []
    for item in signatures:
        weighted = bool(item.get("weighted", False))
        source_path = Path(item["source"])
        if not source_path.is_absolute() and item["source"] != "phase13_prespecification":
            source_path = (ROOT / source_path) if item["signature_id"] == "f1_top100_weighted" else (SOURCE / source_path)
        source_hash = "registry_self" if item["source"] == "phase13_prespecification" else sha256(source_path)
        minimum = 70 if item["signature_id"] == "f1_top100_weighted" else min(len(item["genes"]), max(1 if len(item["genes"]) == 1 else 2, (len(item["genes"]) * 3 + 4) // 5))
        scoring = (
            "sum(log1p(counts_per_10000)*loading)/sum(present_loadings)"
            if weighted
            else ("log1p(counts_per_10000)" if len(item["genes"]) == 1 else "mean(gene-wise population-z of log1p(counts_per_10000))")
        )
        registry_rows.append(
            {
                "signature_id": item["signature_id"],
                "display_name": item["display_name"],
                "family": item["family"],
                "compartment": item["compartment"],
                "role": item["role"],
                "n_genes": len(item["genes"]),
                "minimum_present_genes": minimum,
                "primary_scoring_rule": scoring,
                "sensitivity_scoring_rule": "mean(log1p(counts_per_10000)) for multi-gene modules; primary rule retained for single-gene TYMP",
                "source_file": item["source"],
                "source_lines_or_filter": item["source_lines"],
                "source_sha256": source_hash,
                "outcome_blind_lock": "yes",
            }
        )
        f1_by_gene = {row["gene"]: row for row in f1} if weighted else {}
        for order, gene in enumerate(item["genes"], start=1):
            gene_rows.append(
                {
                    "signature_id": item["signature_id"],
                    "gene_order": order,
                    "gene": gene,
                    "raw_weight": f1_by_gene.get(gene, {}).get("loading", "1"),
                    "source_gene_rank": f1_by_gene.get(gene, {}).get("gene_rank", order),
                    "direction": "positive",
                }
            )

    analysis_rules = {
        "registry_version": "phase13_signature_registry_v1",
        "locked_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "outcome_blind": True,
        "datasets_are_separate": True,
        "input_layer": "raw counts only; library-size normalize each cell to 10000 then log1p",
        "primary_unit": "dataset-qualified patient; paired biopsy for verified longitudinal contrasts",
        "compartments": {
            "malignant": "final refined malignant cells after annotation/CNV consensus; provisional scoring may run but confirmatory synthesis waits for CNV finalization",
            "myeloid": "Monocyte, macrophage and dendritic/myeloid states; SPP1-positive malignant cells are excluded",
            "endothelial": "Endothelial",
            "t_nk": "T cell and NK cell",
            "fibroblast": "Fibroblast/CAF candidate",
        },
        "minimum_cells": {"per_patient_compartment_primary": 20, "per_patient_compartment_sensitivity": 10},
        "f1": {
            "minimum_present_genes": 70,
            "weights": "non-negative source loadings renormalized over present genes",
            "high_threshold": "above pooled normal-epithelial 95th percentile within dataset",
            "threshold_sensitivity": "normal-epithelial 90th percentile",
            "phase8_quartile_label": "descriptive only; not used as evidence because its prevalence is definitional",
        },
        "patient_summaries": ["mean", "median", "fraction_detected", "cell_count"],
        "primary_association": "Spearman correlation across patients within each dataset; require at least 5 eligible patients",
        "uncertainty": "patient bootstrap 95% confidence interval with fixed seed 20260812; 9999 replicates when computationally feasible",
        "multiple_testing": "Benjamini-Hochberg within each declared hypothesis family",
        "longitudinal": "within-patient change from verified baseline to cycle2/EOT; no unverified timepoint inference",
        "sensitivity": [
            "unweighted mean log-normalized score",
            "minimum 10 versus 20 cells",
            "F1 normal-epithelial 90th versus 95th percentile",
            "leave-one-patient-out within eligible datasets",
            "leave-one-study-out for comparable study effects",
            "GSE178481 full versus previously included subset when membership can be proven",
            "GSE207493 old versus newly added patients when membership can be proven",
        ],
        "interpretation_limits": [
            "no cell-level pooled p-values",
            "MTA/polyamine scores are transcriptional proxies, not metabolites or flux",
            "cross-sectional myeloid ordering is not lineage",
            "association and colocalization are not mediation or causality",
            "HLA-II supporting module is prespecified in Phase13 but is not recovered from a project-defined signature",
        ],
    }

    atomic_tsv(OUT / "frozen_signature_registry_v1.tsv", list(registry_rows[0]), registry_rows)
    atomic_tsv(OUT / "frozen_signature_genes_v1.tsv", list(gene_rows[0]), gene_rows)
    atomic_text(OUT / "frozen_analysis_rules_v1.json", json.dumps(analysis_rules, ensure_ascii=False, indent=2) + "\n")

    source_paths = sorted({path for path in required})
    source_rows = [
        {
            "path_from_project_root": str(path.relative_to(ROOT)).replace("\\", "/"),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in source_paths
    ]
    atomic_tsv(OUT / "source_manifest_v1.tsv", list(source_rows[0]), source_rows)

    runtime = {
        "command": " ".join([sys.executable, *sys.argv]),
        "python": sys.version,
        "platform": platform.platform(),
        "cwd": str(Path.cwd()),
        "outputs": [
            "frozen_signature_registry_v1.tsv",
            "frozen_signature_genes_v1.tsv",
            "frozen_analysis_rules_v1.json",
            "source_manifest_v1.tsv",
        ],
    }
    atomic_text(OUT / "run_metadata_v1.json", json.dumps(runtime, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
