from __future__ import annotations

import hashlib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(r"C:\Users\31937\project26210_phase109_targeted_gap_analysis_v1")
ANALYSIS = ROOT / "analysis_v3"
OUT = ROOT / "final_synthesis_v1"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


OUT.mkdir(parents=True, exist_ok=False)
ab = pd.read_csv(ANALYSIS / "abundance_primary_results_v1.tsv", sep="\t")
cc = pd.read_csv(ANALYSIS / "cellchat_fine_state_specificity_summary_v1.tsv", sep="\t")
sp = pd.read_csv(ANALYSIS / "spatial_frozen_edge_synthesis_v1.tsv", sep="\t")

edges = pd.DataFrame(
    [
        {
            "edge_id": "E01",
            "edge": "tumour F1 to TYMP RNA",
            "phase109_test": "not retested",
            "phase109_result": "prior evidence retained",
            "integrated_status": "SUPPORTED_WITH_HETEROGENEITY",
            "allowed_claim": "F1 and TYMP RNA are associated in the original Part1 data, with cohort heterogeneity.",
        },
        {
            "edge_id": "E03",
            "edge": "TYMP to measured MTA",
            "phase109_test": "not testable in available D01 data",
            "phase109_result": "no direct measured MTA",
            "integrated_status": "CURRENT_DATA_DO_NOT_SUPPORT",
            "allowed_claim": "Public baseline data do not establish TYMP-driven MTA production; direct perturbational measurement remains client-owned.",
        },
        {
            "edge_id": "E04",
            "edge": "MTA to SPP1 induction",
            "phase109_test": "not testable",
            "phase109_result": "no MTA stimulation plus SPP1 readout",
            "integrated_status": "PENDING_CLIENT_VALIDATION",
            "allowed_claim": "Mechanistic compatibility only; direct induction requires MTA stimulation, SPP1 RNA/protein and blockade/rescue.",
        },
        {
            "edge_id": "E05",
            "edge": "tumour TYMP/F1 state to SPP1 macrophage state",
            "phase109_test": "20-patient abundance/polarization plus 14-patient fine-state communication",
            "phase109_result": "0/4 abundance endpoints supported; 0/2 malignant-SPP1 communication endpoints supported",
            "integrated_status": "CURRENT_DATA_DO_NOT_SUPPORT_GENERAL_EDGE",
            "allowed_claim": "The original Part1 data do not support a general TYMP/F1-driven SPP1-state abundance, polarization or tumour-SPP1 communication edge.",
        },
        {
            "edge_id": "E06",
            "edge": "SPP1 macrophage state to endothelial state",
            "phase109_test": "14-patient fine-state CellChat reaggregation",
            "phase109_result": "SPP1-to-endothelial supported: 14/14 positive, median 0.21065, q=0.00305",
            "integrated_status": "SUPPORTED_COMPUTATIONAL_COMPATIBILITY",
            "allowed_claim": "SPP1-state macrophages show reproducible endothelial-directed communication compatibility; this is not binding or causality.",
        },
        {
            "edge_id": "E07",
            "edge": "endothelial state to SPP1 macrophage state",
            "phase109_test": "14-patient fine-state CellChat reaggregation",
            "phase109_result": "endothelial-to-SPP1 supported: 13/14 positive, median 0.14642, q=0.00305",
            "integrated_status": "SUPPORTED_COMPUTATIONAL_COMPATIBILITY",
            "allowed_claim": "Endothelial cells show reproducible SPP1-state-directed communication compatibility; this is not binding or causality.",
        },
        {
            "edge_id": "E08",
            "edge": "SPP1 macrophage state to vascular remodeling",
            "phase109_test": "frozen Phase70 spatial synthesis; no new testing",
            "phase109_result": "supported in FFPE and frozen GSE175540 preparations",
            "integrated_status": "SUPPORTED_DATASET_SPECIFIC_SPATIAL",
            "allowed_claim": "The SPP1-state-to-vascular edge is supported within GSE175540 and triangulates with, but is not pooled with, CellChat.",
        },
        {
            "edge_id": "E09",
            "edge": "SPP1/vascular ecology to immune exclusion or clinical outcome",
            "phase109_test": "not retested",
            "phase109_result": "prior negative/context-dependent evidence retained",
            "integrated_status": "CURRENT_DATA_DO_NOT_SUPPORT_GENERAL_TRANSLATION",
            "allowed_claim": "No general immune-exclusion, ICB-response, PFS or OS conclusion is established.",
        },
    ]
)
edges.to_csv(OUT / "PHASE109_CHAIN_EDGE_UPDATE_v1.tsv", sep="\t", index=False)

summary = """# Phase109 targeted supplementation: decision summary

## Result

The targeted analysis did not close the unsupported upstream TYMP/F1-to-SPP1 edge. In 20 eligible patients from the original five-study Part1 universe, none of four preregistered tumour-axis-to-SPP1 abundance/polarization endpoints passed the support rule. The two malignant-to-SPP1 fine-state communication endpoints were also not supported.

The analysis did identify a reproducible downstream interface. In the 14 Phase35 CellChat patients, SPP1-state-to-endothelial specificity was positive in 14/14 patients (median 0.21065; BH q=0.00305), and endothelial-to-SPP1 specificity was positive in 13/14 patients (median 0.14642; BH q=0.00305). These are computational communication-compatibility results, not ligand binding, induction or causality.

Frozen Phase70 spatial evidence independently retains dataset-specific support for the SPP1-state-to-vascular-remodeling edge in both FFPE and frozen GSE175540 preparations. The CellChat and spatial results triangulate the downstream SPP1/endothelial-vascular interface without pooling samples or statistics.

## Consequence for the unified chain

- Retain F1-to-TYMP as supported with heterogeneity.
- Mark direct TYMP-to-measured-MTA as unsupported by current public baseline data.
- Keep MTA-to-SPP1 as pending client perturbational validation.
- Mark the general tumour TYMP/F1-to-SPP1-state edge as current-data-do-not-support.
- Add bidirectional SPP1-state/endothelial communication compatibility as supported computational evidence.
- Retain GSE175540 SPP1-state-to-vascular remodeling as dataset-specific spatial support.
- Do not upgrade immune exclusion or clinical translation.

## Validation

Formal analysis status: `PASS_ANALYSIS_COMPLETE`. Independent recomputation: `PASS_PHASE109_INDEPENDENT_VALIDATION`, 54/54 checks passed. Two pre-formal implementation failures are preserved in `analysis_v1` and `analysis_v2`; neither is an accepted result directory. `analysis_v3` is the sole formal result set.
"""
(OUT / "PHASE109_TARGETED_SUPPLEMENT_DECISION_v1.md").write_text(summary, encoding="utf-8")

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)

labels_a = [
    "TYMP axis -> SPP1 fraction",
    "TYMP axis -> SPP1:C1QC",
    "TYMP RNA -> SPP1 fraction",
    "TYMP RNA -> SPP1:C1QC",
]
y = np.arange(len(ab))[::-1]
colors_a = ["#b94a48" if not x else "#2f7d32" for x in ab["support"]]
axes[0].axvline(0, color="#777777", linewidth=1)
axes[0].scatter(ab["estimate"], y, s=70, c=colors_a, zorder=3)
for yy, est, q in zip(y, ab["estimate"], ab["bh_q_4"]):
    axes[0].text(est + (0.02 if est >= 0 else -0.02), yy + 0.16, f"q={q:.3f}", ha="left" if est >= 0 else "right", fontsize=8)
axes[0].set_yticks(y, labels_a)
axes[0].set_xlim(-0.35, 0.35)
axes[0].set_xlabel("Study-adjusted rank estimate")
axes[0].set_title("A. Upstream abundance bridge: 0/4 supported", loc="left", weight="bold")

labels_b = [
    "Malignant -> SPP1",
    "SPP1 -> malignant",
    "SPP1 -> endothelial",
    "Endothelial -> SPP1",
]
y2 = np.arange(len(cc))[::-1]
colors_b = ["#2f7d32" if x else "#b0b0b0" for x in cc["support"]]
axes[1].axvline(0, color="#777777", linewidth=1)
axes[1].scatter(cc["median"], y2, s=70, c=colors_b, zorder=3)
for yy, med, q, npos, n in zip(y2, cc["median"], cc["bh_q_4"], cc["n_positive"], cc["n_patients"]):
    axes[1].text(med + 0.008, yy + 0.14, f"{npos}/{n} positive; q={q:.3g}", ha="left", fontsize=8)
axes[1].set_yticks(y2, labels_b)
axes[1].set_xlim(-0.025, 0.285)
axes[1].set_xlabel("Median fine-state specificity")
axes[1].set_title("B. Fine-state communication: endothelial interface supported", loc="left", weight="bold")

fig.suptitle("Phase109 targeted gap analysis", fontsize=14, weight="bold")
fig.savefig(OUT / "PHASE109_TARGETED_GAP_SUMMARY_v1.png", dpi=220, bbox_inches="tight")
fig.savefig(OUT / "PHASE109_TARGETED_GAP_SUMMARY_v1.pdf", bbox_inches="tight")
plt.close(fig)

rows = []
for path in sorted(ROOT.rglob("*")):
    if not path.is_file() or "__pycache__" in path.parts or path.name == "FILE_MANIFEST_SHA256_v1.tsv":
        continue
    rows.append(
        {
            "relative_path": path.relative_to(ROOT).as_posix(),
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
    )
pd.DataFrame(rows).to_csv(ROOT / "FILE_MANIFEST_SHA256_v1.tsv", sep="\t", index=False)

print(f"PASS_PHASE109_FINAL_SYNTHESIS files_hashed={len(rows)}")
