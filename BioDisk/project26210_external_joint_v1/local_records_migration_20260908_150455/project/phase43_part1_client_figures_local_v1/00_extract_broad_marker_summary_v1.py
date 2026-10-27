from __future__ import annotations

import anndata as ad
import numpy as np


H5AD = (
    "/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/"
    "02_integration_annotation/cell_annotation_v1/"
    "part1_scvi_integrated_cell_annotated_v1.h5ad"
)

MARKERS = [
    "CD3D",
    "NKG7",
    "MS4A1",
    "MZB1",
    "C1QA",
    "S100A8",
    "CLEC9A",
    "TPSAB1",
    "EPCAM",
    "PECAM1",
    "RGS5",
    "COL1A1",
    "HBB",
]

CELL_TYPES = [
    "T cell",
    "NK cell",
    "B cell",
    "Plasma cell",
    "Macrophage",
    "Monocyte",
    "Dendritic cell",
    "Mast cell",
    "Epithelial cell",
    "Endothelial",
    "Pericyte/SMC",
    "Fibroblast",
    "Erythrocyte",
    "Epithelial/Myeloid mixed",
]


adata = ad.read_h5ad(H5AD, backed="r")
obs = adata.obs["cell_type_broad_v1"].astype(str).to_numpy()
gene_to_index = {gene: idx for idx, gene in enumerate(adata.var_names)}
marker_matrix = adata[:, [gene_to_index[marker] for marker in MARKERS]].X
if hasattr(marker_matrix, "toarray"):
    marker_matrix = marker_matrix.toarray()
marker_matrix = np.asarray(marker_matrix)

print("cell_type\tmarker\tmean_log1p\tfraction_expressing\tn_cells")
for cell_type in CELL_TYPES:
    cell_index = np.flatnonzero(obs == cell_type)
    for marker_index, marker in enumerate(MARKERS):
        values = marker_matrix[cell_index, marker_index]
        print(
            f"{cell_type}\t{marker}\t{values.mean():.8g}\t"
            f"{(values > 0).mean():.8g}\t{len(cell_index)}"
        )
