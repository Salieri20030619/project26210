#!/usr/bin/env python3
"""Minimal CPU training smoke test for the installed cell2location/scvi stack."""

import json

import anndata as ad
import numpy as np
import pandas as pd
import scvi
import torch
from cell2location.models import Cell2location, RegressionModel
from scipy import sparse


def main() -> None:
    rng = np.random.default_rng(26210)
    genes = pd.Index([f"GENE_{index:02d}" for index in range(24)])

    reference_counts = sparse.csr_matrix(rng.poisson(2.0, size=(40, len(genes))).astype(np.int32))
    reference = ad.AnnData(
        X=reference_counts.copy(),
        obs=pd.DataFrame(
            {
                "cell_type": pd.Categorical(["type_a"] * 20 + ["type_b"] * 20),
                "batch": pd.Categorical(["batch_1", "batch_2"] * 20),
            },
            index=[f"cell_{index:02d}" for index in range(40)],
        ),
        var=pd.DataFrame(index=genes),
    )
    reference.layers["counts"] = reference_counts.copy()
    scvi.settings.seed = 26210
    torch.set_num_threads(2)
    RegressionModel.setup_anndata(reference, layer="counts", batch_key="batch", labels_key="cell_type")
    reference_model = RegressionModel(reference)
    reference_model.train(
        max_epochs=1,
        batch_size=20,
        train_size=1,
        lr=0.002,
        accelerator="cpu",
        enable_progress_bar=False,
    )

    spatial_counts = sparse.csr_matrix(rng.poisson(3.0, size=(12, len(genes))).astype(np.int32))
    spatial = ad.AnnData(
        X=spatial_counts.copy(),
        obs=pd.DataFrame(
            {"section_batch": pd.Categorical(["section_1"] * 12)},
            index=[f"spot_{index:02d}" for index in range(12)],
        ),
        var=pd.DataFrame(index=genes),
    )
    spatial.layers["counts"] = spatial_counts.copy()
    signatures = pd.DataFrame(
        rng.uniform(0.5, 4.0, size=(len(genes), 2)),
        index=genes,
        columns=["type_a", "type_b"],
    )
    Cell2location.setup_anndata(spatial, layer="counts", batch_key="section_batch")
    spatial_model = Cell2location(
        spatial,
        cell_state_df=signatures,
        N_cells_per_location=10,
        detection_alpha=20,
    )
    spatial_model.train(
        max_epochs=1,
        batch_size=None,
        train_size=1,
        lr=0.002,
        accelerator="cpu",
        enable_progress_bar=False,
    )

    print(
        json.dumps(
            {
                "status": "PASS",
                "reference_train_epochs": 1,
                "spatial_train_epochs": 1,
                "device": "cpu",
                "scvi_tools": scvi.__version__,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
