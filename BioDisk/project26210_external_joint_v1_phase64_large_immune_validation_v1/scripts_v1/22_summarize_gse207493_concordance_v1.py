#!/usr/bin/env python3
"""Compare reproduced broad cell composition with counts reported by the authors."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


AUTHOR = {
    "T cells": 30339,
    "Tumor cells": 19819,
    "Monocyte_macrophages": 16444,
    "NK cells": 15115,
    "Endothelial cells": 8720,
    "CAF": 7742,
    "Dendritic cells": 1728,
    "B cells": 1979,
    "Mast cells": 837,
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("celltype_summary", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    frame = pd.read_csv(args.celltype_summary, sep="\t")
    observed = dict(zip(frame["cell_type"], frame["n_cells"]))
    reproduced = {
        "T cells": sum(observed.get(x, 0) for x in [
            "Treg", "Exhausted_T", "Proliferating_T", "CD8_T", "CD4_T"
        ]),
        "Tumor cells": observed.get("Tumor", 0),
        "Monocyte_macrophages": observed.get("Myeloid", 0),
        "NK cells": observed.get("NK", 0),
        "Endothelial cells": observed.get("Endothelial", 0),
        "CAF": observed.get("CAF", 0),
        "Dendritic cells": observed.get("Dendritic", 0),
        "B cells": observed.get("B_cell", 0) + observed.get("Plasma", 0),
        "Mast cells": observed.get("Mast", 0),
    }
    total = sum(AUTHOR.values())
    rows = []
    for cell_type, author_count in AUTHOR.items():
        reproduced_count = int(reproduced[cell_type])
        rows.append({
            "cell_type": cell_type,
            "author_count": author_count,
            "reproduced_count": reproduced_count,
            "count_difference": reproduced_count - author_count,
            "author_fraction": author_count / total,
            "reproduced_fraction": reproduced_count / total,
            "percentage_point_difference": 100 * (reproduced_count - author_count) / total,
        })
    comparison = pd.DataFrame(rows)
    comparison.to_csv(
        args.output_dir / "gse207493_author_composition_concordance_v1.tsv",
        sep="\t", index=False,
    )

    immune_types = {
        "T cells", "Monocyte_macrophages", "NK cells", "Dendritic cells",
        "B cells", "Mast cells",
    }
    author_immune = comparison.loc[comparison.cell_type.isin(immune_types), "author_count"].sum()
    reproduced_immune = comparison.loc[
        comparison.cell_type.isin(immune_types), "reproduced_count"
    ].sum()
    metrics = pd.DataFrame([
        {"metric": "total_cells", "value": total},
        {"metric": "pearson_r_across_9_broad_types", "value": stats.pearsonr(
            comparison.author_count, comparison.reproduced_count
        ).statistic},
        {"metric": "spearman_rho_across_9_broad_types", "value": stats.spearmanr(
            comparison.author_count, comparison.reproduced_count
        ).statistic},
        {"metric": "mean_absolute_percentage_point_difference", "value": np.mean(
            np.abs(comparison.percentage_point_difference)
        )},
        {"metric": "author_immune_fraction", "value": author_immune / total},
        {"metric": "reproduced_immune_fraction", "value": reproduced_immune / total},
        {"metric": "immune_fraction_percentage_point_difference", "value": 100 * (
            reproduced_immune - author_immune
        ) / total},
    ])
    metrics.to_csv(
        args.output_dir / "gse207493_author_composition_metrics_v1.tsv",
        sep="\t", index=False,
    )
    print(comparison.to_string(index=False))
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
