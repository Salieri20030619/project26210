#!/usr/bin/env python3
"""用 scFEA 官方 M168 网络推断恶性上皮细胞代谢通量。

上游 scFEA 代码依赖 pandas<2 的 ``DataFrame.append``，并把同一个基因矩阵
为 168 个模块重复展开。本实现保留官方模块、化学计量矩阵、神经网络结构、
超参数、损失项和随机种子，仅用逐模块基因遮罩替代显式重复矩阵。运行前会将
本实现和上游 ``FLUX`` 前向计算逐元素比较，确认数学等价。
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from torch import nn
from tqdm import tqdm


SEED = 16  # scFEA 官方代码固定的 torch seed
PROJECT_SEED = 26210
LEARN_RATE = 0.008
LAMB_BA = 1.0
LAMB_NG = 1.0
LAMB_CELL = 1.0
LAMB_MOD = 1e-2


class MaskedFlux(nn.Module):
    """与官方 FLUX 等价，但不显式构造 module x gene 的重复输入。"""

    def __init__(self, n_genes: int, n_modules: int, module_mask: torch.Tensor):
        super().__init__()
        self.register_buffer("module_mask", module_mask)
        self.m_encoder = nn.ModuleList(
            [
                nn.Sequential(
                    nn.Linear(n_genes, 8, bias=False),
                    nn.Tanhshrink(),
                    nn.Linear(8, 1),
                    nn.Tanhshrink(),
                )
                for _ in range(n_modules)
            ]
        )

    def forward(self, expression: torch.Tensor, stoichiometry: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        modules = [
            encoder(expression * self.module_mask[index])
            for index, encoder in enumerate(self.m_encoder)
        ]
        flux = torch.cat(modules, dim=1)
        balance = flux @ stoichiometry.T
        return flux, balance


def scfea_loss(
    flux: torch.Tensor,
    balance: torch.Tensor,
    gene_scale: torch.Tensor,
    module_scale: torch.Tensor,
) -> tuple[torch.Tensor, ...]:
    """官方 myLoss 的向量化等价实现，增加零方差数值保护。"""
    loss_balance = torch.sum(torch.sum(balance**2, dim=1) * LAMB_BA)
    loss_nonnegative = torch.sum(torch.sum(torch.abs(flux) - flux, dim=1) * LAMB_NG)

    difference = (torch.sum(flux, dim=1) - gene_scale) ** 2
    cell_penalty = torch.sqrt(difference) if torch.all(difference > 0) else difference
    loss_cell = torch.sum(cell_penalty * LAMB_CELL)

    flux_centered = flux - flux.mean(dim=1, keepdim=True)
    module_centered = module_scale - module_scale.mean(dim=1, keepdim=True)
    numerator = torch.sum(flux_centered * module_centered, dim=1)
    denominator = torch.linalg.vector_norm(flux_centered, dim=1) * torch.linalg.vector_norm(
        module_centered, dim=1
    )
    correlation = numerator / torch.clamp(denominator, min=1e-12)
    loss_module = torch.sum((1.0 - torch.abs(correlation)) * LAMB_MOD)

    total = loss_balance + loss_nonnegative + loss_cell + loss_module
    return total, loss_balance, loss_nonnegative, loss_cell, loss_module


def build_module_mask(
    module_gene: pd.DataFrame, gene_names: list[str]
) -> tuple[np.ndarray, np.ndarray]:
    gene_index = {gene: index for index, gene in enumerate(gene_names)}
    mask = np.zeros((module_gene.shape[0], len(gene_names)), dtype=np.float32)
    module_lengths = module_gene.notna().sum(axis=1).to_numpy(dtype=np.float32)
    for module_index, (_, row) in enumerate(module_gene.iterrows()):
        for gene in row.dropna().astype(str):
            if gene in gene_index:
                mask[module_index, gene_index[gene]] = 1.0
    return mask, module_lengths


def validate_forward_equivalence(
    model: MaskedFlux,
    expression: torch.Tensor,
    stoichiometry: torch.Tensor,
    source_dir: Path,
) -> float:
    """用同一权重比较官方显式展开前向传播和遮罩前向传播。"""
    sys.path.insert(0, str(source_dir / "src"))
    from ClassFlux import FLUX as UpstreamFlux  # type: ignore

    probe = expression[: min(8, expression.shape[0])].detach().cpu()
    mask = model.module_mask.detach().cpu()
    flattened = torch.cat([probe * mask[index] for index in range(mask.shape[0])], dim=1)
    upstream = UpstreamFlux(
        flattened,
        n_modules=mask.shape[0],
        f_in=probe.shape[1],
        f_out=1,
    ).cpu()
    upstream.load_state_dict(model.cpu().state_dict(), strict=False)
    with torch.no_grad():
        expected_flux, expected_balance = upstream(
            flattened,
            mask.shape[0],
            probe.shape[1],
            stoichiometry.shape[0],
            stoichiometry.cpu(),
        )
        observed_flux, observed_balance = model(probe, stoichiometry.cpu())
    difference = max(
        float(torch.max(torch.abs(expected_flux - observed_flux))),
        float(torch.max(torch.abs(expected_balance - observed_balance))),
    )
    if difference > 1e-5:
        raise RuntimeError(f"scFEA compatibility forward test failed: max abs diff={difference}")
    return difference


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/09_metabolism_trajectory/metabolism_v1"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    project = Path("/mnt/BioDisk/ludingyu/project26210")
    env_dir = project / "envs" / "part2_6_2_7_20260726_v1"
    source = env_dir / "sources" / "scFEA"
    out = args.out_dir
    input_file = out / "scFEA_input" / "malignant_all_cells_scFEA_M168_raw_counts.csv"
    flux_file = out / "source_data" / "scFEA_M168_cell_flux.tsv.gz"
    balance_file = out / "source_data" / "scFEA_M168_cell_metabolite_balance.tsv.gz"
    model_file = out / "objects" / "scFEA_M168_masked_compat_state.pt"
    qc_file = out / "qc" / "scFEA_M168_run_qc.json"
    loss_file = out / "qc" / "scFEA_M168_training_loss.tsv"
    figure_file = out / "figures" / "F_scFEA_training_loss"
    for path in (flux_file.parent, model_file.parent, qc_file.parent, figure_file.parent):
        path.mkdir(parents=True, exist_ok=True)

    if not input_file.exists():
        raise FileNotFoundError(input_file)
    torch.set_num_threads(min(args.threads, os.cpu_count() or 1))
    torch.manual_seed(SEED)
    np.random.seed(PROJECT_SEED)
    device = torch.device("cpu")

    started = time.time()
    print(f"Reading raw scFEA input: {input_file}", flush=True)
    gene_by_cell = pd.read_csv(input_file, index_col=0)
    raw_values = gene_by_cell.to_numpy()
    integer_fraction = float(np.mean(np.isclose(raw_values, np.round(raw_values))))
    if raw_values.min() < 0 or integer_fraction < 1.0:
        raise ValueError("Input is not a non-negative integer raw UMI matrix")

    expression_frame = gene_by_cell.T.astype(np.float32)
    # 与官方 scFEA 一致：只要最大值超过 50，所有值做 log2(count+1)。
    upstream_log2_applied = bool(expression_frame.to_numpy().max() > 50)
    if upstream_log2_applied:
        expression_frame = np.log2(expression_frame + 1.0)

    module_gene = pd.read_csv(source / "data" / "module_gene_m168.csv", index_col=0)
    stoichiometry_np = pd.read_csv(
        source / "data" / "cmMat_c70_m168.csv", header=None
    ).to_numpy(dtype=np.float32)
    compound_names = pd.read_csv(source / "data" / "cName_c70_m168.csv").columns.tolist()

    module_gene_union = sorted(set(module_gene.stack().astype(str)))
    gene_names = sorted(set(expression_frame.columns).intersection(module_gene_union))
    expression_frame = expression_frame.loc[:, gene_names]
    cell_names = expression_frame.index.astype(str).tolist()
    module_mask_np, module_lengths = build_module_mask(module_gene, gene_names)
    if np.any(module_mask_np.sum(axis=1) == 0):
        empty = module_gene.index[module_mask_np.sum(axis=1) == 0].tolist()
        raise ValueError(f"No input genes for scFEA modules: {empty}")

    expression = torch.tensor(expression_frame.to_numpy(), dtype=torch.float32, device=device)
    module_mask = torch.tensor(module_mask_np, dtype=torch.float32, device=device)
    stoichiometry = torch.tensor(stoichiometry_np, dtype=torch.float32, device=device)
    module_scale = (expression @ module_mask.T) / torch.tensor(
        module_lengths, dtype=torch.float32, device=device
    )
    expression_sum = expression.sum(dim=1)
    gene_scale = expression_sum / expression_sum.mean()

    torch.manual_seed(SEED)
    model = MaskedFlux(len(gene_names), module_gene.shape[0], module_mask).to(device)
    compatibility_diff = validate_forward_equivalence(
        model, expression, stoichiometry, source
    )
    model = model.to(device)
    print(f"Forward equivalence max abs diff: {compatibility_diff:.3g}", flush=True)

    optimizer = torch.optim.Adam(model.parameters(), lr=LEARN_RATE)
    history: list[dict[str, float]] = []
    model.train()
    for epoch in tqdm(range(1, args.epochs + 1), desc="scFEA M168 epochs"):
        flux, balance = model(expression, stoichiometry)
        losses = scfea_loss(flux, balance, gene_scale, module_scale)
        optimizer.zero_grad(set_to_none=True)
        losses[0].backward()
        optimizer.step()
        history.append(
            {
                "epoch": epoch,
                "total_loss": float(losses[0].detach()),
                "balance_loss": float(losses[1].detach()),
                "nonnegative_loss": float(losses[2].detach()),
                "cell_variation_loss": float(losses[3].detach()),
                "module_variation_loss": float(losses[4].detach()),
            }
        )

    model.eval()
    with torch.no_grad():
        flux, balance = model(expression, stoichiometry)
    flux_np = flux.cpu().numpy()
    balance_np = balance.cpu().numpy()

    flux_table = pd.DataFrame(flux_np, index=cell_names, columns=module_gene.index)
    flux_table.index.name = "cell"
    flux_table.to_csv(flux_file, sep="\t", compression="gzip")
    balance_table = pd.DataFrame(balance_np, index=cell_names, columns=compound_names)
    balance_table.index.name = "cell"
    balance_table.to_csv(balance_file, sep="\t", compression="gzip")
    pd.DataFrame(history).to_csv(loss_file, sep="\t", index=False)

    torch.save(
        {
            "state_dict": model.state_dict(),
            "gene_names": gene_names,
            "module_names": module_gene.index.tolist(),
            "module_mask": module_mask_np,
            "parameters": {
                "epochs": args.epochs,
                "learning_rate": LEARN_RATE,
                "lamb_balance": LAMB_BA,
                "lamb_nonnegative": LAMB_NG,
                "lamb_cell": LAMB_CELL,
                "lamb_module": LAMB_MOD,
                "scfea_seed": SEED,
            },
        },
        model_file,
    )

    history_frame = pd.DataFrame(history)
    fig, ax = plt.subplots(figsize=(4.2, 2.9))
    for column, color in (
        ("total_loss", "#222222"),
        ("balance_loss", "#C94755"),
        ("cell_variation_loss", "#3C78A8"),
        ("module_variation_loss", "#4B9A77"),
    ):
        ax.plot(history_frame["epoch"], history_frame[column], lw=0.75, label=column.replace("_loss", ""), color=color)
    ax.set_yscale("symlog")
    ax.set(xlabel="Epoch", ylabel="Loss", title="scFEA M168 training")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)
    ax.legend(frameon=False, fontsize=6.5, ncol=2)
    fig.tight_layout()
    fig.savefig(str(figure_file) + ".pdf")
    fig.savefig(str(figure_file) + ".png", dpi=600, facecolor="white")
    plt.close(fig)

    qc = {
        "n_cells": len(cell_names),
        "n_input_model_genes": len(gene_names),
        "n_modules": int(module_gene.shape[0]),
        "n_compounds": int(stoichiometry_np.shape[0]),
        "raw_input_integer_fraction": integer_fraction,
        "raw_input_max": float(raw_values.max()),
        "upstream_log2_count_plus_one_applied": upstream_log2_applied,
        "epochs": args.epochs,
        "forward_equivalence_max_abs_diff": compatibility_diff,
        "memory_safe_change": "module masks replace explicit 168-fold gene matrix expansion",
        "flux_min": float(flux_np.min()),
        "flux_max": float(flux_np.max()),
        "flux_nonfinite_values": int(np.size(flux_np) - np.isfinite(flux_np).sum()),
        "runtime_seconds": time.time() - started,
        "python": sys.version,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "torch": torch.__version__,
        "torch_threads": torch.get_num_threads(),
        "device": str(device),
        "project_seed": PROJECT_SEED,
        "scfea_upstream_seed": SEED,
    }
    qc_file.write_text(json.dumps(qc, indent=2), encoding="utf-8")
    print(json.dumps(qc, indent=2), flush=True)


if __name__ == "__main__":
    main()
