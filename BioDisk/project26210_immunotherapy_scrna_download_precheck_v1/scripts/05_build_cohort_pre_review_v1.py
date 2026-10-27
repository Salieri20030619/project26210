#!/usr/bin/env python3
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

BASE = Path("/mnt/BioDisk/wangmingyuan/project26210_immunotherapy_scrna_download_precheck_v1")
META = BASE / "metadata"
PRE = BASE / "precheck"
OUT = PRE / "cohort_pre_review_v1.md"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    g264 = read_tsv(META / "GSE264586_patient_response_v1.tsv")
    g294 = read_tsv(META / "GSE294109_patient_metadata_v1.tsv")
    tenx = read_tsv(PRE / "tenx_sample_precheck_summary_v1.tsv")
    tenx_by_dataset = {row["dataset"]: row for row in tenx}

    g264_response = Counter(row["response_binary"] for row in g264)
    g294_response = Counter(row["response_binary"] for row in g294)
    g264_timepoints = Counter(tp for row in g264 for tp in row["timepoints"].split(";"))
    g264_longitudinal = sum(int(row["n_samples"]) > 1 for row in g264)
    g294_multiregion = sum(int(row["n_libraries"]) > 1 for row in g294)

    text = f"""# Project 26210：免疫治疗单细胞队列下载与预审

## 结论

GSE264586 与 GSE294109 均可作为第一批 RCC/ccRCC 免疫治疗单细胞验证队列。两者应分别建模并在患者层面比较方向，不建议直接合并原始细胞后做单一显著性检验，因为治疗方案、取样时点、取样部位以及多区域设计均不同。

## GSE264586

患者数：{len(g264)}；10x 文库数：{tenx_by_dataset['GSE264586']['libraries']}；10x 条形码总数：{tenx_by_dataset['GSE264586']['total_barcodes']}。

疗效分组：Responder {g264_response['Responder']} 例，Nonresponder {g264_response['Nonresponder']} 例。时点构成为 Baseline {g264_timepoints['Baseline']}、C2 {g264_timepoints['C2']}、EOT {g264_timepoints['EOT']} 个样本。

适用方向：基线预测、治疗前后动态、应答者与非应答者的细胞状态变化。主要限制是应答组仅 3 例，且活检器官存在异质性；推断应以患者为统计单位，并报告效应量和方向一致性，不能依赖细胞数制造伪重复。

## GSE294109

患者数：{len(g294)}；10x 文库数：{tenx_by_dataset['GSE294109']['libraries']}；10x 条形码总数：{tenx_by_dataset['GSE294109']['total_barcodes']}；多区域患者数：{g294_multiregion}。

ICB 疗效可评估患者：Responder {g294_response['Responder']} 例，Nonresponder {g294_response['Nonresponder']} 例；另有 untreated {g294_response['NA']} 例及 cabozantinib 非 ICB 主分析 {g294_response['Not_ICB_primary']} 例。

适用方向：患者层面的 ICB 应答关联、不同免疫方案间的共享方向、GSE264586 候选信号的独立复核。多区域文库必须先在患者内汇总，或在模型中设置患者随机效应；不能把 25 个文库当作 25 个独立患者。

## 建议的首轮分析顺序

1. 两队列分别质控、注释与患者级 pseudobulk，保留各自的疗效定义。
2. GSE264586 优先做 Baseline 3 比 9 的方向性筛查，再做 {g264_longitudinal} 例具有纵向样本患者的配对变化分析。
3. GSE294109 以 12 例 ICB 可评估患者为主分析，先在患者内合并多区域文库，再比较 9 比 3。
4. 对上一轮报告中的核心细胞群、基因集与通路逐队列复核，记录效应方向、置信区间及患者留一法稳定性。
5. 最终只把跨队列方向一致、且不由单一患者或单一区域驱动的结果列为稳健外部证据。

## 文件完整性判定

GSE264586：{tenx_by_dataset['GSE264586']['complete_libraries']}/{tenx_by_dataset['GSE264586']['libraries']} 个文库三件套完整，{tenx_by_dataset['GSE264586']['dimension_match_libraries']}/{tenx_by_dataset['GSE264586']['libraries']} 个矩阵维度与 features/barcodes 行数一致。

GSE294109：{tenx_by_dataset['GSE294109']['complete_libraries']}/{tenx_by_dataset['GSE294109']['libraries']} 个文库三件套完整，{tenx_by_dataset['GSE294109']['dimension_match_libraries']}/{tenx_by_dataset['GSE294109']['libraries']} 个矩阵维度与 features/barcodes 行数一致。
"""
    OUT.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
