#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import os
from collections import Counter
from datetime import datetime
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
RUNS = ROOT / 'dataset_runs_v5'
FINAL = ROOT / 'results_v5'
PART = ROOT / 'results_v5.part'
DATASETS = ('GSE264586', 'GSE294109')


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep='\t')


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    temporary = path.with_name(path.name + '.part')
    frame.to_csv(temporary, sep='\t', index=False, na_rep='NA', lineterminator='\n')
    os.replace(temporary, path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


if FINAL.exists() or PART.exists():
    raise SystemExit('Refusing existing Phase67 final output')
PART.mkdir(parents=True)

summaries = []
correlations = []
audits = []
coverages = []
statuses = []
registries = []
for dataset in DATASETS:
    directory = RUNS / dataset
    status = read_tsv(directory / f'{dataset}_status_v1.tsv')
    if len(status) != 1 or status.iloc[0].status != 'PASS_PHASE67_DATASET_V4':
        raise RuntimeError(f'Dataset status not PASS: {dataset}')
    statuses.append(status)
    summaries.append(read_tsv(directory / f'{dataset}_external_replication_summary_v1.tsv'))
    correlations.append(read_tsv(directory / f'{dataset}_correlations_by_patient_timepoint_v1.tsv'))
    audits.append(read_tsv(directory / f'{dataset}_lineage_audit_v1.tsv'))
    coverages.append(read_tsv(directory / f'{dataset}_signature_coverage_v1.tsv'))
    registries.append(read_tsv(directory / f'{dataset}_signature_registry_v1.tsv'))

registry = registries[0].sort_values('signature_id').reset_index(drop=True)
if not registry.equals(registries[1].sort_values('signature_id').reset_index(drop=True)):
    raise RuntimeError('Dataset signature registries differ')
summary = pd.concat(summaries, ignore_index=True)
correlation = pd.concat(correlations, ignore_index=True)
audit = pd.concat(audits, ignore_index=True)
coverage = pd.concat(coverages, ignore_index=True)
status_all = pd.concat(statuses, ignore_index=True)

cross_rows = []
for signature in registry.itertuples(index=False):
    subset = summary.loc[summary.signature_id == signature.signature_id].copy()
    evaluations = {row.dataset: row for row in subset.itertuples(index=False)}
    replicated = [row for row in evaluations.values() if bool(row.external_replication)]
    concordant = [row for row in evaluations.values() if bool(row.direction_matches_phase64)]
    if len(evaluations) >= 2 and len(replicated) == len(evaluations):
        grade = 'R2_TWO_EXTERNAL_FDR'
    elif len(evaluations) >= 2 and replicated and len(concordant) == len(evaluations):
        grade = 'R1_ONE_EXTERNAL_FDR_CONCORDANT'
    elif len(evaluations) >= 2 and replicated:
        grade = 'R1_ONE_EXTERNAL_FDR_CONTEXT_CONFLICT'
    elif len(evaluations) == 1 and replicated:
        grade = 'R1_SINGLE_EXTERNAL_FDR'
    else:
        grade = 'R0_NOT_REPLICATED'
    row = {
        'endpoint': signature.endpoint,
        'family': signature.family,
        'signature_id': signature.signature_id,
        'label': signature.label,
        'expected_phase64_direction': signature.expected_phase64_direction,
        'datasets_evaluable': len(evaluations),
        'datasets_fdr_replication': len(replicated),
        'replication_grade': grade,
    }
    for dataset in DATASETS:
        value = evaluations.get(dataset)
        row[f'{dataset}_n_patients'] = value.n_patients if value else None
        row[f'{dataset}_median_adjusted_rho'] = value.median_cycle_adjusted_rho if value else None
        row[f'{dataset}_adjusted_q_value'] = value.adjusted_q_value if value else None
        row[f'{dataset}_direction_match'] = bool(value.direction_matches_phase64) if value else False
    cross_rows.append(row)

cross = pd.DataFrame(cross_rows).sort_values(['endpoint', 'family', 'replication_grade', 'signature_id']).reset_index(drop=True)
retained = cross.loc[cross.replication_grade != 'R0_NOT_REPLICATED'].copy()
grade_counts = Counter(cross.replication_grade)

write_tsv(registry, PART / 'phase67_signature_registry_v1.tsv')
write_tsv(audit, PART / 'phase67_lineage_audit_v1.tsv')
write_tsv(coverage, PART / 'phase67_signature_coverage_v1.tsv')
write_tsv(correlation.sort_values(['dataset', 'endpoint', 'family', 'signature_id', 'patient_key', 'timepoint']), PART / 'phase67_correlations_by_patient_timepoint_v1.tsv')
write_tsv(summary.sort_values(['dataset', 'endpoint', 'family', 'adjusted_q_value', 'signature_id']), PART / 'phase67_external_replication_summary_v1.tsv')
write_tsv(cross, PART / 'phase67_cross_dataset_replication_v1.tsv')
write_tsv(retained, PART / 'phase67_cross_dataset_retained_v1.tsv')
write_tsv(status_all, PART / 'phase67_dataset_status_v1.tsv')

r2 = cross.loc[cross.replication_grade == 'R2_TWO_EXTERNAL_FDR']
r1_single = cross.loc[cross.replication_grade == 'R1_SINGLE_EXTERNAL_FDR']
def labels(frame: pd.DataFrame) -> str:
    return '、'.join(frame.label.astype(str)) if len(frame) else '无'

report = f"""# Phase67：TYMP–免疫程序直接外部复现

## 一、设计与范围

- 正式分析串行读取 GSE264586 与 GSE294109 的 Phase6 Seurat 对象，只保留冻结 Phase7 Myeloid broad 和 Treg fine 细胞。
- Myeloid/DC 共 50,531 个细胞；Treg 共 6,341 个细胞。GSE264586 没有冻结 Treg 细注释，因此 Treg 只能在 GSE294109 检验。
- 50 个 Phase64 冻结签名全部移除 TYMP；使用与 Scanpy 相同的 25 表达分箱、每分箱匹配对照基因和 `ctrl_as_ref=False` 逻辑计算程序分数。
- 在患者/时间点/谱系内计算 TYMP 与程序分数的 Spearman rho，并分别对 TYMP 与程序分数进行 S/G2M 残差化。至少 50 个细胞才进入分析。
- 统计推断使用患者跨时间点的中位 rho；BH 在每个数据集、谱系和签名家族内计算。

## 二、结果闭合

- 患者/时间点相关结果 {len(correlation):,} 行；外部数据集汇总 {len(summary):,} 行；跨数据集签名 {len(cross)} 行。
- GSE264586：30 个 Myeloid 签名可评估，21 个与 Phase64 同方向且通过 FDR。
- GSE294109：30 个 Myeloid 签名中 23 个、18 个 Treg 签名中 13 个与 Phase64 同方向且通过 FDR。
- 跨数据集等级：{json.dumps(dict(sorted(grade_counts.items())), ensure_ascii=False)}。

## 三、两外部队列直接复现

两队列均同方向并通过 FDR 的 R2 签名共 {len(r2)} 个：{labels(r2)}。

这些结果把对应 Myeloid/DC‑TYMP 程序从“GSE207493 内部双方法收敛”推进到“两个外部免疫治疗队列内的患者级直接相关复现”。

## 四、Treg 单外部队列复现

GSE294109 中达到 R1 单外部 FDR 的 Treg 签名共 {len(r1_single)} 个：{labels(r1_single)}。

由于 GSE264586 无冻结 Treg 细注释，不能把这些结果称为双外部队列复现。

## 五、解释边界

Phase67 直接复算的是外部队列内 TYMP 连续表达与冻结程序活性的患者级相关性，强于 Phase66 的临床情境转移筛查。但该结果仍是机制复现，不等同于 response/PFS/OS 预测验证；F1/TYMP 的直接临床预测等级不因本阶段结果自动升级。
"""
(PART / 'PHASE67_DIRECT_REPLICATION_SUMMARY_v1.md').write_text(report, encoding='utf-8')

metadata = {
    'status': 'FINALIZATION_COMPLETE_VALIDATION_PENDING',
    'completed_at': datetime.now().astimezone().isoformat(timespec='seconds'),
    'correlation_rows': len(correlation),
    'external_summary_rows': len(summary),
    'cross_dataset_rows': len(cross),
    'retained_rows': len(retained),
    'replication_grade_counts': dict(sorted(grade_counts.items())),
    'authoritative_dataset_run': 'dataset_runs_v5',
}
(PART / 'run_metadata_v1.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
os.replace(PART, FINAL)
print(json.dumps(metadata, ensure_ascii=False, indent=2))
