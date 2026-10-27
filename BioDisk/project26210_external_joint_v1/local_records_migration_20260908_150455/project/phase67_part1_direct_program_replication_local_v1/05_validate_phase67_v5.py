#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / 'results_v5'


def read_tsv(name: str) -> pd.DataFrame:
    return pd.read_csv(RESULTS / name, sep='\t')


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


checks = []
def check(name: str, condition: bool, detail: str) -> None:
    checks.append({'check': name, 'status': 'PASS' if condition else 'FAIL', 'detail': detail})


status = read_tsv('phase67_dataset_status_v1.tsv')
check('dataset_statuses', len(status) == 2 and set(status.status) == {'PASS_PHASE67_DATASET_V4'}, str(status[['dataset', 'status']].to_dict('records')))
check('selected_cells', int(status.selected_cells.sum()) == 56872, str(int(status.selected_cells.sum())))
check('dataset_row_counts', dict(zip(status.dataset, status.correlation_rows.astype(int))) == {'GSE264586': 570, 'GSE294109': 588}, str(dict(zip(status.dataset, status.correlation_rows))))

registry = read_tsv('phase67_signature_registry_v1.tsv')
check('signature_registry', len(registry) == 50 and registry.signature_id.nunique() == 50, str(len(registry)))
coverage = read_tsv('phase67_signature_coverage_v1.tsv')
check('coverage_rows', len(coverage) == 82, str(len(coverage)))
eligible_coverage = coverage.loc[coverage.eligible.astype(bool)]
skipped_coverage = coverage.loc[~coverage.eligible.astype(bool)]
check('coverage_gate_closure', len(eligible_coverage) == 78 and len(skipped_coverage) == 4, f'eligible={len(eligible_coverage)} skipped={len(skipped_coverage)}')
check('coverage_skip_reason', (skipped_coverage.n_genes_present < 10).all(), str(skipped_coverage.n_genes_present.tolist()))
check('coverage_min_tested_genes', int(eligible_coverage.n_genes_present.min()) >= 10, str(int(eligible_coverage.n_genes_present.min())))

correlation = read_tsv('phase67_correlations_by_patient_timepoint_v1.tsv')
summary = read_tsv('phase67_external_replication_summary_v1.tsv')
cross = read_tsv('phase67_cross_dataset_replication_v1.tsv')
check('correlation_rows', len(correlation) == 1158, str(len(correlation)))
key_columns = ['dataset', 'endpoint', 'family', 'signature_id', 'patient_key', 'timepoint']
check('correlation_keys_unique', not correlation.duplicated(key_columns).any(), 'unique')
check('correlations_finite', np.isfinite(correlation[['raw_spearman_rho', 'cycle_adjusted_spearman_rho']].to_numpy(float)).all(), 'all finite')
check('summary_rows', len(summary) == 78, str(len(summary)))
check('cross_rows', len(cross) == 50 and cross.signature_id.nunique() == 50, str(len(cross)))

patient = correlation.groupby(['dataset', 'endpoint', 'family', 'signature_id', 'patient_key'], observed=True).agg(
    median_adjusted=('cycle_adjusted_spearman_rho', 'median')
).reset_index()
recalc_median = patient.groupby(['dataset', 'endpoint', 'family', 'signature_id'], observed=True).median_adjusted.median()
stored_median = summary.set_index(['dataset', 'endpoint', 'family', 'signature_id']).median_cycle_adjusted_rho
joined = pd.concat([recalc_median.rename('recalc'), stored_median.rename('stored')], axis=1)
check('patient_median_recalculation', np.allclose(joined.recalc, joined.stored, atol=1e-12), str(float(np.max(np.abs(joined.recalc - joined.stored)))))

bh_ok = True
for _, group in summary.groupby(['dataset', 'endpoint', 'family'], observed=True):
    values = group.adjusted_p_value.to_numpy(float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    output = np.empty(len(values)); output[order] = np.minimum(adjusted, 1.0)
    if not np.allclose(output, group.adjusted_q_value.to_numpy(float), atol=1e-12):
        bh_ok = False
check('bh_recalculation', bh_ok, 'all dataset-endpoint-family groups')
check('external_replication_logic', all(bool(row.external_replication) == (bool(row.direction_matches_phase64) and row.adjusted_q_value < 0.05) for row in summary.itertuples(index=False)), 'exact')

grade_counts = Counter(cross.replication_grade)
check('cross_grade_counts', sum(grade_counts.values()) == 50 and grade_counts['R2_TWO_EXTERNAL_FDR'] > 0 and grade_counts['R1_SINGLE_EXTERNAL_FDR'] > 0, str(dict(grade_counts)))
retained = read_tsv('phase67_cross_dataset_retained_v1.tsv')
check('retained_closure', len(retained) == sum(value for key, value in grade_counts.items() if key != 'R0_NOT_REPLICATED'), str(len(retained)))

summary_text = (RESULTS / 'PHASE67_DIRECT_REPLICATION_SUMMARY_v1.md').read_text(encoding='utf-8')
check('interpretation_boundary', '仍是机制复现，不等同于 response/PFS/OS 预测验证' in summary_text and '不能把这些结果称为双外部队列复现' in summary_text, 'present')

validation_path = RESULTS / 'PHASE67_VALIDATION_v1.tsv'
with validation_path.open('w', encoding='utf-8', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=['check', 'status', 'detail'], delimiter='\t', lineterminator='\n')
    writer.writeheader(); writer.writerows(checks)

manifest_rows = []
for path in sorted(RESULTS.iterdir()):
    if path.is_file() and path.name not in {'FILE_MANIFEST_SHA256_v1.tsv', 'PHASE67_PASS_v1'}:
        manifest_rows.append({'file': path.name, 'size_bytes': path.stat().st_size, 'sha256': sha256(path)})
with (RESULTS / 'FILE_MANIFEST_SHA256_v1.tsv').open('w', encoding='utf-8', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=['file', 'size_bytes', 'sha256'], delimiter='\t', lineterminator='\n')
    writer.writeheader(); writer.writerows(manifest_rows)

failures = [row for row in checks if row['status'] != 'PASS']
if failures:
    raise SystemExit('Phase67 validation failed: ' + ', '.join(row['check'] for row in failures))
(RESULTS / 'PHASE67_PASS_v1').write_text(f'status=PASS\nchecks={len(checks)}\n', encoding='utf-8')
print(json.dumps({'status': 'PASS', 'checks': len(checks), 'manifest_rows': len(manifest_rows)}, indent=2))
