#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / 'inputs_v1'
RESULTS = ROOT / 'results_v1'


def read_tsv(path: Path) -> list[dict[str, str]]:
    if path.suffix == '.gz':
        handle_context = gzip.open(path, mode='rt', encoding='utf-8-sig', newline='')
    else:
        handle_context = path.open(mode='r', encoding='utf-8-sig', newline='')
    with handle_context as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


checks: list[dict[str, str]] = []


def check(name: str, condition: bool, detail: str) -> None:
    checks.append({'check': name, 'status': 'PASS' if condition else 'FAIL', 'detail': detail})


phase62 = read_tsv(INPUT / 'PHASE62_STATUS_v1.tsv')
check('phase62_status', len(phase62) == 1 and phase62[0]['status'] == 'PASS_PHASE62_IMMUNE_PSEUDOBULK', phase62[0]['status'])
check('phase62_gene_rows', int(phase62[0]['gene_rows']) == 133463, phase62[0]['gene_rows'])

with (INPUT / 'extraction_metadata_v1.json').open(encoding='utf-8') as handle:
    extraction = json.load(handle)
check('regulon_extraction', extraction['status'] == 'PASS_FROZEN_REGULON_EXTRACTION_V1', extraction['status'])
check('regulon_tf_counts', extraction['myeloid_tfs'] == 15 and extraction['treg_tfs'] == 8, str(extraction))
targets = read_tsv(INPUT / 'phase64_frozen_regulon_targets_v1.tsv')
check('regulon_target_rows', len(targets) == extraction['target_rows'] == 15430, str(len(targets)))

registry = read_tsv(RESULTS / 'TRANSFER_SIGNATURE_REGISTRY_v1.tsv')
counts = Counter((row['endpoint'], row['family']) for row in registry)
check('signature_counts', counts == Counter({('myeloid_dc_tymp', 'regulon'): 15, ('myeloid_dc_tymp', 'metabolic'): 17, ('treg_tymp', 'regulon'): 8, ('treg_tymp', 'metabolic'): 10}), str(dict(counts)))
check('signature_ids_unique', len({row['signature_id'] for row in registry}) == len(registry), str(len(registry)))

coverage = read_tsv(RESULTS / 'TRANSFER_SIGNATURE_COVERAGE_v1.tsv')
full = read_tsv(RESULTS / 'TRANSFER_GSEA_FULL_v1.tsv')
summary = read_tsv(RESULTS / 'TRANSFER_MODEL_SUMMARY_v1.tsv')
check('five_models', len({row['model_id'] for row in coverage}) == 5, str(len({row['model_id'] for row in coverage})))
check('coverage_rows', len(coverage) == 146, str(len(coverage)))
eligible = [row for row in coverage if row['eligible'].lower() == 'true']
skipped = [row for row in coverage if row['eligible'].lower() != 'true']
check('coverage_gate_closure', len(eligible) == 123 and len(skipped) == 23, f'eligible={len(eligible)} skipped={len(skipped)}')
check('coverage_skip_reason', all(int(row['n_genes_in_model']) < 10 or int(row['n_genes_in_model']) > 500 for row in skipped), 'all skipped outside 10-500 genes')
check('gsea_row_closure', len(full) == len(eligible), f'{len(full)}/{len(eligible)}')
check('summary_rows', len(summary) == 10, str(len(summary)))

check('finite_p_q', all(0 <= float(row['pval']) <= 1 and 0 <= float(row['transfer_q_value']) <= 1 for row in full), 'all finite in [0,1]')
groups = {}
for row in full:
    groups.setdefault((row['model_id'], row['family']), []).append(row)
bh_ok = True
for rows in groups.values():
    ordered = sorted(rows, key=lambda row: float(row['pval']))
    adjusted = [0.0] * len(ordered)
    running = 1.0
    for index in range(len(ordered) - 1, -1, -1):
        running = min(running, float(ordered[index]['pval']) * len(ordered) / (index + 1))
        adjusted[index] = min(running, 1.0)
    if any(abs(value - float(row['transfer_q_value'])) > 1e-12 for value, row in zip(adjusted, ordered)):
        bh_ok = False
check('bh_recalculation', bh_ok, f'groups={len(groups)}')
check('grade_consistency', all(
    (row['transfer_grade'] == 'T1_FDR_CONTEXT_TRANSFER') == (float(row['transfer_q_value']) < 0.05)
    for row in full
), 'T1 iff q<0.05')
check('status_counts', sum(float(row['transfer_q_value']) < 0.05 for row in full) == 42 and sum(float(row['pval']) < 0.05 <= float(row['transfer_q_value']) for row in full) == 10, 'FDR=42 nominal_only=10')
check('no_tymp_in_metabolic_leading_edge', all(
    row['family'] != 'metabolic' or 'TYMP' not in row.get('leadingEdge', '').split(';')
    for row in full
), 'TYMP excluded')

summary_text = (RESULTS / 'PHASE66_TRANSFER_SUMMARY_v1.md').read_text(encoding='utf-8')
check('interpretation_boundary', '不是 TYMP 连续表达相关性的直接复算' in summary_text and '不能单独把程序升级为疗效预测标志物' in summary_text, 'boundary present')
interpretation_text = (RESULTS / 'PHASE66_INTERPRETATION_v1.md').read_text(encoding='utf-8')
check('interpretation_result_counts', 'regulon 35 个、代谢通路 7 个' in interpretation_text and '另有 10 个仅名义显著结果' in interpretation_text, '42 FDR plus 10 nominal-only documented')
check('interpretation_context_boundary', '不存在统一、方向一致的“疗效获益程序”' in interpretation_text and '直接 response/PFS/OS 预测等级保持不变' in interpretation_text, 'context dependence preserved')

validation_path = RESULTS / 'PHASE66_VALIDATION_v1.tsv'
with validation_path.open('w', encoding='utf-8', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=['check', 'status', 'detail'], delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(checks)

failures = [row for row in checks if row['status'] != 'PASS']
manifest_rows = []
for path in sorted(INPUT.iterdir()):
    if path.is_file():
        manifest_rows.append({'scope': 'input', 'file': path.name, 'size_bytes': path.stat().st_size, 'sha256': sha256(path)})
for path in sorted(RESULTS.iterdir()):
    if path.is_file() and path.name not in {'FILE_MANIFEST_SHA256_v1.tsv', 'PHASE66_PASS_v1'}:
        manifest_rows.append({'scope': 'output', 'file': path.name, 'size_bytes': path.stat().st_size, 'sha256': sha256(path)})
with (RESULTS / 'FILE_MANIFEST_SHA256_v1.tsv').open('w', encoding='utf-8', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=['scope', 'file', 'size_bytes', 'sha256'], delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(manifest_rows)

if failures:
    raise SystemExit('Phase66 validation failed: ' + ', '.join(row['check'] for row in failures))
(RESULTS / 'PHASE66_PASS_v1').write_text(f'status=PASS\nchecks={len(checks)}\n', encoding='utf-8')
print(json.dumps({'status': 'PASS', 'checks': len(checks), 'manifest_rows': len(manifest_rows)}, indent=2))
