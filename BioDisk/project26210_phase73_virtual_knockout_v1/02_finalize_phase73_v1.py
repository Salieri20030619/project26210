#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path('/mnt/BioDisk/wangmingyuan/project26210_phase73_virtual_knockout_v1')
RESULTS = ROOT / 'results_v1'
DATASETS = ('GSE207493', 'GSE264586', 'GSE294109')
TARGETS = ('STAT3', 'SPP1', 'TYMP', 'MTAP')
PROGRAMS = (
    'MTA_A2B_STAT3_RESPONSE', 'SPP1_STATE_NO_SPP1',
    'IMMUNOSUPPRESSION_NO_IL10', 'MTA_PRODUCTION_CONTEXT',
    'C1QC_STATE', 'MONOCYTE_PROGRAM')


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open('r', encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    tmp = path.with_suffix(path.suffix + '.part')
    with tmp.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict[str, str]] = []
    for dataset in DATASETS:
        dataset_dir = ROOT / 'dataset_runs_v1' / dataset
        marker = dataset_dir / 'DATASET_COMPLETE'
        source = dataset_dir / f'{dataset}_program_impact_all_v1.tsv'
        if not marker.exists() or not source.exists():
            raise RuntimeError(f'incomplete dataset: {dataset}')
        rows = read_tsv(source)
        if len(rows) != len(TARGETS) * len(PROGRAMS):
            raise RuntimeError(f'unexpected program rows for {dataset}: {len(rows)}')
        all_rows.extend(rows)

    all_path = RESULTS / 'phase73_program_impact_all_v1.tsv'
    fields = list(all_rows[0])
    write_tsv(all_path, all_rows, fields)
    index = {(r['dataset'], r['knockout_gene'], r['program']): r for r in all_rows}

    cross_rows: list[dict[str, object]] = []
    for target in TARGETS:
        for program in PROGRAMS:
            d = {dataset: index[(dataset, target, program)] for dataset in DATASETS}
            ext_supported = [x for x in ('GSE264586', 'GSE294109')
                             if d[x]['decision'] == 'PROGRAM_IMPACT_SUPPORTED']
            if len(ext_supported) == 2:
                decision = 'TWO_EXTERNAL_COMPUTATIONAL_KO_SUPPORTED'
            elif len(ext_supported) == 1:
                decision = 'ONE_EXTERNAL_COMPUTATIONAL_KO_SUPPORTED'
            else:
                decision = 'NOT_EXTERNALLY_SUPPORTED'
            cross_rows.append({
                'knockout_gene': target,
                'program': program,
                'hypothesis_tier': d['GSE207493']['hypothesis_tier'],
                'external_decision': decision,
                'GSE207493_decision': d['GSE207493']['decision'],
                'GSE207493_median_percentile': d['GSE207493']['median_impact_percentile'],
                'GSE207493_q_value': d['GSE207493']['q_value'],
                'GSE264586_decision': d['GSE264586']['decision'],
                'GSE264586_median_percentile': d['GSE264586']['median_impact_percentile'],
                'GSE264586_q_value': d['GSE264586']['q_value'],
                'GSE294109_decision': d['GSE294109']['decision'],
                'GSE294109_median_percentile': d['GSE294109']['median_impact_percentile'],
                'GSE294109_q_value': d['GSE294109']['q_value'],
                'external_supported_cohorts': ';'.join(ext_supported),
                'direction_available': 'FALSE',
                'boundary': 'Unsigned observational-GRN perturbation impact; not measured KO, expression direction, MTA abundance, protein activity, or clinical causality.'
            })
    cross_fields = list(cross_rows[0])
    cross_path = RESULTS / 'phase73_cross_cohort_decisions_v1.tsv'
    write_tsv(cross_path, cross_rows, cross_fields)

    primary = [r for r in cross_rows if r['hypothesis_tier'] == 'PRIMARY']
    two = [r for r in primary if r['external_decision'] == 'TWO_EXTERNAL_COMPUTATIONAL_KO_SUPPORTED']
    one = [r for r in primary if r['external_decision'] == 'ONE_EXTERNAL_COMPUTATIONAL_KO_SUPPORTED']
    summary = [
        '# Phase73 virtual knockout summary v1', '',
        '## Scope', '',
        'Four genes (STAT3, SPP1, TYMP and MTAP) were virtually deleted from independently reconstructed cohort-level myeloid GRNs using the frozen Phase71 cells and raw RNA counts.', '',
        '## Decision counts', '',
        f'- Primary KO-program hypotheses: {len(primary)}.',
        f'- Two-external computational KO supported: {len(two)}.',
        f'- One-external computational KO supported: {len(one)}.',
        f'- Not externally supported: {len(primary)-len(two)-len(one)}.', '',
        '## Supported primary pairs', ''
    ]
    if two:
        summary.extend([f"- TWO_EXTERNAL: {r['knockout_gene']} KO -> {r['program']} network impact." for r in two])
    else:
        summary.append('- No primary pair reached the two-external rule.')
    if one:
        summary.extend([f"- ONE_EXTERNAL: {r['knockout_gene']} KO -> {r['program']} network impact ({r['external_supported_cohorts']})." for r in one])
    summary.extend([
        '', '## Interpretation boundary', '',
        'scTenifoldKnk differential regulation is an unsigned network-embedding distance. These results cannot be described as program increase or decrease. They do not measure 5′-MTA, STAT3 phosphorylation, SPP1 protein secretion, cell-cell functional response, tumor progression or clinical outcome. Phase69 clinical grade D is unchanged.', ''
    ])
    summary_path = RESULTS / 'PHASE73_VIRTUAL_KNOCKOUT_SUMMARY_v1.md'
    tmp_summary = summary_path.with_suffix('.md.part')
    tmp_summary.write_text('\n'.join(summary), encoding='utf-8')
    os.replace(tmp_summary, summary_path)

    manifest_rows = []
    for path in sorted(RESULTS.rglob('*')):
        if path.is_file() and path.name != 'FILE_MANIFEST_SHA256_v1.tsv':
            manifest_rows.append({'path': str(path.relative_to(ROOT)), 'bytes': path.stat().st_size, 'sha256': sha256(path)})
    write_tsv(RESULTS / 'FILE_MANIFEST_SHA256_v1.tsv', manifest_rows, ['path', 'bytes', 'sha256'])
    metadata = {
        'status': 'PASS_PHASE73_FINALIZATION_VALIDATION_PENDING',
        'completed_at': datetime.now(timezone.utc).astimezone().isoformat(),
        'datasets': list(DATASETS), 'targets': list(TARGETS),
        'program_rows': len(all_rows), 'cross_rows': len(cross_rows),
        'primary_rows': len(primary), 'primary_two_external': len(two),
        'primary_one_external': len(one), 'direction_available': False
    }
    metadata_path = RESULTS / 'run_metadata_v1.json'
    tmp_metadata = metadata_path.with_suffix('.json.part')
    tmp_metadata.write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    os.replace(tmp_metadata, metadata_path)
    print(json.dumps(metadata, ensure_ascii=False))


if __name__ == '__main__':
    main()
