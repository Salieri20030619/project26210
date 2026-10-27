#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path

ROOT = Path('/mnt/BioDisk/wangmingyuan/project26210_phase73_virtual_knockout_v1')
RESULTS = ROOT / 'results_v1'
VALID = ROOT / 'independent_validation_v1'
DATASETS = ('GSE207493', 'GSE264586', 'GSE294109')
TARGETS = ('STAT3', 'SPP1', 'TYMP', 'MTAP')
PROGRAMS = ('MTA_A2B_STAT3_RESPONSE','SPP1_STATE_NO_SPP1','IMMUNOSUPPRESSION_NO_IL10','MTA_PRODUCTION_CONTEXT','C1QC_STATE','MONOCYTE_PROGRAM')


def rows(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


checks: list[dict[str, str]] = []


def check(name: str, passed: bool, detail: str) -> None:
    checks.append({'check': name, 'status': 'PASS' if passed else 'FAIL', 'detail': detail})


def main() -> None:
    VALID.mkdir(parents=True, exist_ok=True)
    preflight = rows(ROOT / 'preflight_v1' / 'input_preflight_v1.tsv')
    check('preflight_rows', len(preflight) == 12, f'rows={len(preflight)}')
    check('preflight_targets_present', all(r['target_present'] == 'TRUE' for r in preflight), 'all targets present')
    check('preflight_counts_layer', all(r['counts_layer'] == 'counts' for r in preflight), 'raw counts used')

    for dataset in DATASETS:
        ddir = ROOT / 'dataset_runs_v1' / dataset
        check(f'{dataset}_complete_marker', (ddir / 'DATASET_COMPLETE').exists(), str(ddir / 'DATASET_COMPLETE'))
        status_rows = rows(ddir / f'{dataset}_status_v1.tsv')
        check(f'{dataset}_status', len(status_rows) == 1 and status_rows[0]['status'] == 'PASS_PHASE73_DATASET_V1', str(status_rows))
        program_rows = rows(ddir / f'{dataset}_program_impact_all_v1.tsv')
        check(f'{dataset}_program_rows', len(program_rows) == 24, f'rows={len(program_rows)}')
        keys = {(r['knockout_gene'], r['program']) for r in program_rows}
        check(f'{dataset}_program_keys', keys == {(t,p) for t in TARGETS for p in PROGRAMS}, f'keys={len(keys)}')
        numeric_ok = True
        for r in program_rows:
            for field in ('n_modeled_genes','median_impact_percentile','p_value','q_value'):
                try:
                    value = float(r[field])
                    numeric_ok &= math.isfinite(value)
                except Exception:
                    numeric_ok = False
        check(f'{dataset}_program_numeric', numeric_ok, 'finite program statistics')
        for target in TARGETS:
            diff_path = ddir / f'KO_{target}' / f'{dataset}_{target}_diff_regulation_v1.tsv.gz'
            diff = rows(diff_path)
            genes = [r['affected_gene'] for r in diff]
            check(f'{dataset}_{target}_diff_rows', len(diff) >= 2500, f'rows={len(diff)}')
            check(f'{dataset}_{target}_unique_genes', len(genes) == len(set(genes)), f'unique={len(set(genes))}')
            check(f'{dataset}_{target}_target_present', target in set(genes), f'target={target}')

    cross = rows(RESULTS / 'phase73_cross_cohort_decisions_v1.tsv')
    check('cross_rows', len(cross) == 24, f'rows={len(cross)}')
    check('cross_unique_keys', len({(r['knockout_gene'],r['program']) for r in cross}) == 24, 'unique KO-program pairs')
    check('primary_rows', sum(r['hypothesis_tier'] == 'PRIMARY' for r in cross) == 10, 'expected 10 primary hypotheses')
    allowed = {'TWO_EXTERNAL_COMPUTATIONAL_KO_SUPPORTED','ONE_EXTERNAL_COMPUTATIONAL_KO_SUPPORTED','NOT_EXTERNALLY_SUPPORTED'}
    check('cross_decisions', all(r['external_decision'] in allowed for r in cross), 'allowed decisions only')
    check('direction_boundary', all(r['direction_available'] == 'FALSE' for r in cross), 'direction unavailable recorded for every pair')

    manifest = rows(RESULTS / 'FILE_MANIFEST_SHA256_v1.tsv')
    hash_ok = True
    for r in manifest:
        path = ROOT / r['path']
        hash_ok &= path.exists() and path.stat().st_size == int(r['bytes']) and sha256(path) == r['sha256']
    check('manifest_hashes', hash_ok and len(manifest) >= 3, f'entries={len(manifest)}')
    metadata = json.loads((RESULTS / 'run_metadata_v1.json').read_text(encoding='utf-8'))
    check('metadata_status', metadata.get('status') == 'PASS_PHASE73_FINALIZATION_VALIDATION_PENDING', str(metadata.get('status')))
    check('metadata_dimensions', metadata.get('program_rows') == 72 and metadata.get('cross_rows') == 24, str(metadata))

    failures = [x for x in checks if x['status'] != 'PASS']
    checks_path = VALID / 'validation_checks_v1.tsv'
    tmp = checks_path.with_suffix('.tsv.part')
    with tmp.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['check','status','detail'], delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(checks)
    os.replace(tmp, checks_path)
    summary = {'status': 'PASS' if not failures else 'FAIL', 'checks': len(checks), 'passed': len(checks)-len(failures), 'failed': len(failures)}
    summary_path = VALID / 'validation_summary_v1.json'
    summary_path.write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    if failures:
        raise SystemExit('PHASE73_VALIDATION_FAIL ' + json.dumps(failures))
    (VALID / 'INDEPENDENT_VALIDATION_PASS_v1').write_text('PHASE73_INDEPENDENT_VALIDATION_PASS\n', encoding='utf-8')
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
