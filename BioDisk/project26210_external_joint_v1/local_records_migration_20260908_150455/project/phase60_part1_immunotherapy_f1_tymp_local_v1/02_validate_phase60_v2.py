#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import os
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    part = path.with_name(path.name + '.part')
    frame.to_csv(part, sep='\t', index=False, lineterminator='\n')
    os.replace(part, path)


def atomic_text(text: str, path: Path) -> None:
    part = path.with_name(path.name + '.part')
    part.write_text(text, encoding='utf-8')
    os.replace(part, path)


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit('Usage: 02_validate_phase60_v2.py <project_root> <output_dir>')
    project = Path(sys.argv[1]).resolve()
    out = Path(sys.argv[2]).resolve()
    checks: list[dict[str, str]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({'check': name, 'result': 'PASS' if passed else 'FAIL', 'detail': detail})

    required = [
        'clinical_mapping_audit_v1.tsv',
        'score_method_definitions_v1.tsv',
        'patient_f1_tymp_correlations_v1.tsv',
        'malignancy_tier_cell_counts_v1.tsv',
        'immune_composition_by_unit_v1.tsv',
        'immune_state_scores_patient_timepoint_v1.tsv',
        'GSE264586_longitudinal_delta_rho_v1.tsv',
        'GSE264586_longitudinal_summary_v1.tsv',
        'GSE294109_icb_patient_rho_v1.tsv',
        'GSE294109_icb_response_comparison_v1.tsv',
        'malignant_rho_immune_state_associations_v1.tsv',
        'malignant_rho_immune_composition_associations_v1.tsv',
        'PHASE60_STATUS_v1.tsv',
    ]
    missing = [name for name in required if not (out / name).is_file() or (out / name).stat().st_size == 0]
    check('required_outputs', not missing, 'missing=' + ';'.join(missing))
    if missing:
        raise RuntimeError(f'Missing required outputs: {missing}')

    corr = pd.read_csv(out / 'patient_f1_tymp_correlations_v1.tsv', sep='\t')
    key = ['dataset', 'unit', 'malignancy_tier', 'method', 'patient_key', 'sample_id', 'timepoint']
    check('correlation_keys_unique', not corr.duplicated(key).any(), f'rows={len(corr)}')
    check('datasets_exact', set(corr.dataset) == {'GSE264586', 'GSE294109'}, str(sorted(corr.dataset.unique())))
    check(
        'tiers_exact',
        set(corr.malignancy_tier) == {'strict_high_confidence', 'expanded_any_method_positive'},
        str(sorted(corr.malignancy_tier.unique())),
    )
    eligible = corr[(corr.unit == 'patient_timepoint') & (corr.status == 'ELIGIBLE_GE20')]
    check('eligible_rho_finite', eligible.spearman_rho.notna().all(), f'rows={len(eligible)}')

    count = pd.read_csv(out / 'malignancy_tier_cell_counts_v1.tsv', sep='\t')
    wide = count.pivot_table(
        index=['dataset', 'patient_key', 'sample_id', 'timepoint'],
        columns='malignancy_tier', values='n_cells', aggfunc='sum', fill_value=0,
    )
    strict = wide.get('strict_high_confidence', pd.Series(0, index=wide.index))
    expanded = wide.get('expanded_any_method_positive', pd.Series(0, index=wide.index))
    check('strict_nested_in_expanded', bool((strict <= expanded).all()), f'units={len(wide)}')

    mapping = pd.read_csv(out / 'clinical_mapping_audit_v1.tsv', sep='\t')
    expected = {'GSE264586': (19, 12), 'GSE294109': (25, 15)}
    mapping_ok = all(
        int(row.objects) == expected[row.dataset][0]
        and int(row.patients) == expected[row.dataset][1]
        and bool(row.mapping_pass)
        for row in mapping.itertuples()
    )
    check('clinical_mapping_conserved', mapping_ok, mapping.to_json(orient='records'))

    crosswalk = pd.read_csv(
        project / 'phase6_metadata_confirmation_local_v1' / 'tables' / 'confirmed_object_patient_crosswalk.tsv',
        sep='\t',
    )
    roster = crosswalk[
        (crosswalk.dataset == 'GSE294109')
        & (crosswalk.timepoint_confirmed == 'Surgery_after_treatment')
        & (crosswalk.treatment_confirmed != 'Cabozantinib')
        & crosswalk.response_confirmed.isin(['CR', 'PR', 'No_effect'])
    ].copy()
    roster['response_group'] = roster.response_confirmed.map({'CR': 'benefit', 'PR': 'benefit', 'No_effect': 'no_effect'})
    roster = roster.groupby(
        ['dataset', 'patient_key', 'patient_id_confirmed', 'treatment_confirmed', 'response_confirmed', 'response_group'],
        as_index=False,
    ).agg(n_region_samples=('sample_id_confirmed', 'nunique'))
    roster = roster.rename(columns={
        'patient_id_confirmed': 'patient_id',
        'treatment_confirmed': 'treatment',
        'response_confirmed': 'response',
    })
    atomic_tsv(roster, out / 'GSE294109_icb_response_roster_v1.tsv')
    group_counts = roster.groupby('response_group').patient_id.nunique().to_dict()
    check('gse294109_response_mapping', group_counts == {'benefit': 9, 'no_effect': 3}, str(group_counts))

    response = pd.read_csv(out / 'GSE294109_icb_patient_rho_v1.tsv', sep='\t')
    primary_strict = response[
        (response.malignancy_tier == 'strict_high_confidence') & (response.method == 'primary_gene_z')
    ]
    eligible_counts = primary_strict.groupby('response_group').patient_id.nunique().to_dict()
    check('strict_response_eligibility_recorded', eligible_counts == {'benefit': 4, 'no_effect': 1}, str(eligible_counts))
    check('gse264586_response_not_imputed', 'GSE264586' not in response.dataset.unique(), 'response table only GSE294109')

    status = pd.read_csv(out / 'PHASE60_STATUS_v1.tsv', sep='\t').iloc[0]
    check('analysis_status', status.status == 'PASS_PHASE60_F1_TYMP', status.status)
    check('no_partial_files', not list(out.glob('*.part')), 'no .part files')

    validation = pd.DataFrame(checks)
    atomic_tsv(validation, out / 'PHASE60_VALIDATION_v2.tsv')
    passed = bool((validation.result == 'PASS').all())

    manifest_rows = []
    excluded = {'OUTPUT_MANIFEST_SHA256_v2.tsv', 'PHASE60_VALIDATION_STATUS_v2.json'}
    for path in sorted(out.glob('*')):
        if path.is_file() and not path.name.endswith('.part') and path.name not in excluded:
            manifest_rows.append({
                'file': path.name,
                'bytes': path.stat().st_size,
                'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            })
    manifest = pd.DataFrame(manifest_rows)
    atomic_tsv(manifest, out / 'OUTPUT_MANIFEST_SHA256_v2.tsv')
    payload = {
        'status': 'PASS_PHASE60_VALIDATION' if passed else 'FAIL_PHASE60_VALIDATION',
        'completed_at': datetime.now().astimezone().isoformat(timespec='seconds'),
        'checks': len(validation),
        'passed': int((validation.result == 'PASS').sum()),
        'failed': int((validation.result == 'FAIL').sum()),
        'output_files': len(manifest),
    }
    atomic_text(json.dumps(payload, ensure_ascii=False, indent=2) + '\n', out / 'PHASE60_VALIDATION_STATUS_v2.json')
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    if not passed:
        raise SystemExit(2)


if __name__ == '__main__':
    main()
