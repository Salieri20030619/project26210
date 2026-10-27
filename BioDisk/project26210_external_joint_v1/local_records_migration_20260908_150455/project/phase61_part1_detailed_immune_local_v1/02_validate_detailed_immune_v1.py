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


def atomic_json(payload: dict, path: Path) -> None:
    part = path.with_name(path.name + '.part')
    part.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    os.replace(part, path)


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit('Usage: 02_validate_detailed_immune_v1.py <output_dir>')
    out = Path(sys.argv[1]).resolve()
    required = [
        'immune_features_long_v1.tsv',
        'GSE294109_icb_immune_response_comparison_v1.tsv',
        'GSE264586_immune_longitudinal_patient_deltas_v1.tsv',
        'GSE264586_immune_longitudinal_summary_v1.tsv',
        'GSE294109_multiregion_immune_heterogeneity_v1.tsv',
        'immune_feature_vs_f1_tymp_rho_associations_v1.tsv',
        'PHASE61_STATUS_v1.json',
    ]
    checks = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({'check': name, 'result': 'PASS' if passed else 'FAIL', 'detail': detail})

    missing = [name for name in required if not (out / name).is_file() or (out / name).stat().st_size == 0]
    check('required_outputs', not missing, 'missing=' + ';'.join(missing))
    if missing:
        raise RuntimeError(missing)

    features = pd.read_csv(out / 'immune_features_long_v1.tsv', sep='\t')
    key = ['dataset', 'unit', 'patient_key', 'timepoint', 'feature_family', 'feature_id']
    sample_key = key + ['sample_id']
    patient = features[features.unit == 'patient_timepoint']
    sample = features[features.unit == 'sample_region']
    check('patient_feature_keys_unique', not patient.duplicated(key).any(), f'rows={len(patient)}')
    check('sample_feature_keys_unique', not sample.duplicated(sample_key).any(), f'rows={len(sample)}')
    fractions = features[features.feature_family.isin(['immune_abundance', 'broad_composition', 'fine_composition'])]
    check('composition_values_in_unit_interval', fractions.value.dropna().between(0, 1).all(), f'rows={len(fractions)}')

    response = pd.read_csv(out / 'GSE294109_icb_immune_response_comparison_v1.tsv', sep='\t')
    anchor = response[response.feature_id == 'immune_fraction_of_all_cells']
    response_ok = len(anchor) == 1 and int(anchor.iloc[0].n_benefit) == 9 and int(anchor.iloc[0].n_no_effect) == 3
    check('response_patient_unit_9_vs_3', response_ok, anchor.to_json(orient='records'))
    check('response_fdr_valid', response.p_fdr_bh_within_family.dropna().between(0, 1).all(), f'tests={len(response)}')

    longitudinal = pd.read_csv(out / 'GSE264586_immune_longitudinal_patient_deltas_v1.tsv', sep='\t')
    anchor_long = longitudinal[longitudinal.feature_id == 'immune_fraction_of_all_cells']
    pair_counts = anchor_long.groupby('comparison').patient_id.nunique().to_dict()
    check(
        'longitudinal_pairs_conserved',
        pair_counts == {'Cycle2_vs_Baseline': 4, 'End_of_treatment_vs_Baseline': 3},
        str(pair_counts),
    )

    multiregion = pd.read_csv(out / 'GSE294109_multiregion_immune_heterogeneity_v1.tsv', sep='\t')
    patients = set(multiregion.patient_id.unique())
    check('multiregion_patients_exact', patients == {'sc02', 'sc03', 'sc06', 'sc08', 'sc10'}, str(sorted(patients)))
    check('multiregion_counts_ge2', (multiregion.n_regions >= 2).all(), f'rows={len(multiregion)}')

    associations = pd.read_csv(out / 'immune_feature_vs_f1_tymp_rho_associations_v1.tsv', sep='\t')
    check('association_datasets_exact', set(associations.dataset) == {'GSE264586', 'GSE294109'}, str(sorted(associations.dataset.unique())))
    check(
        'association_tiers_exact',
        set(associations.malignancy_tier) == {'strict_high_confidence', 'expanded_any_method_positive'},
        str(sorted(associations.malignancy_tier.unique())),
    )
    check(
        'association_fdr_valid',
        associations.p_fdr_bh_within_dataset_tier_family.dropna().between(0, 1).all(),
        f'tests={len(associations)}',
    )
    check(
        'association_minimum_n_for_inference',
        associations.loc[associations.n_patients < 5, 'p_value'].isna().all()
        and associations.loc[associations.p_value.notna(), 'n_patients'].ge(5).all(),
        f'tested={(associations.p_value.notna()).sum()}',
    )

    status = json.loads((out / 'PHASE61_STATUS_v1.json').read_text(encoding='utf-8'))
    check('analysis_status', status.get('status') == 'PASS_PHASE61_DETAILED_IMMUNE', str(status.get('status')))
    check('no_partial_files', not list(out.glob('*.part')), 'no .part files')

    validation = pd.DataFrame(checks)
    atomic_tsv(validation, out / 'PHASE61_VALIDATION_v1.tsv')
    passed = bool((validation.result == 'PASS').all())
    excluded = {'OUTPUT_MANIFEST_SHA256_v1.tsv', 'PHASE61_VALIDATION_STATUS_v1.json'}
    manifest_rows = []
    for path in sorted(out.glob('*')):
        if path.is_file() and path.name not in excluded and not path.name.endswith('.part'):
            manifest_rows.append({
                'file': path.name,
                'bytes': path.stat().st_size,
                'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            })
    atomic_tsv(pd.DataFrame(manifest_rows), out / 'OUTPUT_MANIFEST_SHA256_v1.tsv')
    payload = {
        'status': 'PASS_PHASE61_VALIDATION' if passed else 'FAIL_PHASE61_VALIDATION',
        'completed_at': datetime.now().astimezone().isoformat(timespec='seconds'),
        'checks': len(validation),
        'passed': int((validation.result == 'PASS').sum()),
        'failed': int((validation.result == 'FAIL').sum()),
        'output_files': len(manifest_rows),
    }
    atomic_json(payload, out / 'PHASE61_VALIDATION_STATUS_v1.json')
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    if not passed:
        raise SystemExit(2)


if __name__ == '__main__':
    main()
