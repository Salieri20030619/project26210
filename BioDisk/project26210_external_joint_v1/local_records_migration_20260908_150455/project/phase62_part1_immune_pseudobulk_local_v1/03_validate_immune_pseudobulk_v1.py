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
        raise SystemExit('Usage: 03_validate_immune_pseudobulk_v1.py <output_dir>')
    out = Path(sys.argv[1]).resolve()
    required = [
        'pseudobulk_input_audit_v1.tsv',
        'pseudobulk_unit_audit_v1.tsv',
        'pseudobulk_model_summary_v1.tsv',
        'pseudobulk_gene_differential_expression_v1.tsv.gz',
        'pseudobulk_hallmark_fgsea_v1.tsv',
        'PHASE62_STATUS_v1.tsv',
    ]
    checks = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({'check': name, 'result': 'PASS' if passed else 'FAIL', 'detail': detail})

    missing = [name for name in required if not (out / name).is_file() or (out / name).stat().st_size == 0]
    check('required_outputs', not missing, 'missing=' + ';'.join(missing))
    if missing:
        raise RuntimeError(missing)

    inputs = pd.read_csv(out / 'pseudobulk_input_audit_v1.tsv', sep='\t')
    check('datasets_exact', set(inputs.dataset) == {'GSE264586', 'GSE294109'}, str(sorted(inputs.dataset)))
    mapping_ok = (inputs.object_cells >= inputs.mapped_cells).all() and (inputs.mapped_cells == inputs.compact_cells).all()
    check('cell_mapping_complete', bool(mapping_ok), inputs.to_json(orient='records'))

    models = pd.read_csv(out / 'pseudobulk_model_summary_v1.tsv', sep='\t')
    check('all_models_tested', len(models) == 12 and models.status.eq('TESTED').all(), f'models={len(models)}')
    response = models[models.dataset == 'GSE294109']
    check(
        'response_patient_counts_bounded',
        response.n_benefit.between(6, 9).all() and response.n_no_effect.eq(3).all()
        and (response.n_units == response.n_benefit + response.n_no_effect).all(),
        response[['fine_annotation', 'n_benefit', 'n_no_effect']].to_json(orient='records'),
    )
    cd8 = response[response.fine_annotation == 'CD8_cytotoxic']
    check('cd8_full_9_vs_3', len(cd8) == 1 and int(cd8.iloc[0].n_benefit) == 9 and int(cd8.iloc[0].n_no_effect) == 3, cd8.to_json(orient='records'))
    longitudinal = models[models.dataset == 'GSE264586']
    check('paired_design_units', (longitudinal.n_units == 2 * longitudinal.n_pairs).all(), longitudinal.to_json(orient='records'))
    check(
        'paired_counts_expected',
        set(longitudinal.loc[longitudinal.contrast == 'Cycle2_vs_Baseline', 'n_pairs']) == {4}
        and set(longitudinal.loc[longitudinal.contrast == 'End_of_treatment_vs_Baseline', 'n_pairs']) == {3},
        longitudinal[['contrast', 'fine_annotation', 'n_pairs']].to_json(orient='records'),
    )

    units = pd.read_csv(out / 'pseudobulk_unit_audit_v1.tsv', sep='\t')
    unit_key = ['dataset', 'contrast', 'fine_annotation', 'patient_key', 'timepoint']
    check('pseudobulk_units_unique', not units.duplicated(unit_key).any(), f'units={len(units)}')
    check('minimum_50_cells', units.n_cells.ge(50).all(), f'min={int(units.n_cells.min())}')

    genes = pd.read_csv(out / 'pseudobulk_gene_differential_expression_v1.tsv.gz', sep='\t')
    gene_key = ['dataset', 'contrast', 'fine_annotation', 'gene']
    check('gene_keys_unique', not genes.duplicated(gene_key).any(), f'rows={len(genes)}')
    check('gene_fdr_valid', genes.FDR.dropna().between(0, 1).all(), f'rows={len(genes)}')
    modeled_keys = set(map(tuple, models[['dataset', 'contrast', 'fine_annotation']].to_numpy()))
    gene_model_keys = set(map(tuple, genes[['dataset', 'contrast', 'fine_annotation']].drop_duplicates().to_numpy()))
    check('gene_models_complete', gene_model_keys == modeled_keys, f'gene_models={len(gene_model_keys)}')

    pathways = pd.read_csv(out / 'pseudobulk_hallmark_fgsea_v1.tsv', sep='\t')
    check('pathway_padj_valid', pathways.padj.dropna().between(0, 1).all(), f'rows={len(pathways)}')
    pathway_model_keys = set(map(tuple, pathways[['dataset', 'contrast', 'fine_annotation']].drop_duplicates().to_numpy()))
    check('pathway_models_complete', pathway_model_keys == modeled_keys, f'pathway_models={len(pathway_model_keys)}')

    status = pd.read_csv(out / 'PHASE62_STATUS_v1.tsv', sep='\t').iloc[0]
    check(
        'analysis_status',
        status.status == 'PASS_PHASE62_IMMUNE_PSEUDOBULK' and int(status.models_tested) == 12,
        status.to_json(),
    )
    check('no_partial_files', not list(out.glob('*.part')) and not list(out.glob('*.part.gz')), 'no partial files')

    validation = pd.DataFrame(checks)
    atomic_tsv(validation, out / 'PHASE62_VALIDATION_v1.tsv')
    passed = bool(validation.result.eq('PASS').all())
    excluded = {'OUTPUT_MANIFEST_SHA256_v1.tsv', 'PHASE62_VALIDATION_STATUS_v1.json'}
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
        'status': 'PASS_PHASE62_VALIDATION' if passed else 'FAIL_PHASE62_VALIDATION',
        'completed_at': datetime.now().astimezone().isoformat(timespec='seconds'),
        'checks': len(validation),
        'passed': int(validation.result.eq('PASS').sum()),
        'failed': int(validation.result.eq('FAIL').sum()),
        'output_files': len(manifest_rows),
    }
    atomic_json(payload, out / 'PHASE62_VALIDATION_STATUS_v1.json')
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    if not passed:
        raise SystemExit(2)


if __name__ == '__main__':
    main()
