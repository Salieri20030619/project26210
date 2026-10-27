#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import subprocess
import traceback
import zipfile
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats


SEED = 26210
MIN_CELLS = 50
MEMORY_GATE = Path('/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh')
DATASETS = ('GSE264586', 'GSE294109')
ENDPOINTS = ('myeloid_dc_tymp', 'treg_tymp')
LINEAGE = {'myeloid_dc_tymp': 'Myeloid_DC', 'treg_tymp': 'Treg'}


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec='seconds')


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    temporary = path.with_name(path.name + '.part')
    frame.to_csv(temporary, sep='\t', index=False, na_rep='NA', lineterminator='\n')
    os.replace(temporary, path)


def memory_gate(log_path: Path, stage: str) -> None:
    result = subprocess.run([str(MEMORY_GATE)], text=True, capture_output=True)
    line = (result.stdout + result.stderr).strip()
    with log_path.open('a', encoding='utf-8') as handle:
        handle.write(f'{now()}\t{stage}\treturncode={result.returncode}\t{line}\n')
    if result.returncode != 0:
        raise RuntimeError(f'Memory gate blocked {stage}: {line}')


def bh_adjust(values: pd.Series) -> np.ndarray:
    array = values.to_numpy(float)
    order = np.argsort(array)
    ranked = array[order]
    adjusted = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    output = np.empty(len(array), dtype=float)
    output[order] = np.minimum(adjusted, 1.0)
    return output


def signed_rank(values: pd.Series) -> tuple[float, float]:
    array = values.to_numpy(float)
    array = array[np.isfinite(array)]
    if not len(array) or np.allclose(array, 0):
        return 0.0, 1.0
    result = stats.wilcoxon(array, alternative='two-sided', zero_method='wilcox', method='auto')
    return float(result.statistic), float(result.pvalue)


def residualize(values: np.ndarray, s_score: np.ndarray, g2m_score: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(len(values)), s_score, g2m_score])
    coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
    return values - design @ coefficients


def parse_reactome(path: Path) -> dict[str, list[str]]:
    with zipfile.ZipFile(path) as archive:
        lines = archive.read('ReactomePathways.gmt').decode('utf-8').splitlines()
    result = {}
    for line in lines:
        fields = line.split('\t')
        result[fields[1]] = list(dict.fromkeys(gene.upper() for gene in fields[2:]))
    return result


def load_signatures(input_dir: Path) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    targets = pd.read_csv(input_dir / 'phase64_frozen_regulon_targets_v1.tsv', sep='\t')
    directions = pd.read_csv(input_dir / 'de_novo_tf_cross_method_FDR_v1.tsv', sep='\t')
    directions = directions.loc[directions.canonical_metric == 'continuous'].copy()
    directions['tf'] = directions['TF']
    directions['expected_direction'] = np.where(
        (directions.effect_SCENIC > 0) & (directions.effect_VIPER > 0), 'positive', 'negative'
    )
    direction_map = directions.set_index(['endpoint', 'tf'])['expected_direction'].to_dict()

    registry_rows = []
    genes_by_signature = {}
    targets = targets.sort_values(['endpoint', 'tf_order', 'tf', 'weight', 'target_gene'], ascending=[True, True, True, False, True])
    for (endpoint, tf_order, tf), group in targets.groupby(['endpoint', 'tf_order', 'tf'], sort=False):
        genes = list(dict.fromkeys(group.target_gene.astype(str).str.upper()))[:200]
        genes = [gene for gene in genes if gene != 'TYMP']
        signature_id = f'{endpoint}__REGULON__{int(tf_order):02d}__{tf}'
        genes_by_signature[signature_id] = genes
        registry_rows.append({
            'endpoint': endpoint,
            'family': 'regulon',
            'signature_id': signature_id,
            'label': f'{tf} positive regulon top-weighted targets',
            'expected_phase64_direction': direction_map[(endpoint, tf)],
            'n_frozen_genes': len(genes),
            'source': 'Phase64 de novo SCENIC frozen cross-method TF list',
        })

    reactome = parse_reactome(input_dir / 'ReactomePathways.gmt.zip')
    metabolic = pd.read_csv(input_dir / 'gse207493_metabolic_contrast_summary_retained_v1.tsv', sep='\t')
    metabolic = metabolic.loc[metabolic.adjusted_q_value < 0.05].copy()
    endpoint_map = {
        'TYMP_high_vs_low_MyeloidDC': 'myeloid_dc_tymp',
        'TYMP_high_vs_low_Treg': 'treg_tymp',
    }
    metabolic = metabolic.loc[metabolic.contrast.isin(endpoint_map)].copy()
    for row in metabolic.itertuples(index=False):
        endpoint = endpoint_map[row.contrast]
        genes = [gene for gene in reactome[row.stable_id] if gene != 'TYMP']
        signature_id = f'{endpoint}__REACTOME__{row.stable_id}'
        genes_by_signature[signature_id] = genes
        registry_rows.append({
            'endpoint': endpoint,
            'family': 'metabolic',
            'signature_id': signature_id,
            'label': row.pathway,
            'expected_phase64_direction': 'positive' if row.direction == 'higher_in_feature_high' else 'negative',
            'n_frozen_genes': len(genes),
            'source': 'Phase64 cycle-adjusted TYMP-excluded Reactome contrast FDR<0.05',
        })
    registry = pd.DataFrame(registry_rows).sort_values(['endpoint', 'family', 'signature_id']).reset_index(drop=True)
    if len(registry) != 50 or registry.signature_id.duplicated().any():
        raise RuntimeError(f'Unexpected signature registry: rows={len(registry)}')
    return registry, genes_by_signature


def build_summary(correlations: pd.DataFrame, registry: pd.DataFrame) -> pd.DataFrame:
    patient = correlations.groupby(
        ['dataset', 'endpoint', 'family', 'signature_id', 'patient_key'], observed=True
    ).agg(
        patient_median_raw_rho=('raw_spearman_rho', 'median'),
        patient_median_adjusted_rho=('cycle_adjusted_spearman_rho', 'median'),
        n_timepoints=('timepoint', 'nunique'),
        total_cells=('n_cells', 'sum'),
    ).reset_index()
    rows = []
    for keys, group in patient.groupby(['dataset', 'endpoint', 'family', 'signature_id'], observed=True):
        raw_stat, raw_p = signed_rank(group.patient_median_raw_rho)
        adj_stat, adj_p = signed_rank(group.patient_median_adjusted_rho)
        rows.append({
            'dataset': keys[0], 'endpoint': keys[1], 'family': keys[2], 'signature_id': keys[3],
            'n_patients': group.patient_key.nunique(),
            'n_patient_timepoints': int(group.n_timepoints.sum()),
            'total_cells': int(group.total_cells.sum()),
            'median_raw_rho': group.patient_median_raw_rho.median(),
            'raw_wilcoxon_statistic': raw_stat, 'raw_p_value': raw_p,
            'median_cycle_adjusted_rho': group.patient_median_adjusted_rho.median(),
            'adjusted_wilcoxon_statistic': adj_stat, 'adjusted_p_value': adj_p,
        })
    summary = pd.DataFrame(rows)
    summary['raw_q_value'] = summary.groupby(['dataset', 'endpoint', 'family'], observed=True)['raw_p_value'].transform(bh_adjust)
    summary['adjusted_q_value'] = summary.groupby(['dataset', 'endpoint', 'family'], observed=True)['adjusted_p_value'].transform(bh_adjust)
    summary = summary.merge(registry, on=['endpoint', 'family', 'signature_id'], how='left', validate='many_to_one')
    summary['observed_direction'] = np.where(summary.median_cycle_adjusted_rho > 0, 'positive', 'negative')
    summary['direction_matches_phase64'] = summary.observed_direction == summary.expected_phase64_direction
    summary['external_replication'] = summary.direction_matches_phase64 & (summary.adjusted_q_value < 0.05)
    return summary.sort_values(['dataset', 'endpoint', 'family', 'adjusted_q_value', 'signature_id']).reset_index(drop=True)


def build_cross_dataset(summary: pd.DataFrame, registry: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for signature in registry.itertuples(index=False):
        subset = summary.loc[summary.signature_id == signature.signature_id].copy()
        evaluations = {row.dataset: row for row in subset.itertuples(index=False)}
        replicated = [row for row in evaluations.values() if row.external_replication]
        concordant = [row for row in evaluations.values() if row.direction_matches_phase64]
        if len(evaluations) >= 2 and len(replicated) == len(evaluations):
            grade = 'R2_TWO_EXTERNAL_FDR'
        elif replicated and len(concordant) == len(evaluations):
            grade = 'R1_ONE_EXTERNAL_FDR_CONCORDANT'
        elif replicated:
            grade = 'R1_ONE_EXTERNAL_FDR_CONTEXT_CONFLICT'
        elif len(evaluations) == 1 and replicated:
            grade = 'R1_SINGLE_EXTERNAL_FDR'
        else:
            grade = 'R0_NOT_REPLICATED'
        row = {
            'endpoint': signature.endpoint, 'family': signature.family,
            'signature_id': signature.signature_id, 'label': signature.label,
            'expected_phase64_direction': signature.expected_phase64_direction,
            'datasets_evaluable': len(evaluations), 'datasets_fdr_replication': len(replicated),
            'replication_grade': grade,
        }
        for dataset in DATASETS:
            value = evaluations.get(dataset)
            row[f'{dataset}_n_patients'] = getattr(value, 'n_patients', np.nan) if value else np.nan
            row[f'{dataset}_median_adjusted_rho'] = getattr(value, 'median_cycle_adjusted_rho', np.nan) if value else np.nan
            row[f'{dataset}_adjusted_q_value'] = getattr(value, 'adjusted_q_value', np.nan) if value else np.nan
            row[f'{dataset}_direction_match'] = getattr(value, 'direction_matches_phase64', False) if value else False
        rows.append(row)
    return pd.DataFrame(rows).sort_values(['endpoint', 'family', 'replication_grade', 'signature_id']).reset_index(drop=True)


parser = argparse.ArgumentParser()
parser.add_argument('--h5ad', required=True, type=Path)
parser.add_argument('--metadata-dir', required=True, type=Path)
parser.add_argument('--input-dir', required=True, type=Path)
parser.add_argument('--cell-cycle-genes', required=True, type=Path)
parser.add_argument('--output', required=True, type=Path)
args = parser.parse_args()

part = args.output.with_name(args.output.name + '.part')
if args.output.exists() or part.exists():
    raise SystemExit('Refusing existing output directory')
part.mkdir(parents=True)
gate_log = part / 'memory_gate_v1.log'

try:
    memory_gate(gate_log, 'load_backed_h5ad')
    backed = sc.read_h5ad(args.h5ad, backed='r')
    if backed.n_obs != 258120 or backed.n_vars != 37795:
        raise RuntimeError(f'Unexpected H5AD dimensions {backed.shape}')
    obs_names = pd.Index(backed.obs_names.astype(str))
    if obs_names.has_duplicates:
        raise RuntimeError('Duplicate H5AD obs names')
    position = pd.Series(np.arange(len(obs_names), dtype=int), index=obs_names)
    source_var_names = pd.Index(backed.var_names.astype(str))
    upper_var_names = source_var_names.str.upper()
    collision_keys = set(upper_var_names[upper_var_names.duplicated(keep=False)])
    expected_collisions = {'C2ORF15', 'C6ORF165', 'ZHX1-C8ORF76'}
    if collision_keys != expected_collisions:
        raise RuntimeError(f'Unexpected uppercase gene collisions: {sorted(collision_keys)}')
    var_keep = ~upper_var_names.isin(collision_keys)
    var_positions = np.flatnonzero(var_keep)
    var_names = pd.Index(upper_var_names[var_keep])
    if var_names.has_duplicates or 'TYMP' not in var_names:
        raise RuntimeError('Invalid post-exclusion H5AD gene names or missing TYMP')

    registry, genes_by_signature = load_signatures(args.input_dir)
    cycle = pd.read_csv(args.cell_cycle_genes, sep='\t')
    s_genes = cycle.loc[cycle.phase == 'S', 'gene'].astype(str).str.upper().tolist()
    g2m_genes = cycle.loc[cycle.phase == 'G2M', 'gene'].astype(str).str.upper().tolist()
    if len(s_genes) != 43 or len(g2m_genes) != 54:
        raise RuntimeError('Unexpected frozen cell-cycle gene counts')

    write_tsv(registry, part / 'phase67_signature_registry_v1.tsv')
    correlations = []
    audit_rows = []
    score_index = 0
    sc.settings.n_jobs = 8
    np.random.seed(SEED)

    for dataset in DATASETS:
        metadata_path = args.metadata_dir / f'{dataset}_phase67_metadata_v1.tsv.gz'
        metadata = pd.read_csv(metadata_path, sep='\t', dtype=str)
        if metadata.cell_id.duplicated().any():
            raise RuntimeError(f'Duplicate metadata cell IDs: {dataset}')
        gsm = metadata.cell_id.str.extract(r'(GSM\d+)', expand=False)
        barcode = metadata.cell_id.str.rsplit('__', n=1).str[-1]
        metadata['h5ad_cell_id'] = gsm + '_' + barcode
        if metadata.h5ad_cell_id.isna().any() or metadata.h5ad_cell_id.duplicated().any():
            raise RuntimeError(f'Invalid derived H5AD IDs: {dataset}')
        missing = pd.Index(metadata.h5ad_cell_id).difference(obs_names)
        audit_rows.append({'dataset': dataset, 'endpoint': 'ALL', 'status': 'PASS' if len(missing) == 0 else 'FAIL', 'metadata_cells': len(metadata), 'matched_cells': len(metadata) - len(missing), 'eligible_groups': np.nan, 'note': f'missing={len(missing)}'})
        if len(missing):
            raise RuntimeError(f'Metadata/H5AD mismatch {dataset}: {len(missing)} missing')

        for endpoint in ENDPOINTS:
            selected_meta = metadata.loc[metadata.lineage_v1 == LINEAGE[endpoint]].copy()
            if selected_meta.empty:
                audit_rows.append({'dataset': dataset, 'endpoint': endpoint, 'status': 'UNEVALUABLE_NO_LINEAGE', 'metadata_cells': 0, 'matched_cells': 0, 'eligible_groups': 0, 'note': LINEAGE[endpoint]})
                continue
            memory_gate(gate_log, f'{dataset}_{endpoint}_score')
            selected_meta['h5ad_position'] = position.loc[selected_meta.h5ad_cell_id].to_numpy()
            selected_meta = selected_meta.sort_values('h5ad_position').reset_index(drop=True)
            matrix = backed.X[selected_meta.h5ad_position.to_numpy(), :][:, var_positions]
            matrix = sparse.csr_matrix(matrix)
            if matrix.data.size and (not np.isfinite(matrix.data).all() or matrix.data.min() < 0):
                raise RuntimeError(f'Invalid expression matrix {dataset} {endpoint}')
            subset = ad.AnnData(
                X=matrix,
                obs=selected_meta.set_index('h5ad_cell_id').drop(columns=['h5ad_position']),
                var=pd.DataFrame(index=var_names),
            )
            sc.tl.score_genes_cell_cycle(subset, s_genes=s_genes, g2m_genes=g2m_genes)
            tymp = np.asarray(subset[:, ['TYMP']].X.toarray()).ravel()
            subset.obs['TYMP_log1p'] = tymp

            endpoint_registry = registry.loc[registry.endpoint == endpoint]
            coverage = []
            score_columns = []
            for signature in endpoint_registry.itertuples(index=False):
                present = [gene for gene in genes_by_signature[signature.signature_id] if gene in var_names]
                present = list(dict.fromkeys(present))
                eligible = len(present) >= 10
                coverage.append({
                    'dataset': dataset, 'endpoint': endpoint, 'family': signature.family,
                    'signature_id': signature.signature_id, 'n_frozen_genes': signature.n_frozen_genes,
                    'n_genes_present': len(present), 'eligible': eligible,
                })
                if not eligible:
                    continue
                score_name = f'SIG_{score_index:03d}'
                score_index += 1
                sc.tl.score_genes(
                    subset, gene_list=present, ctrl_size=max(20, min(50, len(present))),
                    n_bins=25, score_name=score_name, random_state=SEED,
                    ctrl_as_ref=False, use_raw=False,
                )
                score_columns.append((signature.signature_id, signature.family, score_name))
            coverage_frame = pd.DataFrame(coverage)
            coverage_frame.to_csv(part / f'{dataset}_{endpoint}_signature_coverage_v1.tsv', sep='\t', index=False, lineterminator='\n')

            eligible_groups = 0
            for (patient_key, timepoint), group in subset.obs.groupby(['patient_key', 'timepoint'], observed=True):
                if len(group) < MIN_CELLS:
                    continue
                predictor = group.TYMP_log1p.to_numpy(float)
                if np.unique(predictor).size < 3:
                    continue
                s_score = group.S_score.to_numpy(float)
                g2m_score = group.G2M_score.to_numpy(float)
                predictor_residual = residualize(predictor, s_score, g2m_score)
                eligible_groups += 1
                for signature_id, family, score_name in score_columns:
                    score = group[score_name].to_numpy(float)
                    if np.unique(score).size < 3:
                        continue
                    score_residual = residualize(score, s_score, g2m_score)
                    raw_rho = stats.spearmanr(predictor, score).statistic
                    adjusted_rho = stats.spearmanr(predictor_residual, score_residual).statistic
                    correlations.append({
                        'dataset': dataset, 'endpoint': endpoint, 'family': family,
                        'signature_id': signature_id, 'patient_key': patient_key,
                        'patient_id': group.patient_id.iloc[0], 'timepoint': timepoint,
                        'treatment': group.treatment.iloc[0], 'response': group.response.iloc[0],
                        'n_cells': len(group), 'raw_spearman_rho': raw_rho,
                        'cycle_adjusted_spearman_rho': adjusted_rho,
                    })
            audit_rows.append({'dataset': dataset, 'endpoint': endpoint, 'status': 'PASS', 'metadata_cells': len(selected_meta), 'matched_cells': len(selected_meta), 'eligible_groups': eligible_groups, 'note': f'signatures_scored={len(score_columns)}'})
            del subset, matrix, selected_meta
            gc.collect()

    correlation_frame = pd.DataFrame(correlations)
    if correlation_frame.empty:
        raise RuntimeError('No patient-timepoint correlations produced')
    summary = build_summary(correlation_frame, registry)
    cross = build_cross_dataset(summary, registry)
    retained = cross.loc[cross.replication_grade != 'R0_NOT_REPLICATED'].copy()
    audit = pd.DataFrame(audit_rows)

    write_tsv(audit, part / 'phase67_input_and_lineage_audit_v1.tsv')
    write_tsv(correlation_frame.sort_values(['dataset', 'endpoint', 'family', 'signature_id', 'patient_key', 'timepoint']), part / 'phase67_correlations_by_patient_timepoint_v1.tsv')
    write_tsv(summary, part / 'phase67_external_replication_summary_v1.tsv')
    write_tsv(cross, part / 'phase67_cross_dataset_replication_v1.tsv')
    write_tsv(retained, part / 'phase67_cross_dataset_retained_v1.tsv')

    grade_counts = cross.replication_grade.value_counts().to_dict()
    record = [
        '# Phase67 direct TYMP-program external replication v1', '',
        '## Design', '',
        '- Frozen Phase64 regulon and TYMP-excluded Reactome signatures were scored with the same Scanpy matched-control method used in Phase64.',
        '- TYMP and every program score were residualized against S and G2M within patient/timepoint and lineage.',
        '- Spearman correlations were calculated within patient/timepoint; inference used patient-level median correlations and BH within dataset, endpoint and signature family.',
        '- GSE294109 regions are pooled through the shared patient/timepoint unit. At least 50 cells are required per unit.', '',
        '## Scope', '',
        f'- Patient-timepoint correlation rows: {len(correlation_frame):,}.',
        f'- External summary rows: {len(summary):,}; cross-dataset signature rows: {len(cross):,}.',
        f'- Replication grades: {json.dumps(grade_counts, sort_keys=True)}.', '',
        '## Boundary', '',
        'A replicated TYMP-program correlation strengthens lineage-specific mechanism evidence. It does not by itself validate response, PFS or OS prediction, and Treg is evaluable only in GSE294109 because GSE264586 lacks a frozen Treg fine annotation.',
    ]
    (part / 'PHASE67_REPLICATION_RECORD_v1.md').write_text('\n'.join(record) + '\n', encoding='utf-8')
    status = {
        'status': 'PASS_PHASE67_DIRECT_REPLICATION_V1',
        'completed_at': now(),
        'h5ad_cells': backed.n_obs,
        'h5ad_genes': backed.n_vars,
        'analysis_genes_after_collision_exclusion': len(var_names),
        'uppercase_collision_keys_excluded': sorted(collision_keys),
        'source_columns_excluded': int((~var_keep).sum()),
        'metadata_cells': int(sum(len(pd.read_csv(args.metadata_dir / f'{dataset}_phase67_metadata_v1.tsv.gz', sep='\t')) for dataset in DATASETS)),
        'correlation_rows': len(correlation_frame),
        'external_summary_rows': len(summary),
        'cross_dataset_rows': len(cross),
        'retained_rows': len(retained),
        'replication_grade_counts': grade_counts,
    }
    (part / 'status_v1.json').write_text(json.dumps(status, indent=2) + '\n', encoding='utf-8')
    manifest_rows = []
    for path in sorted(part.iterdir()):
        if path.is_file():
            manifest_rows.append({'file': path.name, 'size_bytes': path.stat().st_size, 'sha256': sha256(path)})
    write_tsv(pd.DataFrame(manifest_rows), part / 'FILE_MANIFEST_SHA256_v1.tsv')
    (part / 'RUN_COMPLETE').write_text('PASS\n', encoding='utf-8')
    backed.file.close()
    os.replace(part, args.output)
    print(json.dumps(status, indent=2))
except Exception:
    (part / 'RUN_FAILURE_v1.txt').write_text(traceback.format_exc(), encoding='utf-8')
    raise
