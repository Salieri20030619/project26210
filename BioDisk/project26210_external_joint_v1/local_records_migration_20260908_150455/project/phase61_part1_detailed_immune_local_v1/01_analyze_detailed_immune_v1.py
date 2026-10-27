#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


IMMUNE_BROAD = ['T_cell', 'NK_cell', 'Myeloid']


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    part = path.with_name(path.name + '.part')
    frame.to_csv(part, sep='\t', index=False, na_rep='NA', lineterminator='\n')
    os.replace(part, path)


def atomic_json(payload: dict, path: Path) -> None:
    part = path.with_name(path.name + '.part')
    part.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    os.replace(part, path)


def bh_adjust(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.notna()
    p = values.loc[valid].astype(float)
    if p.empty:
        return result
    order = np.argsort(p.to_numpy())
    ranked = p.to_numpy()[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    original = np.empty_like(adjusted)
    original[order] = adjusted
    result.loc[p.index] = original
    return result


def composition_features(comp: pd.DataFrame, unit: str) -> pd.DataFrame:
    current = comp[comp.unit == unit].copy()
    if unit == 'sample_region':
        keys = [
            'dataset', 'patient_key', 'patient_id', 'sample_id', 'object_name',
            'timepoint', 'treatment', 'response',
        ]
    else:
        keys = ['dataset', 'patient_key', 'patient_id', 'timepoint', 'treatment', 'response']
    base = current[keys].drop_duplicates()
    total = current.groupby(keys, dropna=False).n_cells.sum().rename('total_cells').reset_index()
    immune = current[current.broad_annotation.isin(IMMUNE_BROAD)].copy()
    immune_total = immune.groupby(keys, dropna=False).n_cells.sum().rename('immune_cells').reset_index()
    base = base.merge(total, on=keys, validate='one_to_one').merge(immune_total, on=keys, how='left', validate='one_to_one')
    base['immune_cells'] = base.immune_cells.fillna(0)

    frames = []
    overall = base.copy()
    overall['feature_family'] = 'immune_abundance'
    overall['feature_id'] = 'immune_fraction_of_all_cells'
    overall['numerator_cells'] = overall.immune_cells
    overall['denominator_cells'] = overall.total_cells
    overall['value'] = overall.numerator_cells / overall.denominator_cells
    frames.append(overall[keys + ['feature_family', 'feature_id', 'numerator_cells', 'denominator_cells', 'value']])

    for family, column, prefix in [
        ('broad_composition', 'broad_annotation', 'broad_fraction_within_immune__'),
        ('fine_composition', 'fine_annotation', 'fine_fraction_within_immune__'),
    ]:
        catalog = sorted(immune[column].dropna().unique())
        counts = immune.groupby(keys + [column], dropna=False).n_cells.sum().rename('numerator_cells').reset_index()
        grid = base[keys + ['immune_cells']].assign(_join=1).merge(
            pd.DataFrame({column: catalog, '_join': 1}), on='_join', validate='many_to_many'
        ).drop(columns='_join')
        grid = grid.merge(counts, on=keys + [column], how='left', validate='one_to_one')
        grid['numerator_cells'] = grid.numerator_cells.fillna(0)
        grid['denominator_cells'] = grid.immune_cells
        grid['value'] = np.where(grid.denominator_cells > 0, grid.numerator_cells / grid.denominator_cells, np.nan)
        grid['feature_family'] = family
        grid['feature_id'] = prefix + grid[column].astype(str)
        frames.append(grid[keys + ['feature_family', 'feature_id', 'numerator_cells', 'denominator_cells', 'value']])
    result = pd.concat(frames, ignore_index=True)
    result['unit'] = unit
    return result


def response_comparison(features: pd.DataFrame) -> pd.DataFrame:
    source = features[
        (features.dataset == 'GSE294109')
        & (features.unit == 'patient_timepoint')
        & (features.timepoint == 'Surgery_after_treatment')
        & (features.treatment != 'Cabozantinib')
    ].copy()
    source['response_group'] = np.select(
        [source.response.isin(['CR', 'PR']), source.response.eq('No_effect')],
        ['benefit', 'no_effect'], default=None,
    )
    source = source[source.response_group.notna() & source.value.notna()]
    rows = []
    for (family, feature), group in source.groupby(['feature_family', 'feature_id'], sort=True):
        benefit = group.loc[group.response_group == 'benefit', 'value'].to_numpy(float)
        no_effect = group.loc[group.response_group == 'no_effect', 'value'].to_numpy(float)
        if len(benefit) and len(no_effect):
            test = stats.mannwhitneyu(benefit, no_effect, alternative='two-sided', method='asymptotic')
            p_value = float(test.pvalue)
            rank_biserial = float(2 * test.statistic / (len(benefit) * len(no_effect)) - 1)
        else:
            p_value = np.nan
            rank_biserial = np.nan
        rows.append({
            'feature_family': family,
            'feature_id': feature,
            'n_benefit': len(benefit),
            'n_no_effect': len(no_effect),
            'median_benefit': np.median(benefit) if len(benefit) else np.nan,
            'median_no_effect': np.median(no_effect) if len(no_effect) else np.nan,
            'median_difference_benefit_minus_no_effect': np.median(benefit) - np.median(no_effect) if len(benefit) and len(no_effect) else np.nan,
            'rank_biserial_benefit_vs_no_effect': rank_biserial,
            'p_value': p_value,
        })
    result = pd.DataFrame(rows)
    result['p_fdr_bh_within_family'] = result.groupby('feature_family', group_keys=False).p_value.apply(bh_adjust)
    result['reporting_status'] = np.where(
        (result.n_benefit == 9) & (result.n_no_effect == 3), 'FULL_9_VS_3', 'AVAILABLE_CASE'
    )
    return result.sort_values(['feature_family', 'p_value', 'feature_id'])


def longitudinal_analysis(features: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    source = features[(features.dataset == 'GSE264586') & (features.unit == 'patient_timepoint')].copy()
    patient_rows = []
    summary_rows = []
    for followup in ['Cycle2', 'End_of_treatment']:
        before = source[source.timepoint == 'Baseline'][
            ['patient_key', 'patient_id', 'feature_family', 'feature_id', 'value']
        ].rename(columns={'value': 'value_baseline'})
        after = source[source.timepoint == followup][
            ['patient_key', 'patient_id', 'feature_family', 'feature_id', 'value']
        ].rename(columns={'value': 'value_followup'})
        paired = before.merge(
            after, on=['patient_key', 'patient_id', 'feature_family', 'feature_id'], validate='one_to_one'
        ).dropna(subset=['value_baseline', 'value_followup'])
        paired['comparison'] = followup + '_vs_Baseline'
        paired['delta_followup_minus_baseline'] = paired.value_followup - paired.value_baseline
        patient_rows.append(paired)
        for (family, feature), group in paired.groupby(['feature_family', 'feature_id'], sort=True):
            delta = group.delta_followup_minus_baseline.to_numpy(float)
            if len(delta) >= 2 and np.any(delta != 0):
                p_value = float(stats.wilcoxon(group.value_followup, group.value_baseline, alternative='two-sided', method='approx').pvalue)
            elif len(delta) >= 2:
                p_value = 1.0
            else:
                p_value = np.nan
            summary_rows.append({
                'comparison': followup + '_vs_Baseline',
                'feature_family': family,
                'feature_id': feature,
                'n_pairs': len(group),
                'median_baseline': group.value_baseline.median(),
                'median_followup': group.value_followup.median(),
                'median_delta_followup_minus_baseline': group.delta_followup_minus_baseline.median(),
                'p_value': p_value,
            })
    patient = pd.concat(patient_rows, ignore_index=True)
    summary = pd.DataFrame(summary_rows)
    summary['p_fdr_bh_within_comparison_family'] = summary.groupby(
        ['comparison', 'feature_family'], group_keys=False
    ).p_value.apply(bh_adjust)
    return patient.sort_values(['comparison', 'feature_family', 'feature_id', 'patient_id']), summary.sort_values(
        ['comparison', 'feature_family', 'p_value', 'feature_id']
    )


def multiregion_summary(features: pd.DataFrame) -> pd.DataFrame:
    source = features[(features.dataset == 'GSE294109') & (features.unit == 'sample_region')].dropna(subset=['value']).copy()
    rows = []
    for (patient_key, patient_id, timepoint, treatment, response, family, feature), group in source.groupby(
        ['patient_key', 'patient_id', 'timepoint', 'treatment', 'response', 'feature_family', 'feature_id'], sort=True
    ):
        if group.sample_id.nunique() < 2:
            continue
        values = group.value.to_numpy(float)
        rows.append({
            'patient_key': patient_key,
            'patient_id': patient_id,
            'timepoint': timepoint,
            'treatment': treatment,
            'response': response,
            'feature_family': family,
            'feature_id': feature,
            'n_regions': group.sample_id.nunique(),
            'region_mean': values.mean(),
            'region_sd': values.std(ddof=1),
            'region_minimum': values.min(),
            'region_maximum': values.max(),
            'region_range': values.max() - values.min(),
        })
    return pd.DataFrame(rows).sort_values(['patient_id', 'feature_family', 'feature_id'])


def rho_associations(features: pd.DataFrame, correlations: pd.DataFrame) -> pd.DataFrame:
    rho = correlations[
        (correlations.unit == 'patient_timepoint')
        & (correlations.method == 'primary_gene_z')
        & (correlations.status == 'ELIGIBLE_GE20')
    ].copy()
    rho = rho[
        ((rho.dataset == 'GSE264586') & (rho.timepoint == 'Baseline'))
        | ((rho.dataset == 'GSE294109') & (rho.timepoint == 'Surgery_after_treatment') & (rho.treatment != 'Cabozantinib'))
    ]
    immune = features[features.unit == 'patient_timepoint'].dropna(subset=['value'])
    merged = rho.merge(
        immune,
        on=['dataset', 'patient_key', 'patient_id', 'timepoint'],
        suffixes=('_rho', '_immune'), validate='many_to_many',
    )
    rows = []
    for (dataset, tier, family, feature), group in merged.groupby(
        ['dataset', 'malignancy_tier', 'feature_family', 'feature_id'], sort=True
    ):
        x = group.spearman_rho.to_numpy(float)
        y = group.value.to_numpy(float)
        if len(group) >= 5 and np.unique(x).size > 1 and np.unique(y).size > 1:
            test = stats.spearmanr(x, y)
            rho_value, p_value = float(test.statistic), float(test.pvalue)
            testing_status = 'TESTED_GE5'
        elif len(group) >= 3 and np.unique(x).size > 1 and np.unique(y).size > 1:
            rho_value = float(stats.spearmanr(x, y).statistic)
            p_value = np.nan
            testing_status = 'DIRECTION_ONLY_LT5'
        else:
            rho_value, p_value = np.nan, np.nan
            testing_status = 'UNTESTABLE'
        rows.append({
            'dataset': dataset,
            'malignancy_tier': tier,
            'feature_family': family,
            'feature_id': feature,
            'n_patients': len(group),
            'spearman_rho': rho_value,
            'p_value': p_value,
            'testing_status': testing_status,
        })
    result = pd.DataFrame(rows)
    result['p_fdr_bh_within_dataset_tier_family'] = result.groupby(
        ['dataset', 'malignancy_tier', 'feature_family'], group_keys=False
    ).p_value.apply(bh_adjust)
    return result.sort_values(['dataset', 'malignancy_tier', 'feature_family', 'p_value', 'feature_id'])


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit('Usage: 01_analyze_detailed_immune_v1.py <phase60_results> <output_dir>')
    phase60 = Path(sys.argv[1]).resolve()
    out = Path(sys.argv[2]).resolve()
    if out.exists():
        raise FileExistsError(out)
    out.mkdir(parents=True)

    comp = pd.read_csv(phase60 / 'immune_composition_by_unit_v1.tsv', sep='\t')
    states = pd.read_csv(phase60 / 'immune_state_scores_patient_timepoint_v1.tsv', sep='\t')
    correlations = pd.read_csv(phase60 / 'patient_f1_tymp_correlations_v1.tsv', sep='\t')

    composition = pd.concat(
        [composition_features(comp, 'sample_region'), composition_features(comp, 'patient_timepoint')],
        ignore_index=True,
    )
    state_features = states.rename(columns={'mean_score': 'value'}).copy()
    state_features['feature_family'] = 'immune_state'
    state_features['feature_id'] = state_features.analysis_compartment + '__' + state_features.score_id
    state_features['unit'] = 'patient_timepoint'
    state_features['numerator_cells'] = state_features.n_cells
    state_features['denominator_cells'] = state_features.n_cells
    state_features = state_features[
        ['dataset', 'patient_key', 'patient_id', 'timepoint', 'treatment', 'response', 'unit',
         'feature_family', 'feature_id', 'numerator_cells', 'denominator_cells', 'value']
    ]
    features = pd.concat([composition, state_features], ignore_index=True, sort=False)
    features = features.sort_values(['dataset', 'unit', 'patient_id', 'timepoint', 'feature_family', 'feature_id'])
    atomic_tsv(features, out / 'immune_features_long_v1.tsv')

    response = response_comparison(features)
    atomic_tsv(response, out / 'GSE294109_icb_immune_response_comparison_v1.tsv')
    long_patient, long_summary = longitudinal_analysis(features)
    atomic_tsv(long_patient, out / 'GSE264586_immune_longitudinal_patient_deltas_v1.tsv')
    atomic_tsv(long_summary, out / 'GSE264586_immune_longitudinal_summary_v1.tsv')
    multiregion = multiregion_summary(features)
    atomic_tsv(multiregion, out / 'GSE294109_multiregion_immune_heterogeneity_v1.tsv')
    associations = rho_associations(features, correlations)
    atomic_tsv(associations, out / 'immune_feature_vs_f1_tymp_rho_associations_v1.tsv')

    status = {
        'status': 'PASS_PHASE61_DETAILED_IMMUNE',
        'completed_at': datetime.now().astimezone().isoformat(timespec='seconds'),
        'feature_rows': len(features),
        'response_tests': len(response),
        'longitudinal_tests': len(long_summary),
        'multiregion_rows': len(multiregion),
        'rho_association_tests': len(associations),
        'response_unit': 'patient after pooling regions',
        'longitudinal_unit': 'paired patient',
        'multiple_testing': 'BH within analysis family',
    }
    atomic_json(status, out / 'PHASE61_STATUS_v1.json')
    print(json.dumps(status, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
