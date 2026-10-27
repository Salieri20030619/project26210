import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(r'D:\project26210_external_joint_v1\phase114_unweighted_seven_feature_analysis_local_v1')
DOWNLOADS = ROOT / 'downloads_v1'
INPUTS = ROOT / 'inputs_v1'
RESULTS = ROOT / 'results_v1'
FIGURES = ROOT / 'figures_v1'
SEED = 262114
N_BOOT = 10000
N_PERM = 100000

DATASETS = {
    'RC18': ('matched_Harmonized_RawData_RC18.csv', 'matched_tpm_RC18.csv'),
    'RC20': ('matched_Harmonized_RawData_RC20.csv', 'matched_tpm_RC20.csv'),
    'CPTAC': ('matched_Harmonized_Met_CPTAC.csv', 'matched_tpm_CPTAC.csv'),
    'CPTAC_val': ('matched_Harmonized_Met_CPTAC_val.csv', 'matched_tpm_CPTAC_val.csv'),
}

FEATURES = {
    '5-MTA': ('5-methylthioadenosine (MTA)', -1),
    'Thymidine': ('thymidine', 1),
    "2'-Deoxyuridine": ("2'-deoxyuridine", 1),
    'Methionine': ('methionine', -1),
    'SAM': ('S-adenosylmethionine (SAM)', -1),
    'Thymine': ('thymine', -1),
    'Spermidine': ('spermidine', -1),
}

FULL_FEATURES = list(FEATURES)
CPTAC_FEATURES = ['5-MTA', 'Methionine', 'Thymine']
NCI_FEATURES = ['5-MTA', "2'-Deoxyuridine", 'Methionine', 'SAM', 'Spermidine']
NCI_ORIGINAL = {
    '5-MTA': "5'-s-methyl-5'-thioadenosine",
    "2'-Deoxyuridine": "2'-deoxyuridine",
    'Methionine': 'methionine',
    'SAM': "S-(5'-adenosyl)-L-methionine",
    'Spermidine': 'Spermidine',
}


def rank_z(values):
    ranks = stats.rankdata(np.asarray(values, dtype=float), method='average')
    return (ranks - ranks.mean()) / ranks.std(ddof=0)


def raw_z(values):
    values = np.asarray(values, dtype=float)
    return (values - values.mean()) / values.std(ddof=0)


def bh_adjust(values):
    p = np.asarray(values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def bootstrap_spearman(x, y, rng, n_boot=N_BOOT):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    draws = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        rho = stats.spearmanr(x[idx], y[idx]).statistic
        if np.isfinite(rho):
            draws.append(float(rho))
    if len(draws) < n_boot * 0.95:
        raise RuntimeError('too many invalid bootstrap correlations')
    return float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975)), len(draws)


def partial_rank_corr(x, y, covariates):
    xr = rank_z(x)
    yr = rank_z(y)
    covariates = np.asarray(covariates, dtype=float)
    if covariates.ndim == 1:
        covariates = covariates[:, None]
    design = np.column_stack([np.ones(len(xr)), covariates])
    xres = xr - design @ np.linalg.lstsq(design, xr, rcond=None)[0]
    yres = yr - design @ np.linalg.lstsq(design, yr, rcond=None)[0]
    rho, p = stats.pearsonr(xres, yres)
    return float(rho), float(p), xres, yres


def permutation_spearman(x, y, rng, n_perm=N_PERM, chunk=5000):
    xr = stats.rankdata(np.asarray(x, dtype=float), method='average')
    yr = stats.rankdata(np.asarray(y, dtype=float), method='average')
    xr = xr - xr.mean()
    yr = yr - yr.mean()
    denom = math.sqrt(float(np.sum(xr ** 2) * np.sum(yr ** 2)))
    observed = float(np.dot(xr, yr) / denom)
    exceed = 0
    completed = 0
    while completed < n_perm:
        current = min(chunk, n_perm - completed)
        order = np.argsort(rng.random((current, len(xr))), axis=1)
        permuted = xr[order]
        values = permuted @ yr / denom
        exceed += int(np.sum(np.abs(values) >= abs(observed) - 1e-15))
        completed += current
    return observed, (exceed + 1.0) / (n_perm + 1.0)


def within_group_permutation(xres, yres, groups, rng, n_perm=N_PERM, chunk=2000):
    xres = np.asarray(xres, dtype=float)
    yres = np.asarray(yres, dtype=float)
    groups = np.asarray(groups)
    denom = math.sqrt(float(np.sum(xres ** 2) * np.sum(yres ** 2)))
    observed = float(np.dot(xres, yres) / denom)
    indices = [np.flatnonzero(groups == group) for group in np.unique(groups)]
    exceed = 0
    completed = 0
    while completed < n_perm:
        current = min(chunk, n_perm - completed)
        matrix = np.tile(xres, (current, 1))
        for idx in indices:
            order = np.argsort(rng.random((current, len(idx))), axis=1)
            matrix[:, idx] = xres[idx][order]
        values = matrix @ yres / denom
        exceed += int(np.sum(np.abs(values) >= abs(observed) - 1e-15))
        completed += current
    return observed, (exceed + 1.0) / (n_perm + 1.0)


def score_frame(frame, feature_names, score_name):
    scored = frame.copy()
    signed_columns = []
    raw_signed_columns = []
    for feature in feature_names:
        sign = FEATURES[feature][1]
        z_col = f'zrank__{feature}'
        signed_col = f'signed_zrank__{feature}'
        raw_col = f'zraw__{feature}'
        raw_signed_col = f'signed_zraw__{feature}'
        scored[z_col] = rank_z(scored[feature])
        scored[signed_col] = sign * scored[z_col]
        scored[raw_col] = raw_z(scored[feature])
        scored[raw_signed_col] = sign * scored[raw_col]
        signed_columns.append(signed_col)
        raw_signed_columns.append(raw_signed_col)
    scored[score_name] = scored[signed_columns].mean(axis=1)
    scored[f'{score_name}_raw_z'] = scored[raw_signed_columns].mean(axis=1)
    scored[f'{score_name}_rank_for_spearman'] = stats.rankdata(scored[score_name], method='average')
    return scored


def prepare_unitedmet(dataset, feature_names, score_name):
    met_name, rna_name = DATASETS[dataset]
    mapping = pd.read_csv(INPUTS / 'MasterMapping_ccRCC.csv', dtype=str).fillna('')
    mapd = mapping.loc[mapping['Dataset'].eq(dataset) & mapping['TN'].eq('Tumor')].copy()
    met = pd.read_csv(DOWNLOADS / met_name, index_col=0)
    rna = pd.read_csv(DOWNLOADS / rna_name, index_col=0)
    met.index = met.index.astype(str)
    rna.index = rna.index.astype(str)
    met = met.loc[mapd['MetabID']]
    rna = rna.loc[mapd['RNAID']]
    common_ids = mapd['CommonID'].where(mapd['CommonID'].ne(''), mapd['MetabID']).to_numpy()
    frame = pd.DataFrame(index=common_ids)
    frame.index.name = 'sample_id'
    for feature in feature_names:
        frame[feature] = pd.to_numeric(met[FEATURES[feature][0]], errors='coerce').to_numpy()
    frame['TYMP'] = np.log2(pd.to_numeric(rna['TYMP'], errors='coerce').to_numpy() + 1.0)
    frame['MTAP'] = np.log2(pd.to_numeric(rna['MTAP'], errors='coerce').to_numpy() + 1.0)
    frame = frame.replace([np.inf, -np.inf], np.nan).dropna()
    scored = score_frame(frame, feature_names, score_name)
    scored.insert(0, 'dataset', dataset)
    return scored


def meta_correlation(rows):
    rho = np.clip(np.asarray([row['rho'] for row in rows], dtype=float), -0.999999, 0.999999)
    n = np.asarray([row['n'] for row in rows], dtype=float)
    z = np.arctanh(rho)
    var = 1.0 / (n - 3.0)
    weights = 1.0 / var
    fixed_z = float(np.sum(weights * z) / np.sum(weights))
    fixed_se = math.sqrt(1.0 / float(np.sum(weights)))
    q_stat = float(np.sum(weights * (z - fixed_z) ** 2))
    df = len(rows) - 1
    q_p = float(stats.chi2.sf(q_stat, df))
    c_value = float(np.sum(weights) - np.sum(weights ** 2) / np.sum(weights))
    tau2 = max(0.0, (q_stat - df) / c_value) if c_value > 0 else 0.0
    random_weights = 1.0 / (var + tau2)
    random_z = float(np.sum(random_weights * z) / np.sum(random_weights))
    random_se = math.sqrt(1.0 / float(np.sum(random_weights)))
    i2 = max(0.0, (q_stat - df) / q_stat) * 100.0 if q_stat > 0 else 0.0
    return {
        'k': len(rows),
        'total_complete_case_n': int(np.sum(n)),
        'fixed_rho': math.tanh(fixed_z),
        'fixed_ci_low': math.tanh(fixed_z - 1.96 * fixed_se),
        'fixed_ci_high': math.tanh(fixed_z + 1.96 * fixed_se),
        'fixed_p': float(2 * stats.norm.sf(abs(fixed_z / fixed_se))),
        'random_rho': math.tanh(random_z),
        'random_ci_low': math.tanh(random_z - 1.96 * random_se),
        'random_ci_high': math.tanh(random_z + 1.96 * random_se),
        'random_p': float(2 * stats.norm.sf(abs(random_z / random_se))),
        'q_stat': q_stat,
        'q_p': q_p,
        'i2_percent': i2,
        'tau2': tau2,
        'directions': ';'.join('negative' if value < 0 else 'positive' if value > 0 else 'zero' for value in rho),
    }


def save_tsv(frame, name):
    frame.to_csv(RESULTS / name, sep='\t', index=False, float_format='%.17g')


def analyse_unitedmet_primary():
    rng = np.random.default_rng(SEED)
    score_frames = []
    primary = []
    components = []
    partials = []
    leave_feature = []
    leave_sample = []
    permutations = []
    raw_sensitivity = []
    contributions = []
    correlation_long = []

    for dataset in ['RC18', 'RC20']:
        frame = prepare_unitedmet(dataset, FULL_FEATURES, 'S7')
        score_frames.append(frame.reset_index())
        rho, p = stats.spearmanr(frame['TYMP'], frame['S7_rank_for_spearman'])
        low, high, valid_boot = bootstrap_spearman(frame['TYMP'], frame['S7_rank_for_spearman'], rng)
        primary.append({
            'dataset': dataset, 'n': len(frame), 'rho': float(rho), 'bootstrap_ci_low': low,
            'bootstrap_ci_high': high, 'valid_bootstrap_draws': valid_boot, 'p': float(p),
            'expected_direction': 'negative',
        })
        partial_rho, partial_p, _, _ = partial_rank_corr(frame['TYMP'], frame['S7_rank_for_spearman'], rank_z(frame['MTAP']))
        partials.append({'dataset': dataset, 'n': len(frame), 'partial_rho_controlling_MTAP': partial_rho, 'p': partial_p})

        component_rows = []
        for feature in FULL_FEATURES:
            value_rho, value_p = stats.spearmanr(frame['TYMP'], frame[feature])
            expected = -FEATURES[feature][1]
            component_rows.append({
                'dataset': dataset, 'feature': feature, 'n': len(frame), 'rho': float(value_rho),
                'p': float(value_p), 'expected_rho_sign': int(expected),
                'expected_direction_hit': bool(np.sign(value_rho) == expected),
            })
        q_values = bh_adjust([row['p'] for row in component_rows])
        for row, q_value in zip(component_rows, q_values):
            row['q_within_cohort_7'] = float(q_value)
            components.append(row)

        for omitted in FULL_FEATURES:
            retained = [feature for feature in FULL_FEATURES if feature != omitted]
            temp = score_frame(frame[retained + ['TYMP', 'MTAP']].copy(), retained, 'S6')
            value_rho, value_p = stats.spearmanr(temp['TYMP'], temp['S6_rank_for_spearman'])
            leave_feature.append({
                'dataset': dataset, 'omitted_feature': omitted, 'n': len(temp),
                'rho': float(value_rho), 'p_descriptive': float(value_p),
                'direction_negative': bool(value_rho < 0),
            })

        for sample_id in frame.index:
            subset = frame.drop(index=sample_id)
            raw = subset[FULL_FEATURES + ['TYMP', 'MTAP']].copy()
            temp = score_frame(raw, FULL_FEATURES, 'S7_LOO')
            value_rho = stats.spearmanr(temp['TYMP'], temp['S7_LOO_rank_for_spearman']).statistic
            leave_sample.append({
                'dataset': dataset, 'omitted_sample_id': sample_id, 'n': len(temp),
                'rho': float(value_rho), 'direction_negative': bool(value_rho < 0),
            })

        perm_rho, perm_p = permutation_spearman(frame['TYMP'], frame['S7_rank_for_spearman'], rng)
        permutations.append({'dataset': dataset, 'n': len(frame), 'rho': perm_rho, 'permutations': N_PERM, 'permutation_p': perm_p})
        raw_rho, raw_p = stats.spearmanr(frame['TYMP'], frame['S7_raw_z'])
        raw_sensitivity.append({'dataset': dataset, 'n': len(frame), 'rho_raw_z_score': float(raw_rho), 'p': float(raw_p)})

        signed = [f'signed_zrank__{feature}' for feature in FULL_FEATURES]
        matrix = frame[signed]
        score_var = float(np.var(frame['S7'], ddof=0))
        for feature, column in zip(FULL_FEATURES, signed):
            cov = float(np.cov(matrix[column] / len(FULL_FEATURES), frame['S7'], ddof=0)[0, 1])
            contributions.append({
                'dataset': dataset, 'feature': feature,
                'correlation_with_total_score': float(stats.pearsonr(matrix[column], frame['S7']).statistic),
                'variance_contribution_fraction': cov / score_var,
            })
        corr = matrix.corr(method='spearman')
        for row_feature, row_col in zip(FULL_FEATURES, signed):
            for col_feature, col_col in zip(FULL_FEATURES, signed):
                correlation_long.append({
                    'dataset': dataset, 'feature_1': row_feature, 'feature_2': col_feature,
                    'signed_spearman_rho': float(corr.loc[row_col, col_col]),
                })

    q_primary = bh_adjust([row['p'] for row in primary])
    for row, q_value in zip(primary, q_primary):
        row['q_primary_family_2'] = float(q_value)
    q_partial = bh_adjust([row['p'] for row in partials])
    for row, q_value in zip(partials, q_partial):
        row['q_partial_family_2'] = float(q_value)
    meta = meta_correlation(primary)
    both_negative = all(row['rho'] < 0 for row in primary)
    both_q = all(row['q_primary_family_2'] < 0.05 for row in primary)
    if both_negative and both_q:
        decision = 'COMPLETE_MIGRATION_SUPPORT'
    elif both_negative and meta['fixed_rho'] < 0 and meta['fixed_p'] < 0.05:
        decision = 'DIRECTIONALLY_CONSISTENT_PARTIAL_SUPPORT'
    elif both_negative:
        decision = 'EXPLORATORY_DIRECTIONAL_CONSISTENCY'
    else:
        decision = 'CURRENT_DATA_DO_NOT_SUPPORT_GENERAL_SEVEN_FEATURE_TRANSFER'

    save_tsv(pd.concat(score_frames, ignore_index=True), 'unitedmet_full_score_sample_scores_v1.tsv')
    save_tsv(pd.DataFrame(primary), 'full_score_primary_associations_v1.tsv')
    save_tsv(pd.DataFrame(components), 'full_score_component_associations_v1.tsv')
    save_tsv(pd.DataFrame(partials), 'full_score_mtap_partial_v1.tsv')
    save_tsv(pd.DataFrame([meta]), 'full_score_meta_analysis_v1.tsv')
    save_tsv(pd.DataFrame(leave_feature), 'robustness_leave_one_feature_v1.tsv')
    save_tsv(pd.DataFrame(leave_sample), 'robustness_leave_one_sample_v1.tsv')
    save_tsv(pd.DataFrame(permutations), 'robustness_permutation_v1.tsv')
    save_tsv(pd.DataFrame(raw_sensitivity), 'robustness_raw_z_sensitivity_v1.tsv')
    save_tsv(pd.DataFrame(contributions), 'feature_contribution_audit_v1.tsv')
    save_tsv(pd.DataFrame(correlation_long), 'signed_feature_correlation_long_v1.tsv')
    decision_record = {
        'decision': decision,
        'primary_complete_case_n': {row['dataset']: row['n'] for row in primary},
        'both_cohort_rho_negative': both_negative,
        'both_primary_q_below_0_05': both_q,
        'fixed_meta_rho': meta['fixed_rho'],
        'fixed_meta_p': meta['fixed_p'],
        'conditional_steps_8_to_11_gate': 'OPEN' if decision in ['COMPLETE_MIGRATION_SUPPORT', 'DIRECTIONALLY_CONSISTENT_PARTIAL_SUPPORT'] else 'CLOSED',
        'evidence_boundary': 'Observational transfer-state evidence only; not replication of client TYMP knockdown and not causal closure.',
    }
    (RESULTS / 'full_score_evidence_decision_v1.json').write_text(json.dumps(decision_record, indent=2), encoding='utf-8')
    return primary, meta, decision_record


def analyse_cptac():
    rng = np.random.default_rng(SEED + 10)
    score_frames = []
    associations = []
    components = []
    for dataset in ['CPTAC', 'CPTAC_val']:
        frame = prepare_unitedmet(dataset, CPTAC_FEATURES, 'S3')
        score_frames.append(frame.reset_index())
        rho, p = stats.spearmanr(frame['TYMP'], frame['S3_rank_for_spearman'])
        low, high, valid_boot = bootstrap_spearman(frame['TYMP'], frame['S3_rank_for_spearman'], rng)
        partial_rho, partial_p, _, _ = partial_rank_corr(frame['TYMP'], frame['S3_rank_for_spearman'], rank_z(frame['MTAP']))
        associations.append({
            'dataset': dataset, 'coverage': '3/7', 'n': len(frame), 'rho': float(rho),
            'bootstrap_ci_low': low, 'bootstrap_ci_high': high,
            'valid_bootstrap_draws': valid_boot, 'p': float(p),
            'partial_rho_controlling_MTAP': partial_rho, 'partial_p': partial_p,
        })
        cohort_components = []
        for feature in CPTAC_FEATURES:
            value_rho, value_p = stats.spearmanr(frame['TYMP'], frame[feature])
            expected = -FEATURES[feature][1]
            cohort_components.append({
                'dataset': dataset, 'coverage': '3/7', 'feature': feature, 'n': len(frame),
                'rho': float(value_rho), 'p': float(value_p), 'expected_rho_sign': int(expected),
                'expected_direction_hit': bool(np.sign(value_rho) == expected),
            })
        q = bh_adjust([row['p'] for row in cohort_components])
        for row, value in zip(cohort_components, q):
            row['q_within_cohort_3'] = float(value)
            components.append(row)
    primary_q = bh_adjust([row['p'] for row in associations])
    partial_q = bh_adjust([row['partial_p'] for row in associations])
    for row, q1, q2 in zip(associations, primary_q, partial_q):
        row['q_s3_primary_family_2'] = float(q1)
        row['partial_q_family_2'] = float(q2)
    save_tsv(pd.concat(score_frames, ignore_index=True), 'cptac_s3_sample_scores_v1.tsv')
    save_tsv(pd.DataFrame(associations), 'cptac_s3_associations_v1.tsv')
    save_tsv(pd.DataFrame(components), 'cptac_s3_component_associations_v1.tsv')
    return associations


def normalize_cell(value):
    text = str(value).upper().strip()
    if text.endswith('/ATCC'):
        text = text[:-5]
    return ''.join(char for char in text if char.isalnum())


def analyse_nci60():
    rng = np.random.default_rng(SEED + 20)
    raw = pd.read_csv(INPUTS / 'NCI60_WEB_DATA_METABOLON_ALL.TXT')
    selected = raw.loc[raw['TITLE'].isin(NCI_ORIGINAL.values())].copy()
    selected['cell_key'] = selected['cellname'].map(normalize_cell)
    selected['feature'] = selected['TITLE'].map({value: key for key, value in NCI_ORIGINAL.items()})
    selected['Value'] = pd.to_numeric(selected['Value'], errors='coerce')
    agg = selected.groupby(['cell_key', 'feature'], as_index=False)['Value'].mean()
    wide = agg.pivot(index='cell_key', columns='feature', values='Value')
    existing = pd.read_csv(INPUTS / 'nci60_all_lineages_preprocessed_v1.tsv', sep='\t')
    existing = existing.set_index('cell_key')
    frame = existing[['metabolomics_cell_line', 'tissue', 'is_renal', 'tymp_log2_fpkm_plus1', 'mtap_log2_fpkm_plus1']].join(wide[NCI_FEATURES])
    frame = frame.rename(columns={'tymp_log2_fpkm_plus1': 'TYMP', 'mtap_log2_fpkm_plus1': 'MTAP'}).dropna()
    scored = score_frame(frame, NCI_FEATURES, 'S5')
    scored.index.name = 'cell_key'
    rho, p = stats.spearmanr(scored['TYMP'], scored['S5_rank_for_spearman'])
    low, high, valid_boot = bootstrap_spearman(scored['TYMP'], scored['S5_rank_for_spearman'], rng)
    tissues = pd.get_dummies(scored['tissue'], drop_first=True, dtype=float).to_numpy()
    partial_rho, partial_p, xres, yres = partial_rank_corr(scored['TYMP'], scored['S5_rank_for_spearman'], tissues)
    perm_rho, perm_p = within_group_permutation(xres, yres, scored['tissue'].to_numpy(), rng)
    all_result = pd.DataFrame([{
        'dataset': 'NCI60', 'coverage': '5/7', 'n': len(scored), 'rho': float(rho),
        'bootstrap_ci_low': low, 'bootstrap_ci_high': high, 'valid_bootstrap_draws': valid_boot,
        'p': float(p), 'tissue_adjusted_partial_rho': partial_rho,
        'tissue_adjusted_analytic_p': partial_p, 'within_tissue_permutations': N_PERM,
        'within_tissue_permutation_rho': perm_rho, 'within_tissue_permutation_p': perm_p,
    }])
    renal = scored.loc[scored['is_renal'].astype(str).str.upper().eq('TRUE')].copy()
    renal_rho = float(stats.spearmanr(renal['TYMP'], renal['S5_rank_for_spearman']).statistic)
    loo = []
    for cell_key in renal.index:
        subset = renal.drop(index=cell_key)
        value = float(stats.spearmanr(subset['TYMP'], subset['S5_rank_for_spearman']).statistic)
        loo.append({'omitted_cell_key': cell_key, 'n': len(subset), 'rho': value, 'direction_negative': bool(value < 0)})
    renal_rows = renal.reset_index()[['cell_key', 'metabolomics_cell_line', 'TYMP', 'MTAP'] + NCI_FEATURES + ['S5', 'S5_rank_for_spearman']]
    renal_rows.insert(1, 'dataset', 'NCI60_renal')
    renal_summary = pd.DataFrame([{
        'dataset': 'NCI60_renal', 'coverage': '5/7', 'n': len(renal), 'rho_descriptive': renal_rho,
        'loo_rho_min': min(row['rho'] for row in loo), 'loo_rho_max': max(row['rho'] for row in loo),
        'loo_negative_count': sum(row['direction_negative'] for row in loo),
        'inferential_tests': 'NOT_CALCULATED_N_LT_10',
    }])
    save_tsv(scored.reset_index(), 'nci60_s5_sample_scores_v1.tsv')
    save_tsv(all_result, 'nci60_s5_all_lineage_association_v1.tsv')
    save_tsv(renal_rows, 'nci60_s5_renal_line_level_v1.tsv')
    save_tsv(renal_summary, 'nci60_s5_renal_summary_v1.tsv')
    save_tsv(pd.DataFrame(loo), 'nci60_s5_renal_leave_one_line_out_v1.tsv')
    return all_result.iloc[0].to_dict(), renal_summary.iloc[0].to_dict()


def make_figures(primary_rows, meta, cptac_rows, nci_all):
    primary = pd.DataFrame(primary_rows)
    scores = pd.read_csv(RESULTS / 'unitedmet_full_score_sample_scores_v1.tsv', sep='\t')
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), constrained_layout=True)
    colors = {'RC18': '#2563A6', 'RC20': '#D97706'}
    for ax, dataset in zip(axes, ['RC18', 'RC20']):
        sub = scores.loc[scores['dataset'].eq(dataset)]
        ax.scatter(sub['TYMP'], sub['S7'], s=38, alpha=0.82, color=colors[dataset], edgecolor='white', linewidth=0.5)
        slope, intercept = np.polyfit(sub['TYMP'], sub['S7'], 1)
        xline = np.linspace(sub['TYMP'].min(), sub['TYMP'].max(), 100)
        ax.plot(xline, intercept + slope * xline, color='#263238', linewidth=1.2)
        row = primary.loc[primary['dataset'].eq(dataset)].iloc[0]
        ax.set_title(f"{dataset}: n={int(row['n'])}, rho={row['rho']:.3f}, q={row['q_primary_family_2']:.3g}")
        ax.set_xlabel('TYMP log2(TPM+1)')
        ax.set_ylabel('Unweighted shTYMP-like S7')
        ax.axhline(0, color='#B0BEC5', linewidth=0.8)
    fig.suptitle('Complete-case seven-feature primary associations', fontsize=13, fontweight='bold')
    fig.savefig(FIGURES / 'FIG1_FULL_SCORE_PRIMARY_ASSOCIATIONS_v1.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES / 'FIG1_FULL_SCORE_PRIMARY_ASSOCIATIONS_v1.pdf', bbox_inches='tight')
    plt.close(fig)

    components = pd.read_csv(RESULTS / 'full_score_component_associations_v1.tsv', sep='\t')
    pivot = components.pivot(index='feature', columns='dataset', values='rho').reindex(FULL_FEATURES)
    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    image = ax.imshow(pivot.to_numpy(), cmap='RdBu_r', vmin=-0.6, vmax=0.6, aspect='auto')
    ax.set_xticks(range(len(pivot.columns)), pivot.columns)
    ax.set_yticks(range(len(pivot.index)), pivot.index)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            ax.text(j, i, f'{pivot.iloc[i, j]:.2f}', ha='center', va='center', fontsize=9)
    ax.set_title('TYMP correlations with seven score components')
    fig.colorbar(image, ax=ax, label='Spearman rho')
    fig.savefig(FIGURES / 'FIG2_COMPONENT_CORRELATIONS_v1.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES / 'FIG2_COMPONENT_CORRELATIONS_v1.pdf', bbox_inches='tight')
    plt.close(fig)

    loo = pd.read_csv(RESULTS / 'robustness_leave_one_feature_v1.tsv', sep='\t')
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True, sharey=True)
    for ax, dataset in zip(axes, ['RC18', 'RC20']):
        sub = loo.loc[loo['dataset'].eq(dataset)]
        ax.barh(sub['omitted_feature'], sub['rho'], color=colors[dataset], alpha=0.82)
        ax.axvline(0, color='#263238', linewidth=1)
        ax.set_title(dataset)
        ax.set_xlabel('TYMP vs six-feature score rho')
    fig.suptitle('Leave-one-feature robustness', fontsize=13, fontweight='bold')
    fig.savefig(FIGURES / 'FIG3_LEAVE_ONE_FEATURE_ROBUSTNESS_v1.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES / 'FIG3_LEAVE_ONE_FEATURE_ROBUSTNESS_v1.pdf', bbox_inches='tight')
    plt.close(fig)

    sensitivity = pd.DataFrame([
        {'dataset': row['dataset'], 'rho': row['rho'], 'n': row['n'], 'coverage': '3/7'} for row in cptac_rows
    ] + [{'dataset': 'NCI60', 'rho': nci_all['rho'], 'n': nci_all['n'], 'coverage': '5/7'}])
    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    ax.bar(sensitivity['dataset'], sensitivity['rho'], color=['#607D8B', '#90A4AE', '#7B1FA2'])
    ax.axhline(0, color='#263238', linewidth=1)
    for i, row in sensitivity.iterrows():
        ax.text(i, row['rho'] + (0.025 if row['rho'] >= 0 else -0.045), f"n={int(row['n'])}, {row['coverage']}", ha='center', fontsize=9)
    ax.set_ylabel('Spearman rho with TYMP')
    ax.set_title('Reduced-feature sensitivity analyses (not full-score validation)')
    fig.savefig(FIGURES / 'FIG4_REDUCED_FEATURE_SENSITIVITIES_v1.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIGURES / 'FIG4_REDUCED_FEATURE_SENSITIVITIES_v1.pdf', bbox_inches='tight')
    plt.close(fig)


def main():
    audit = json.loads((ROOT / 'qc_v1' / 'INPUT_AUDIT_STATUS_v1.json').read_text(encoding='utf-8'))
    if audit['status'] != 'PASS':
        raise RuntimeError('input audit did not pass')
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    if any(RESULTS.iterdir()) or any(FIGURES.iterdir()):
        raise RuntimeError('formal results/figures directory is not empty')

    primary, meta, decision = analyse_unitedmet_primary()
    cptac = analyse_cptac()
    nci_all, nci_renal = analyse_nci60()
    make_figures(primary, meta, cptac, nci_all)

    status = {
        'status': 'PASS_ANALYSIS_COMPLETE',
        'seed': SEED,
        'bootstrap_replicates': N_BOOT,
        'permutation_replicates': N_PERM,
        'primary_decision': decision['decision'],
        'conditional_steps_8_to_11_gate': decision['conditional_steps_8_to_11_gate'],
        'primary_results': primary,
        'meta_result': meta,
        'cptac_reduced_results': cptac,
        'nci60_all_lineage_result': nci_all,
        'nci60_renal_result': nci_renal,
        'formal_tsv_count': len(list(RESULTS.glob('*.tsv'))),
        'formal_figure_count': len(list(FIGURES.glob('*'))),
        'claims_not_authorized': [
            'public replication of client TYMP knockdown',
            'validated TYMP-MTA-SPP1 causal chain',
            'full seven-feature validation from 3/7 or 5/7 sensitivity scores',
        ],
    }
    (RESULTS / 'ANALYSIS_STATUS_v1.json').write_text(json.dumps(status, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(status, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
