import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from scipy import stats


ROOT = Path(r'D:\project26210_external_joint_v1\phase114_unweighted_seven_feature_analysis_local_v1')
RESULTS = ROOT / 'results_v1'
FIGURES = ROOT / 'figures_v1'
VALIDATION = ROOT / 'validation_v1'


def rank_z(values):
    ranks = stats.rankdata(np.asarray(values, dtype=float), method='average')
    return (ranks - ranks.mean()) / ranks.std(ddof=0)


def bh(values):
    p = np.asarray(values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def close(a, b, tolerance=1e-8):
    return bool(np.isfinite(a) and np.isfinite(b) and abs(float(a) - float(b)) <= tolerance)


def add(checks, check_id, passed, evidence):
    checks.append({'check_id': check_id, 'pass': bool(passed), 'evidence': str(evidence)})


def main():
    VALIDATION.mkdir(parents=True, exist_ok=True)
    checks = []
    audit = json.loads((ROOT / 'qc_v1' / 'INPUT_AUDIT_STATUS_v1.json').read_text(encoding='utf-8'))
    status = json.loads((RESULTS / 'ANALYSIS_STATUS_v1.json').read_text(encoding='utf-8'))
    decision = json.loads((RESULTS / 'full_score_evidence_decision_v1.json').read_text(encoding='utf-8'))
    add(checks, 'input_audit_pass', audit['status'] == 'PASS', audit['status'])
    add(checks, 'analysis_status_pass', status['status'] == 'PASS_ANALYSIS_COMPLETE', status['status'])

    scores = pd.read_csv(RESULTS / 'unitedmet_full_score_sample_scores_v1.tsv', sep='\t')
    primary = pd.read_csv(RESULTS / 'full_score_primary_associations_v1.tsv', sep='\t')
    for dataset, expected_n in [('RC18', 26), ('RC20', 13)]:
        sub = scores.loc[scores['dataset'].eq(dataset)].copy()
        add(checks, f'{dataset}_complete_case_n', len(sub) == expected_n, len(sub))
        signed = [column for column in sub if column.startswith('signed_zrank__')]
        reconstructed = sub[signed].mean(axis=1)
        add(checks, f'{dataset}_score_formula', np.max(np.abs(reconstructed - sub['S7'])) < 1e-9, np.max(np.abs(reconstructed - sub['S7'])))
        row = primary.loc[primary['dataset'].eq(dataset)].iloc[0]
        rho, p = stats.spearmanr(sub['TYMP'], sub['S7_rank_for_spearman'])
        add(checks, f'{dataset}_primary_rho', close(rho, row['rho']), f'{rho} vs {row.rho}')
        add(checks, f'{dataset}_primary_p', close(p, row['p']), f'{p} vs {row.p}')
        add(checks, f'{dataset}_bootstrap_order', row['bootstrap_ci_low'] <= row['rho'] <= row['bootstrap_ci_high'], f"{row['bootstrap_ci_low']},{row['rho']},{row['bootstrap_ci_high']}")

    recomputed_q = bh(primary['p'])
    add(checks, 'primary_bh_recomputed', np.max(np.abs(recomputed_q - primary['q_primary_family_2'])) < 1e-9, recomputed_q)
    add(checks, 'primary_directions_opposite', list(np.sign(primary['rho'])) == [-1.0, 1.0], list(primary['rho']))

    rho = np.clip(primary['rho'].to_numpy(), -0.999999, 0.999999)
    n = primary['n'].to_numpy(dtype=float)
    z = np.arctanh(rho)
    weight = n - 3.0
    fixed_z = float(np.sum(weight * z) / np.sum(weight))
    fixed_rho = math.tanh(fixed_z)
    fixed_se = math.sqrt(1.0 / float(np.sum(weight)))
    fixed_p = float(2 * stats.norm.sf(abs(fixed_z / fixed_se)))
    meta = pd.read_csv(RESULTS / 'full_score_meta_analysis_v1.tsv', sep='\t').iloc[0]
    add(checks, 'fixed_meta_rho_recomputed', close(fixed_rho, meta['fixed_rho']), f'{fixed_rho} vs {meta.fixed_rho}')
    add(checks, 'fixed_meta_p_recomputed', close(fixed_p, meta['fixed_p']), f'{fixed_p} vs {meta.fixed_p}')
    add(checks, 'formal_decision_negative', decision['decision'] == 'CURRENT_DATA_DO_NOT_SUPPORT_GENERAL_SEVEN_FEATURE_TRANSFER', decision['decision'])
    add(checks, 'conditional_gate_closed', decision['conditional_steps_8_to_11_gate'] == 'CLOSED', decision['conditional_steps_8_to_11_gate'])

    components = pd.read_csv(RESULTS / 'full_score_component_associations_v1.tsv', sep='\t')
    add(checks, 'component_rows_14', len(components) == 14, len(components))
    for dataset in ['RC18', 'RC20']:
        subscore = scores.loc[scores['dataset'].eq(dataset)].set_index('sample_id')
        subcomp = components.loc[components['dataset'].eq(dataset)]
        values = []
        for _, row in subcomp.iterrows():
            value_rho, value_p = stats.spearmanr(subscore['TYMP'], subscore[row['feature']])
            values.append(value_p)
            add(checks, f"{dataset}_{row['feature']}_component_rho", close(value_rho, row['rho']), value_rho)
        add(checks, f'{dataset}_component_bh', np.max(np.abs(bh(values) - subcomp['q_within_cohort_7'])) < 1e-8, max(abs(bh(values) - subcomp['q_within_cohort_7'])))

    partials = pd.read_csv(RESULTS / 'full_score_mtap_partial_v1.tsv', sep='\t')
    add(checks, 'partial_rows_2', len(partials) == 2, len(partials))
    add(checks, 'partial_bh_recomputed', np.max(np.abs(bh(partials['p']) - partials['q_partial_family_2'])) < 1e-9, bh(partials['p']))

    feature_loo = pd.read_csv(RESULTS / 'robustness_leave_one_feature_v1.tsv', sep='\t')
    sample_loo = pd.read_csv(RESULTS / 'robustness_leave_one_sample_v1.tsv', sep='\t')
    add(checks, 'leave_one_feature_rows_14', len(feature_loo) == 14, len(feature_loo))
    add(checks, 'leave_one_sample_rows_39', len(sample_loo) == 39, len(sample_loo))
    add(checks, 'RC18_all_sample_loo_negative', bool(sample_loo.loc[sample_loo['dataset'].eq('RC18'), 'direction_negative'].all()), sample_loo.loc[sample_loo['dataset'].eq('RC18'), 'rho'].min())
    add(checks, 'RC20_sample_loo_mixed', 0 < int(sample_loo.loc[sample_loo['dataset'].eq('RC20'), 'direction_negative'].sum()) < 13, int(sample_loo.loc[sample_loo['dataset'].eq('RC20'), 'direction_negative'].sum()))

    permutation = pd.read_csv(RESULTS / 'robustness_permutation_v1.tsv', sep='\t')
    add(checks, 'permutation_rows_2', len(permutation) == 2, len(permutation))
    add(checks, 'permutation_count_100000', bool((permutation['permutations'] == 100000).all()), list(permutation['permutations']))
    add(checks, 'permutation_null_concordant', bool((permutation['permutation_p'] > 0.45).all()), list(permutation['permutation_p']))

    raw_sensitivity = pd.read_csv(RESULTS / 'robustness_raw_z_sensitivity_v1.tsv', sep='\t')
    add(checks, 'raw_z_near_null_both', bool((raw_sensitivity['rho_raw_z_score'].abs() < 0.06).all()), list(raw_sensitivity['rho_raw_z_score']))
    contribution = pd.read_csv(RESULTS / 'feature_contribution_audit_v1.tsv', sep='\t')
    for dataset in ['RC18', 'RC20']:
        total = contribution.loc[contribution['dataset'].eq(dataset), 'variance_contribution_fraction'].sum()
        add(checks, f'{dataset}_contributions_sum_one', close(total, 1.0), total)

    cptac_scores = pd.read_csv(RESULTS / 'cptac_s3_sample_scores_v1.tsv', sep='\t')
    cptac = pd.read_csv(RESULTS / 'cptac_s3_associations_v1.tsv', sep='\t')
    for dataset, expected_n in [('CPTAC', 50), ('CPTAC_val', 55)]:
        sub = cptac_scores.loc[cptac_scores['dataset'].eq(dataset)]
        row = cptac.loc[cptac['dataset'].eq(dataset)].iloc[0]
        rho, p = stats.spearmanr(sub['TYMP'], sub['S3_rank_for_spearman'])
        add(checks, f'{dataset}_S3_n', len(sub) == expected_n, len(sub))
        add(checks, f'{dataset}_S3_rho', close(rho, row['rho']), rho)
        add(checks, f'{dataset}_S3_p', close(p, row['p']), p)
    add(checks, 'CPTAC_S3_primary_bh', np.max(np.abs(bh(cptac['p']) - cptac['q_s3_primary_family_2'])) < 1e-9, bh(cptac['p']))
    add(checks, 'S3_not_full_score_label', bool((cptac['coverage'] == '3/7').all()), list(cptac['coverage']))

    nci_scores = pd.read_csv(RESULTS / 'nci60_s5_sample_scores_v1.tsv', sep='\t')
    nci_all = pd.read_csv(RESULTS / 'nci60_s5_all_lineage_association_v1.tsv', sep='\t').iloc[0]
    nci_rho, nci_p = stats.spearmanr(nci_scores['TYMP'], nci_scores['S5_rank_for_spearman'])
    add(checks, 'NCI60_S5_n58', len(nci_scores) == 58, len(nci_scores))
    add(checks, 'NCI60_S5_rho', close(nci_rho, nci_all['rho']), nci_rho)
    add(checks, 'NCI60_S5_p', close(nci_p, nci_all['p']), nci_p)
    add(checks, 'NCI60_within_tissue_permutation_100000', int(nci_all['within_tissue_permutations']) == 100000, nci_all['within_tissue_permutations'])
    renal = pd.read_csv(RESULTS / 'nci60_s5_renal_line_level_v1.tsv', sep='\t')
    renal_summary = pd.read_csv(RESULTS / 'nci60_s5_renal_summary_v1.tsv', sep='\t').iloc[0]
    renal_rho = stats.spearmanr(renal['TYMP'], renal['S5_rank_for_spearman']).statistic
    add(checks, 'NCI60_renal_n8', len(renal) == 8, len(renal))
    add(checks, 'NCI60_renal_rho', close(renal_rho, renal_summary['rho_descriptive']), renal_rho)
    add(checks, 'NCI60_renal_no_inference', renal_summary['inferential_tests'] == 'NOT_CALCULATED_N_LT_10', renal_summary['inferential_tests'])
    add(checks, 'NCI60_7860_once', int((renal['cell_key'].astype(str) == '7860').sum()) == 1, list(renal['cell_key']))
    renal_loo = pd.read_csv(RESULTS / 'nci60_s5_renal_leave_one_line_out_v1.tsv', sep='\t')
    add(checks, 'NCI60_renal_loo_8_negative', len(renal_loo) == 8 and bool(renal_loo['direction_negative'].all()), list(renal_loo['rho']))

    pngs = sorted(FIGURES.glob('*.png'))
    pdfs = sorted(FIGURES.glob('*.pdf'))
    add(checks, 'four_png_figures', len(pngs) == 4, len(pngs))
    add(checks, 'four_pdf_figures', len(pdfs) == 4, len(pdfs))
    dimensions_ok = True
    for path in pngs:
        with Image.open(path) as image:
            dimensions_ok &= image.width >= 1800 and image.height >= 1200
    add(checks, 'figure_dimensions', dimensions_ok, [(p.name, Image.open(p).size) for p in pngs])

    checks_frame = pd.DataFrame(checks)
    checks_frame.to_csv(VALIDATION / 'independent_validation_checks_v1.tsv', sep='\t', index=False)
    passed = int(checks_frame['pass'].sum())
    result = {
        'status': 'PASS' if passed == len(checks) else 'FAIL',
        'checks_passed': passed,
        'checks_total': len(checks),
        'primary_decision': decision['decision'],
        'conditional_gate': decision['conditional_steps_8_to_11_gate'],
    }
    (VALIDATION / 'independent_validation_status_v1.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result['status'] == 'PASS' else 1)


if __name__ == '__main__':
    main()
