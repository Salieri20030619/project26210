import csv
import hashlib
import json
import re
import shutil
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo


PROJECT = Path(r'D:\project26210_external_joint_v1')
ROOT = PROJECT / 'phase114_unweighted_seven_feature_analysis_local_v1'
DELIVERY = Path(r'C:\Users\31937\Downloads\项目26210_七特征分析结果_20260901_v1')
STATE_NAMES = ['PROJECT_STATUS.md', 'DECISIONS.md', 'RESULTS_MANIFEST.tsv', 'CODEX_HANDOFF.md']
NOW = datetime.now(ZoneInfo('Asia/Taipei')).isoformat(timespec='seconds')


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest().upper()


def append_once(path, marker, block):
    text = path.read_text(encoding='utf-8-sig')
    if marker not in text:
        path.write_text(text.rstrip() + '\n\n' + block.strip() + '\n', encoding='utf-8')


def make_delivery_manifest():
    mappings = [
        (ROOT / 'PHASE114_ANALYSIS_SUMMARY_v1.md', DELIVERY / 'PHASE114_ANALYSIS_SUMMARY_v1.md'),
        (ROOT / 'EVIDENCE_SUMMARY_v1.tsv', DELIVERY / 'EVIDENCE_SUMMARY_v1.tsv'),
    ]
    for source in sorted((ROOT / 'results_v1').iterdir()):
        if source.is_file():
            mappings.append((source, DELIVERY / 'tables' / source.name))
    for source in sorted((ROOT / 'figures_v1').iterdir()):
        if source.is_file():
            mappings.append((source, DELIVERY / 'figures' / source.name))
    for source in sorted((ROOT / 'validation_v1').iterdir()):
        if source.is_file():
            mappings.append((source, DELIVERY / 'validation' / source.name))
    rows = []
    for source, delivered in mappings:
        if not source.exists() or not delivered.exists():
            raise FileNotFoundError(f'{source} or {delivered}')
        source_hash = sha256(source)
        delivered_hash = sha256(delivered)
        if source_hash != delivered_hash:
            raise RuntimeError(f'delivery mismatch: {source.name}')
        rows.append([str(source), str(delivered), source_hash, source.stat().st_size, 'PASS_EXACT_COPY'])
    path = ROOT / 'DELIVERY_MANIFEST_v1.tsv'
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['authoritative_path', 'delivery_path', 'sha256', 'size_bytes', 'status'])
        writer.writerows(rows)
    return path, rows


def make_file_manifest(delivery_manifest):
    paths = []
    for path in ROOT.rglob('*'):
        if not path.is_file() or '__pycache__' in path.parts:
            continue
        if path.name == 'FILE_MANIFEST_SHA256_v1.tsv':
            continue
        paths.append(path)
    paths.append(delivery_manifest)
    paths = sorted(set(paths), key=lambda item: str(item).lower())
    manifest = ROOT / 'FILE_MANIFEST_SHA256_v1.tsv'
    with manifest.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['absolute_path', 'sha256', 'size_bytes', 'status'])
        for path in paths:
            writer.writerow([str(path), sha256(path), path.stat().st_size, 'PASS'])
    return manifest


def backup_state():
    backup = ROOT / 'state_backup_before_phase114_v1'
    index = 2
    while backup.exists():
        backup = ROOT / f'state_backup_before_phase114_v{index}'
        index += 1
    backup.mkdir()
    for name in STATE_NAMES:
        source = PROJECT / name
        if not source.exists():
            raise FileNotFoundError(str(source))
        shutil.copy2(source, backup / name)
    return backup


def update_markdown_state(backup):
    status_block = f'''## Phase 114 unweighted seven-feature public-data transfer analysis (2026-09-01)

- Status: **PASS_ANALYSIS_COMPLETE; PASS_INPUT_AUDIT; PASS_INDEPENDENT_VALIDATION 67/67; CONDITIONAL GATE CLOSED**.
- Scope: executed the frozen Phase112 steps 0-7 only. No RNA proxy, single-cell, spatial, survival or treatment-response projection was run.
- Inputs and biological sample sizes: UnitedMet RC18/RC20 contained 97/52 tumour specimens but only 26/13 complete seven-feature cases under the no-imputation rule. CPTAC/CPTAC_val had 50/55 complete 3/7 cases. NCI60 had 58 independent cell lines for 5/7 all-lineage analysis and eight renal lines for descriptive analysis.
- Methods: within-dataset average-rank z scores with fixed client-derived signs and equal weights; two-sided Spearman tests; 10,000 patient bootstraps; separate BH families; MTAP-adjusted partial rank correlations; Fisher-z fixed/random meta-analysis; leave-one-feature, leave-one-sample, 100,000-permutation and raw-z sensitivities. Mathematical score ties are preserved by an explicit stored average-rank column.
- Primary result: RC18 rho=-0.1371, 95% CI [-0.4810,0.2353], P=0.5042, q=0.9715; RC20 rho=0.0110, 95% CI [-0.5607,0.5643], P=q=0.9715. Fixed meta rho=-0.0926, 95% CI [-0.4087,0.2434], P=0.5938. Formal decision: `CURRENT_DATA_DO_NOT_SUPPORT_GENERAL_SEVEN_FEATURE_TRANSFER`.
- Secondary results: MTAP adjustment remained RC18 negative and RC20 near zero. RC18 had 3/7 component directions as expected; RC20 had 4/7, with no reproducible component pattern. RC18 leave-one-sample was 26/26 weakly negative; RC20 was mixed (7/13 negative). Raw-z sensitivities were near zero.
- Reduced sensitivities: CPTAC S3 rho=-0.4269, q=0.0040; CPTAC_val rho=-0.2455, q=0.0709. This is a reduced three-feature candidate, not full-score validation. NCI60 all-lineage S5 rho=0.0479, P=0.7209 and tissue-adjusted permutation P=0.3674. Renal n=8 rho=-0.3095 with 8/8 negative leave-one-line-out estimates; no P/FDR/CI was calculated.
- Evidence boundary: client TYMP-knockdown results and public observational evidence remain separate. The prior single-MTA null result remains unchanged. Downstream experimental plausibility and reduced scores cannot rescue failed full-score public validation.
- Authoritative paths: `{ROOT / 'PHASE114_ANALYSIS_SUMMARY_v1.md'}`, `{ROOT / 'EVIDENCE_SUMMARY_v1.tsv'}`, `{ROOT / 'results_v1'}`, `{ROOT / 'figures_v1'}`, and `{ROOT / 'validation_v1'}`. Downloads package: `{DELIVERY}`. Pre-update state backup: `{backup}`.
- Current breakpoint: do not start Phase112 steps 8-11. A lower-missingness external 7/7 cohort search or a newly preregistered three-feature submodule would be a new task. Active jobs: none.'''
    append_once(PROJECT / 'PROJECT_STATUS.md', '## Phase 114 unweighted seven-feature public-data transfer analysis', status_block)

    decisions_block = '''## 2026-09-01 - Phase 114 seven-feature analysis decisions

- Accept the complete-case sample sizes RC18 n=26 and RC20 n=13. Do not report the total tumour counts 97/52 as if every specimen entered the seven-feature statistic.
- Accept `CURRENT_DATA_DO_NOT_SUPPORT_GENERAL_SEVEN_FEATURE_TRANSFER` as the formal full-score decision because RC18 and RC20 directions differ and the meta effect is near zero.
- Close the preregistered conditional gate. Do not construct an RNA proxy or proceed to single-cell, spatial, prognosis or treatment projection under the current seven-feature plan.
- Preserve the CPTAC/CPTAC_val 3/7 negative directions as a reduced candidate submodule only. The discovery cohort passed FDR; the validation cohort did not. Do not relabel this as full seven-feature validation.
- Preserve NCI60 all-lineage 5/7 as unsupported. Preserve renal n=8 as a descriptive negative-direction candidate without P value, FDR or confidence interval.
- Preserve the prior direct TYMP-MTA null endpoint as a separate result; the composite score does not overwrite it.
- Retain all download, audit and output-precision failures as implementation/recovery audits. Only the final 67/67-validated output is formal.
- Store exact average ranks used for Spearman tests so mathematical ties remain reproducible after TSV round-trip.
- Any search for another 7/7 dataset or analysis of a three-feature submodule requires a new scope and preregistration.'''
    append_once(PROJECT / 'DECISIONS.md', '## 2026-09-01 - Phase 114 seven-feature analysis decisions', decisions_block)

    handoff_block = f'''## Phase 114 seven-feature analysis current breakpoint (2026-09-01)

- Completed: frozen steps 0-7, including input identity/mapping audit, RC18/RC20 full S7 primary analysis, component and MTAP-adjusted tests, fixed/random meta-analysis, all prespecified robustness checks, CPTAC 3/7 and NCI60 5/7 sensitivities, figures and independent recomputation.
- Validation: official Zenodo MD5 8/8; final input audit PASS; analysis status PASS; independent validation 67/67 PASS; downloaded package contains 33 exact-copy files.
- Formal decision: current data do not support general seven-feature transfer. RC18 rho=-0.1371 (n=26) and RC20 rho=0.0110 (n=13); fixed meta rho=-0.0926, P=0.5938.
- Exact resume point: stop under the existing plan. Do not run conditional Phase112 steps 8-11. Ask for a new decision only if the user wants a new lower-missingness 7/7 dataset search or a separately preregistered three-feature candidate analysis.
- Active jobs: none. The `sgs10-codex` Phase114 directory contains the verified download cache and retained segmented-download recovery files; local `D:\project26210_external_joint_v1` remains authoritative, and mirror/cache-specific paths are not directly executable on the other machine.
- Key paths: `{ROOT / 'PHASE114_ANALYSIS_SUMMARY_v1.md'}`, `{ROOT / 'results_v1' / 'ANALYSIS_STATUS_v1.json'}`, `{ROOT / 'validation_v1' / 'independent_validation_status_v1.json'}`, `{DELIVERY}`.
- Claims that must not be made: public replication of client TYMP knockdown; validated full seven-feature transfer; validated TYMP-MTA-SPP1 chain; full-score rescue by CPTAC 3/7; inferential renal NCI60 conclusion; any downstream clinical predictor from this phase.'''
    append_once(PROJECT / 'CODEX_HANDOFF.md', '## Phase 114 seven-feature analysis current breakpoint', handoff_block)


def artifact_type(path):
    name = path.name.lower()
    if name.endswith('.png'):
        return 'figure_png'
    if name.endswith('.pdf'):
        return 'figure_pdf'
    if name.endswith('.json'):
        return 'validation_json' if 'validation' in name or 'audit' in name else 'result_json'
    if name.endswith('.tsv'):
        return 'result_table_tsv'
    if name.endswith('.py') or name.endswith('.sh'):
        return 'formal_script'
    if 'failure_audit' in name:
        return 'failure_audit_md'
    if 'summary' in name:
        return 'analysis_summary_md'
    if 'preregistration' in name:
        return 'preregistration_md'
    return 'formal_artifact'


def result_id_for(path):
    relative = path.relative_to(ROOT)
    stem = re.sub(r'[^a-z0-9]+', '_', str(relative).lower()).strip('_')
    return f'phase114_{stem}'


def update_results_manifest(file_manifest, delivery_manifest):
    paths = [
        ROOT / 'PHASE114_PREREGISTRATION_v1.md',
        ROOT / 'PHASE114_ANALYSIS_SUMMARY_v1.md',
        ROOT / 'EVIDENCE_SUMMARY_v1.tsv',
        ROOT / 'DOWNLOAD_AND_IMPLEMENTATION_FAILURE_AUDIT_v1.md',
        ROOT / 'inputs_v1' / 'ANALYSIS_PARAMETERS_v1.json',
    ]
    paths.extend(sorted((ROOT / 'scripts_v1').glob('*.*')))
    paths.extend(sorted((ROOT / 'qc_v1').glob('*.*')))
    paths.extend(sorted((ROOT / 'results_v1').glob('*.*')))
    paths.extend(sorted((ROOT / 'figures_v1').glob('*.*')))
    paths.extend(sorted((ROOT / 'validation_v1').glob('*.*')))
    paths.extend([file_manifest, delivery_manifest])
    paths = [path for path in dict.fromkeys(paths) if path.is_file() and '__pycache__' not in path.parts]

    manifest = PROJECT / 'RESULTS_MANIFEST.tsv'
    with manifest.open('r', encoding='utf-8-sig', newline='') as handle:
        existing = list(csv.DictReader(handle, delimiter='\t'))
    existing_ids = {row['result_id'] for row in existing}
    rows = []
    for path in paths:
        result_id = result_id_for(path)
        if result_id in existing_ids:
            continue
        evidence = f'SHA256 {sha256(path)}; input audit PASS; analysis complete; independent validation 67/67 PASS; downstream gate CLOSED'
        notes = 'Phase114 formal seven-feature analysis artifact; no conditional RNA/single-cell/spatial/clinical analysis executed'
        rows.append([result_id, 'phase114_unweighted_seven_feature_analysis', artifact_type(path), str(path), NOW, 'PASS', evidence, notes])
    with manifest.open('a', encoding='utf-8', newline='') as handle:
        csv.writer(handle, delimiter='\t', lineterminator='\n').writerows(rows)
    return len(rows)


def validate_state(delivery_rows, new_manifest_rows):
    checks = {}
    for name in STATE_NAMES:
        checks[f'exists_{name}'] = (PROJECT / name).exists()
    with (PROJECT / 'RESULTS_MANIFEST.tsv').open('r', encoding='utf-8-sig', newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    ids = [row['result_id'] for row in rows]
    checks['manifest_schema'] = list(rows[0]) == [
        'result_id', 'stage', 'artifact_type', 'absolute_path', 'generated_at',
        'validation_status', 'validation_evidence', 'notes'
    ]
    checks['manifest_unique_ids'] = len(ids) == len(set(ids))
    phase_rows = [row for row in rows if row['stage'] == 'phase114_unweighted_seven_feature_analysis']
    checks['phase114_rows_registered'] = len(phase_rows) == new_manifest_rows
    checks['phase114_paths_exist'] = all(Path(row['absolute_path']).exists() for row in phase_rows)
    checks['completion_markers'] = all(
        marker in (PROJECT / file_name).read_text(encoding='utf-8-sig')
        for file_name, marker in [
            ('PROJECT_STATUS.md', 'Phase 114 unweighted seven-feature public-data transfer analysis'),
            ('DECISIONS.md', 'Phase 114 seven-feature analysis decisions'),
            ('CODEX_HANDOFF.md', 'Phase 114 seven-feature analysis current breakpoint'),
        ]
    )
    checks['delivery_33_exact_copies'] = len(delivery_rows) == 33 and all(row[-1] == 'PASS_EXACT_COPY' for row in delivery_rows)
    checks['input_audit_pass'] = json.loads((ROOT / 'qc_v1' / 'INPUT_AUDIT_STATUS_v1.json').read_text(encoding='utf-8'))['status'] == 'PASS'
    checks['analysis_status_pass'] = json.loads((ROOT / 'results_v1' / 'ANALYSIS_STATUS_v1.json').read_text(encoding='utf-8'))['status'] == 'PASS_ANALYSIS_COMPLETE'
    checks['independent_validation_67_67'] = json.loads((ROOT / 'validation_v1' / 'independent_validation_status_v1.json').read_text(encoding='utf-8')) == {
        'status': 'PASS', 'checks_passed': 67, 'checks_total': 67,
        'primary_decision': 'CURRENT_DATA_DO_NOT_SUPPORT_GENERAL_SEVEN_FEATURE_TRANSFER',
        'conditional_gate': 'CLOSED'
    }
    result = {
        'status': 'PASS' if all(checks.values()) else 'FAIL',
        'generated_at': NOW,
        'checks': checks,
        'phase114_manifest_rows': len(phase_rows),
        'delivery_directory': str(DELIVERY),
        'authoritative_root': str(PROJECT),
        'active_jobs': 'none',
    }
    out = ROOT / 'validation_v1' / 'PHASE114_STATE_UPDATE_VALIDATION_v1.json'
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    if result['status'] != 'PASS':
        raise RuntimeError(json.dumps(result, ensure_ascii=False))
    manifest = PROJECT / 'RESULTS_MANIFEST.tsv'
    result_id = 'phase114_state_update_validation_v1'
    with manifest.open('r', encoding='utf-8-sig', newline='') as handle:
        existing_ids = {row['result_id'] for row in csv.DictReader(handle, delimiter='\t')}
    if result_id not in existing_ids:
        with manifest.open('a', encoding='utf-8', newline='') as handle:
            csv.writer(handle, delimiter='\t', lineterminator='\n').writerow([
                result_id, 'phase114_project_state_update', 'state_update_validation_json', str(out), NOW, 'PASS',
                'Phase114 mandatory four-file completion gate PASS; local delivery 33/33 exact copies; validation 67/67',
                'State validation generated before self-registration to avoid self-hash recursion',
            ])
    return out, result


def main():
    delivery_manifest, delivery_rows = make_delivery_manifest()
    file_manifest = make_file_manifest(delivery_manifest)
    backup = backup_state()
    update_markdown_state(backup)
    new_rows = update_results_manifest(file_manifest, delivery_manifest)
    state_validation, result = validate_state(delivery_rows, new_rows)
    print(json.dumps({'state_validation': str(state_validation), **result}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
