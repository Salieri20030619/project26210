import csv
import hashlib
import json
import shutil
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo


PROJECT = Path(r'D:\project26210_external_joint_v1')
ROOT = PROJECT / 'phase113_client_concise_analysis_plan_local_v1'
OUTPUT = ROOT / 'output_v1'
VALIDATION = ROOT / 'validation_v1'
DELIVERY = Path(r'C:\Users\31937\Downloads\项目26210_七特征分析思路_概括版_20260901_v1')
DOCX_NAME = '项目26210_无权重七特征公共数据迁移分析思路_概括版_v1.docx'
PDF_NAME = '项目26210_无权重七特征公共数据迁移分析思路_概括版_v1.pdf'
STATE_FILES = ['PROJECT_STATUS.md', 'DECISIONS.md', 'RESULTS_MANIFEST.tsv', 'CODEX_HANDOFF.md']
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
    path = ROOT / 'DELIVERY_MANIFEST_v1.tsv'
    rows = []
    for name in [DOCX_NAME, PDF_NAME]:
        source = OUTPUT / name
        delivered = DELIVERY / name
        if not source.exists() or not delivered.exists():
            raise FileNotFoundError(name)
        source_hash = sha256(source)
        delivered_hash = sha256(delivered)
        if source_hash != delivered_hash:
            raise RuntimeError(f'delivery hash mismatch: {name}')
        rows.append([name, str(source), str(delivered), source_hash, 'PASS_EXACT_COPY'])
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['file_name', 'authoritative_path', 'delivery_path', 'sha256', 'status'])
        writer.writerows(rows)
    return path


def make_file_manifest(delivery_manifest):
    path = ROOT / 'FILE_MANIFEST_SHA256_v1.tsv'
    artifacts = [
        ROOT / 'PCa_reference.docx',
        ROOT / 'artifact.md',
        ROOT / 'scripts' / '01_build_concise_plan_v1.py',
        ROOT / 'scripts' / '02_validate_concise_plan_v1.py',
        ROOT / 'scripts' / '03_complete_phase113_state_v1.py',
        OUTPUT / DOCX_NAME,
        OUTPUT / PDF_NAME,
        VALIDATION / 'validation_summary_v1.json',
        VALIDATION / 'a11y_audit_v1.json',
        OUTPUT / 'render_final' / 'render_status.json',
        delivery_manifest,
        DELIVERY / DOCX_NAME,
        DELIVERY / PDF_NAME,
    ]
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['absolute_path', 'sha256', 'size_bytes', 'status'])
        for artifact in artifacts:
            if not artifact.exists():
                raise FileNotFoundError(str(artifact))
            writer.writerow([str(artifact), sha256(artifact), artifact.stat().st_size, 'PASS'])
    return path


def backup_state():
    backup = ROOT / 'state_backup_before_phase113_v1'
    if backup.exists():
        suffix = 2
        while (ROOT / f'state_backup_before_phase113_v{suffix}').exists():
            suffix += 1
        backup = ROOT / f'state_backup_before_phase113_v{suffix}'
    backup.mkdir(parents=True)
    for name in STATE_FILES:
        source = PROJECT / name
        if not source.exists():
            raise FileNotFoundError(str(source))
        shutil.copy2(source, backup / name)
    return backup


def update_state(file_manifest, delivery_manifest, backup):
    docx_hash = sha256(OUTPUT / DOCX_NAME)
    pdf_hash = sha256(OUTPUT / PDF_NAME)
    validation_hash = sha256(VALIDATION / 'validation_summary_v1.json')
    a11y_hash = sha256(VALIDATION / 'a11y_audit_v1.json')
    render_hash = sha256(OUTPUT / 'render_final' / 'render_status.json')
    reference_hash = sha256(ROOT / 'PCa_reference.docx')
    artifact_hash = sha256(ROOT / 'artifact.md')
    builder_hash = sha256(ROOT / 'scripts' / '01_build_concise_plan_v1.py')
    validator_hash = sha256(ROOT / 'scripts' / '02_validate_concise_plan_v1.py')
    completion_hash = sha256(ROOT / 'scripts' / '03_complete_phase113_state_v1.py')
    delivery_manifest_hash = sha256(delivery_manifest)
    file_manifest_hash = sha256(file_manifest)

    status_block = f'''## Phase 113 client concise seven-feature analysis plan (2026-09-01)

- Status: **PASS_CLIENT_CONCISE_PLAN_COMPLETE; PASS_WORD_RENDER_3/3; PASS_FORMAL_VALIDATION; PASS_A11Y_0/0/0**. This phase created a separate customer-facing concise plan and did not run any biological analysis.
- Scope and independent biological sample size: planning/documentation only; independent biological sample size is not applicable (n=0 analysed samples). The Phase112 16-page detailed execution plan remains unchanged and authoritative for implementation detail.
- Input/design authority: client-provided `PCa(1).docx` was retained as the concise-layout reference (SHA256 {reference_hash}); the user explicitly overrode its original fonts by requiring Microsoft YaHei 12 pt.
- Methods/content: the 3-page plan freezes the unweighted seven-feature shTYMP-like score, then summarizes feature mapping/QC, RC18/RC20 full-score validation, CPTAC/NCI-60 reduced sensitivities, gated RNA transfer, SPP1+ macrophage projection, spatial structure and prognosis/treatment analysis, with explicit stop rules.
- Principal positive deliverable: every visible Word run passed Microsoft YaHei 12 pt validation; A4 geometry matches the reference margins; PDF contains exactly 3 non-sparse pages; pixel boundary audit found no clipping risk; accessibility audit reported high=0, medium=0, low=0.
- Principal negative/boundary: no seven-feature score, association, survival result or public replication claim was computed. Public data cannot replace the client perturbation experiment; reduced 3/7 and 5/7 scores cannot validate the complete seven-feature score; previous single-MTA null evidence remains separate.
- Authoritative deliverables: `{OUTPUT / DOCX_NAME}` and `{OUTPUT / PDF_NAME}`. Downloads delivery: `{DELIVERY}`. Exact-copy hashes are recorded in `{delivery_manifest}`.
- Validation and recovery: `{VALIDATION / 'validation_summary_v1.json'}`, `{VALIDATION / 'a11y_audit_v1.json'}`, `{OUTPUT / 'render_final' / 'render_status.json'}`, and `{file_manifest}`. Pre-update state backup: `{backup}`.
- Current breakpoint: obtain client review of the concise plan. If analysis is authorized, resume from Phase112 steps 0-7 only; conditional RNA/single-cell/spatial/clinical projections remain gated by upstream results. Active background jobs: none.'''
    append_once(PROJECT / 'PROJECT_STATUS.md', '## Phase 113 client concise seven-feature analysis plan', status_block)

    decisions_block = '''## 2026-09-01 - Phase 113 concise client-plan decision

- Keep the Phase112 16-page plan unchanged as the detailed internal/execution document; do not supersede or delete it.
- Add a separate 3-page customer-facing conceptual plan following the information density and paragraph-led structure of the provided `PCa(1).docx`.
- Apply the user's explicit typography rule to all visible text: Microsoft YaHei 12 pt, including title and section labels.
- Retain all seven analysis steps and the fixed seven-feature definition, but omit command-level parameters, field-level mappings and full audit detail from the customer-facing document.
- Preserve the evidence boundary: client perturbation evidence and public database evidence are separate; no downstream experimental plausibility can rescue failed public validation.
- Preserve the prior single-MTA null endpoint. RC18/RC20 are full 7/7 primary cohorts; CPTAC 3/7 and NCI-60 5/7 remain reduced sensitivity analyses only.
- Require upstream gates before RNA, single-cell, spatial or clinical projection. Stopping after an unsupported gate is a valid completed result.
- No biological analysis is authorized or reported by this planning deliverable.'''
    append_once(PROJECT / 'DECISIONS.md', '## 2026-09-01 - Phase 113 concise client-plan decision', decisions_block)

    handoff_block = f'''## Phase 113 concise seven-feature plan current breakpoint (2026-09-01)

- Completed: separate 3-page customer-facing DOCX/PDF plan based on the concise structure of `PCa(1).docx`; the Phase112 16-page detailed plan remains intact.
- Validation: Word render 3/3 pages; every visible run Microsoft YaHei 12 pt; A4/reference margins PASS; formal content/font/page/pixel audit PASS; accessibility high=0, medium=0, low=0; downloaded copies hash-match authoritative files.
- Exact resume point: obtain client review. If the client approves execution, start Phase112 steps 0-7 only with frozen feature membership/signs and RC18/RC20 as full-score primary cohorts.
- Active jobs: none. No seven-feature analysis result exists yet.
- Key paths: `{OUTPUT / DOCX_NAME}`, `{OUTPUT / PDF_NAME}`, `{VALIDATION / 'validation_summary_v1.json'}`, `{DELIVERY}`.
- Claims that must not be made: public replication of client TYMP knockdown; validated TYMP-MTA-SPP1 chain; full seven-feature validation from CPTAC 3/7 or NCI-60 5/7; any result from the client's two unprovided figures; any survival or treatment conclusion before gated execution.'''
    append_once(PROJECT / 'CODEX_HANDOFF.md', '## Phase 113 concise seven-feature plan current breakpoint', handoff_block)

    manifest = PROJECT / 'RESULTS_MANIFEST.tsv'
    with manifest.open('r', encoding='utf-8-sig', newline='') as handle:
        existing = list(csv.DictReader(handle, delimiter='\t'))
    existing_ids = {row['result_id'] for row in existing}
    common_evidence = 'Word render 3/3; formal validation PASS; every visible run Microsoft YaHei 12 pt; a11y 0/0/0; no biological analysis executed'
    rows = [
        ('phase113_reference_pca_docx_v1', 'reference_docx', ROOT / 'PCa_reference.docx', 'PASS_IDENTITY_ONLY', f'SHA256 {reference_hash}; reference rendered 3/3', 'Client-provided design authority retained without modification'),
        ('phase113_template_distillation_v1', 'template_distillation_md', ROOT / 'artifact.md', 'PASS', f'SHA256 {artifact_hash}; {common_evidence}', 'Reference page, style and content-slot contract'),
        ('phase113_concise_plan_builder_v1', 'formal_script', ROOT / 'scripts' / '01_build_concise_plan_v1.py', 'PASS', f'SHA256 {builder_hash}; {common_evidence}', 'Builds the separate concise client plan'),
        ('phase113_concise_plan_validator_v1', 'formal_validation_script', ROOT / 'scripts' / '02_validate_concise_plan_v1.py', 'PASS', f'SHA256 {validator_hash}; {common_evidence}', 'Validates font, content, geometry, page count and pixel boundaries'),
        ('phase113_state_completion_script_v1', 'formal_script', ROOT / 'scripts' / '03_complete_phase113_state_v1.py', 'PASS', f'SHA256 {completion_hash}; state completion gate', 'Creates manifests, backups and four-file state update'),
        ('phase113_client_concise_plan_docx_v1', 'client_analysis_plan_docx', OUTPUT / DOCX_NAME, 'PASS', f'SHA256 {docx_hash}; {common_evidence}', 'Three-page customer-facing conceptual analysis plan'),
        ('phase113_client_concise_plan_pdf_v1', 'rendered_plan_pdf', OUTPUT / PDF_NAME, 'PASS', f'SHA256 {pdf_hash}; {common_evidence}', 'Final Word-rendered PDF'),
        ('phase113_concise_plan_validation_v1', 'validation_json', VALIDATION / 'validation_summary_v1.json', 'PASS', f'SHA256 {validation_hash}; 3 pages; font/content/geometry/pixel checks PASS', 'Formal delivery validation'),
        ('phase113_concise_plan_a11y_v1', 'accessibility_json', VALIDATION / 'a11y_audit_v1.json', 'PASS', f'SHA256 {a11y_hash}; a11y high=0 medium=0 low=0', 'Accessibility audit'),
        ('phase113_concise_plan_render_status_v1', 'render_status_json', OUTPUT / 'render_final' / 'render_status.json', 'PASS', f'SHA256 {render_hash}; Word render 3/3', 'Microsoft Word render status'),
        ('phase113_delivery_manifest_v1', 'delivery_manifest_tsv', delivery_manifest, 'PASS', f'SHA256 {delivery_manifest_hash}; authoritative/download hashes identical', 'Non-overwriting Downloads delivery manifest'),
        ('phase113_file_manifest_v1', 'file_manifest_tsv', file_manifest, 'PASS', f'SHA256 {file_manifest_hash}; all listed paths present and hashed', 'Phase113 recovery manifest'),
        ('phase113_download_docx_v1', 'delivered_docx', DELIVERY / DOCX_NAME, 'PASS', f'SHA256 {docx_hash}; exact copy of authoritative DOCX', 'User Downloads delivery copy'),
        ('phase113_download_pdf_v1', 'delivered_pdf', DELIVERY / PDF_NAME, 'PASS', f'SHA256 {pdf_hash}; exact copy of authoritative PDF', 'User Downloads delivery copy'),
    ]
    with manifest.open('a', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        for result_id, artifact_type, path, status, evidence, notes in rows:
            if result_id not in existing_ids:
                writer.writerow([result_id, 'phase113_client_concise_plan', artifact_type, str(path), NOW, status, evidence, notes])


def validate_state():
    checks = {}
    for name in STATE_FILES:
        checks[f'exists_{name}'] = (PROJECT / name).exists()
    with (PROJECT / 'RESULTS_MANIFEST.tsv').open('r', encoding='utf-8-sig', newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    ids = [row['result_id'] for row in rows]
    checks['manifest_schema_8_columns'] = list(rows[0].keys()) == [
        'result_id', 'stage', 'artifact_type', 'absolute_path', 'generated_at',
        'validation_status', 'validation_evidence', 'notes'
    ]
    checks['manifest_unique_ids'] = len(ids) == len(set(ids))
    phase_rows = [row for row in rows if row['stage'] == 'phase113_client_concise_plan']
    checks['phase113_rows_14'] = len(phase_rows) == 14
    checks['phase113_paths_exist'] = all(Path(row['absolute_path']).exists() for row in phase_rows)
    checks['completion_markers'] = all(
        marker in (PROJECT / file_name).read_text(encoding='utf-8-sig')
        for file_name, marker in [
            ('PROJECT_STATUS.md', 'Phase 113 client concise seven-feature analysis plan'),
            ('DECISIONS.md', 'Phase 113 concise client-plan decision'),
            ('CODEX_HANDOFF.md', 'Phase 113 concise seven-feature plan current breakpoint'),
        ]
    )
    checks['delivery_hashes_match'] = all(
        sha256(OUTPUT / name) == sha256(DELIVERY / name) for name in [DOCX_NAME, PDF_NAME]
    )
    checks['formal_validation_pass'] = json.loads(
        (VALIDATION / 'validation_summary_v1.json').read_text(encoding='utf-8')
    )['status'] == 'PASS'
    result = {
        'status': 'PASS' if all(checks.values()) else 'FAIL',
        'generated_at': NOW,
        'checks': checks,
        'phase113_manifest_rows': len(phase_rows),
        'delivery_directory': str(DELIVERY),
        'active_jobs': 'none',
    }
    out = VALIDATION / 'PHASE113_STATE_UPDATE_VALIDATION_v1.json'
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')

    manifest = PROJECT / 'RESULTS_MANIFEST.tsv'
    with manifest.open('r', encoding='utf-8-sig', newline='') as handle:
        existing_ids = {row['result_id'] for row in csv.DictReader(handle, delimiter='\t')}
    if 'phase113_state_update_validation_v1' not in existing_ids:
        with manifest.open('a', encoding='utf-8', newline='') as handle:
            csv.writer(handle, delimiter='\t', lineterminator='\n').writerow([
                'phase113_state_update_validation_v1', 'phase113_project_state_update',
                'state_update_validation_json', str(out), NOW, result['status'],
                'Phase113 mandatory four-file completion gate; all recorded checks PASS',
                'State validation is written before its manifest registration to avoid self-hash recursion',
            ])
    if result['status'] != 'PASS':
        raise RuntimeError(json.dumps(result, ensure_ascii=False))
    return out, result


def main():
    delivery_manifest = make_delivery_manifest()
    file_manifest = make_file_manifest(delivery_manifest)
    backup = backup_state()
    update_state(file_manifest, delivery_manifest, backup)
    state_validation, result = validate_state()
    print(json.dumps({'state_validation': str(state_validation), **result}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
