import csv
import hashlib
import json
import shutil
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo


PROJECT = Path(r'D:\project26210_external_joint_v1')
ROOT = PROJECT / 'phase115_supervisor_evidence_review_local_v1'
OUTPUT = ROOT / 'output_v1'
VALIDATION = ROOT / 'validation_v1'
DELIVERY = Path(r'C:\Users\31937\Downloads\项目26210_主管会议证据简报_20260901_v1')
DOCX_NAME = '项目26210_现有证据总览与下一步决策简报_v1.docx'
PDF_NAME = '项目26210_现有证据总览与下一步决策简报_v1.pdf'
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
        source_hash = sha256(source)
        delivered_hash = sha256(delivered)
        if source_hash != delivered_hash:
            raise RuntimeError(f'delivery mismatch: {name}')
        rows.append([name, str(source), str(delivered), source_hash, source.stat().st_size, 'PASS_EXACT_COPY'])
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['file_name', 'authoritative_path', 'delivery_path', 'sha256', 'size_bytes', 'status'])
        writer.writerows(rows)
    return path


def formal_paths(delivery_manifest):
    paths = [
        ROOT / 'artifact.md',
        ROOT / 'REVISION_AUDIT_v1.md',
        ROOT / 'assets_v1' / 'current_evidence_chain_status_v1.png',
        ROOT / 'scripts_v1' / '01_build_supervisor_evidence_review_v1.py',
        ROOT / 'scripts_v1' / '02_validate_supervisor_evidence_review_v1.py',
        ROOT / 'scripts_v1' / '03_complete_phase115_state_v1.py',
        OUTPUT / DOCX_NAME,
        OUTPUT / PDF_NAME,
        ROOT / 'render_final_v3' / 'render_status.json',
        VALIDATION / 'DOCUMENT_VALIDATION_v1.json',
        VALIDATION / 'a11y_audit_final_v1.json',
        delivery_manifest,
        DELIVERY / DOCX_NAME,
        DELIVERY / PDF_NAME,
    ]
    return paths


def make_file_manifest(paths):
    path = ROOT / 'FILE_MANIFEST_SHA256_v1.tsv'
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['absolute_path', 'sha256', 'size_bytes', 'status'])
        for item in paths:
            if not item.exists():
                raise FileNotFoundError(str(item))
            writer.writerow([str(item), sha256(item), item.stat().st_size, 'PASS'])
        for png in sorted((ROOT / 'render_final_v3').glob('page-*.png')):
            writer.writerow([str(png), sha256(png), png.stat().st_size, 'PASS_RENDER_PAGE'])
    return path


def backup_state():
    backup = ROOT / 'state_backup_before_phase115_v1'
    index = 2
    while backup.exists():
        backup = ROOT / f'state_backup_before_phase115_v{index}'
        index += 1
    backup.mkdir()
    for name in STATE_FILES:
        shutil.copy2(PROJECT / name, backup / name)
    return backup


def update_state(backup):
    status_block = f'''## Phase 115 supervisor evidence review and decision brief (2026-09-01)

- Status: **PASS_EVIDENCE_SYNTHESIS_COMPLETE; WORD RENDER 11/11; DOCUMENT VALIDATION PASS; ACCESSIBILITY HIGH=0**. No new biological analysis was executed.
- Purpose: create an internal supervisor-meeting brief that fully consolidates supported, partial, negative, contradictory and non-estimable evidence, assesses the latest client scope change and presents explicit next-step options.
- Inherited independent biological units: original D01 single-cell patients up to 40; Phase109 bridge patients n=20 abundance/polarization and n=14 communication; GSE175540 12 FFPE plus 12 frozen sections; UnitedMet tumour cohorts including Phase114 complete S7 cases n=26/13; NCI60 n=58 and renal n=8; small immune-treatment single-cell response/paired subsets; external bulk ICI patient cohorts. This phase creates no new statistical sample.
- Principal synthesis: F1-TYMP RNA association is supported with study heterogeneity. Complete TYMP-MTA, MTA-SPP1, general TYMP/F1-SPP1, immune-exclusion and clinical-translation edges are not closed or not supported. The reproducible positive downstream anchor is the SPP1-state macrophage-endothelial/vascular interface.
- Client-scope assessment: the request to demote metabolism and organize results by responder/non-responder immune ecology is a legitimate new direction, but `TYMP through SPP1 macrophages` remains a hypothesis because prior patient-level, communication and spatial tests were negative.
- Recommended meeting option: use the supported SPP1-endothelial/vascular ecology as the report anchor and, only after supervisor approval, run a minimal response-stratified patient-level module in D02/D03. TYMP-SPP1 must be tested as a bridge rather than assumed. Metabolic UMAP panels may show transcriptional pathway proxies only.
- Authoritative deliverables: `{OUTPUT / DOCX_NAME}` and `{OUTPUT / PDF_NAME}`. Downloads: `{DELIVERY}`. Validation: `{VALIDATION / 'DOCUMENT_VALIDATION_v1.json'}` and `{VALIDATION / 'a11y_audit_final_v1.json'}`. Pre-update backup: `{backup}`.
- Current breakpoint: hold further biological analysis until the supervisor chooses the narrative target and authorizes a route. Active jobs: none.'''
    append_once(PROJECT / 'PROJECT_STATUS.md', '## Phase 115 supervisor evidence review and decision brief', status_block)

    decisions_block = '''## 2026-09-01 - Phase 115 supervisor-meeting evidence decisions

- Preserve every material negative and contradiction in the internal meeting brief; do not optimize the brief around the desired TYMP-SPP1 conclusion.
- Treat the latest client feedback as a proposed scope change from metabolism-gated transfer to response-stratified immune ecology, not as evidence that the new mechanism is already supported.
- Retain the Phase114 seven-feature result as a supplemental negative module. Do not continue generic metabolomic sensitivity expansion unless a new scope is explicitly approved.
- Use SPP1-state macrophage-endothelial/vascular evidence as the strongest current downstream anchor.
- Do not state that TYMP acts through SPP1 macrophages unless a new response-stratified module supports TYMP/F1-response, TYMP/F1-SPP1 and SPP1-T-cell/response edges under patient-level inference.
- If the bridge remains negative, present TYMP/F1 tumour state and SPP1-vascular ecology as parallel modules rather than rescuing the chain.
- Single-cell UMAP can display transcriptional pathway scores, not measured 5-MTA, thymidine or metabolic flux.
- No new route is authorized by the brief itself. Supervisor approval is required before response-stratified reanalysis or additional spatial work.'''
    append_once(PROJECT / 'DECISIONS.md', '## 2026-09-01 - Phase 115 supervisor-meeting evidence decisions', decisions_block)

    handoff_block = f'''## Phase 115 supervisor evidence review current breakpoint (2026-09-01)

- Completed: 11-page internal DOCX/PDF evidence and decision brief with dataset-universe separation, complete chain matrix, exact statistics, client-feedback assessment, allowed/prohibited claims, five route options, a recommended hybrid and a minimum go/no-go analysis block.
- Validation: Word render 11/11 with no blank page; fixed A4 geometry; 19 tables and 12 multirow evidence/decision tables; all multirow tables have repeated header rows; one figure has meaningful alt text; document validation PASS; accessibility high=0, medium=8 layout-table notices, low=0; downloaded DOCX/PDF hashes match authoritative files.
- Exact resume point: meet with the supervisor. Choose whether to approve the recommended `SPP1-vascular anchor + minimal response-stratified immune analysis`, another route, or no further computation.
- Active jobs: none. No new biological result was generated in Phase115.
- Key paths: `{OUTPUT / DOCX_NAME}`, `{OUTPUT / PDF_NAME}`, `{VALIDATION / 'DOCUMENT_VALIDATION_v1.json'}`, `{DELIVERY}`.
- Claims that must not be made: complete TYMP-MTA-SPP1 chain; general TYMP-driven SPP1 polarization; established T-cell killing impairment or physical immune barrier; current ICB/PFS/OS predictor; direct metabolite measurement from transcriptional UMAP scores.'''
    append_once(PROJECT / 'CODEX_HANDOFF.md', '## Phase 115 supervisor evidence review current breakpoint', handoff_block)


def append_manifest(paths, file_manifest):
    manifest = PROJECT / 'RESULTS_MANIFEST.tsv'
    with manifest.open('r', encoding='utf-8-sig', newline='') as handle:
        existing_ids = {row['result_id'] for row in csv.DictReader(handle, delimiter='\t')}
    artifacts = [
        ('phase115_design_contract_v1', 'design_contract_md', ROOT / 'artifact.md'),
        ('phase115_revision_audit_v1', 'revision_audit_md', ROOT / 'REVISION_AUDIT_v1.md'),
        ('phase115_evidence_chain_figure_v1', 'figure_png', ROOT / 'assets_v1' / 'current_evidence_chain_status_v1.png'),
        ('phase115_report_builder_v1', 'formal_script', ROOT / 'scripts_v1' / '01_build_supervisor_evidence_review_v1.py'),
        ('phase115_report_validator_v1', 'formal_validation_script', ROOT / 'scripts_v1' / '02_validate_supervisor_evidence_review_v1.py'),
        ('phase115_state_completion_script_v1', 'formal_script', ROOT / 'scripts_v1' / '03_complete_phase115_state_v1.py'),
        ('phase115_supervisor_brief_docx_v1', 'internal_meeting_brief_docx', OUTPUT / DOCX_NAME),
        ('phase115_supervisor_brief_pdf_v1', 'rendered_meeting_brief_pdf', OUTPUT / PDF_NAME),
        ('phase115_render_status_v1', 'render_status_json', ROOT / 'render_final_v3' / 'render_status.json'),
        ('phase115_document_validation_v1', 'validation_json', VALIDATION / 'DOCUMENT_VALIDATION_v1.json'),
        ('phase115_a11y_audit_v1', 'accessibility_json', VALIDATION / 'a11y_audit_final_v1.json'),
        ('phase115_delivery_manifest_v1', 'delivery_manifest_tsv', ROOT / 'DELIVERY_MANIFEST_v1.tsv'),
        ('phase115_file_manifest_v1', 'file_manifest_tsv', file_manifest),
        ('phase115_download_docx_v1', 'delivered_docx', DELIVERY / DOCX_NAME),
        ('phase115_download_pdf_v1', 'delivered_pdf', DELIVERY / PDF_NAME),
    ]
    new_rows = []
    for result_id, artifact_type, path in artifacts:
        if result_id in existing_ids:
            continue
        evidence = f'SHA256 {sha256(path)}; Word render 11/11; document validation PASS; accessibility high=0; no new biological analysis'
        new_rows.append([result_id, 'phase115_supervisor_evidence_review', artifact_type, str(path), NOW, 'PASS', evidence, 'Internal supervisor evidence synthesis and decision brief'])
    with manifest.open('a', encoding='utf-8', newline='') as handle:
        csv.writer(handle, delimiter='\t', lineterminator='\n').writerows(new_rows)
    return len(new_rows)


def validate_state(delivery_manifest, expected_rows):
    checks = {f'exists_{name}': (PROJECT / name).exists() for name in STATE_FILES}
    with (PROJECT / 'RESULTS_MANIFEST.tsv').open('r', encoding='utf-8-sig', newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    ids = [row['result_id'] for row in rows]
    phase_rows = [row for row in rows if row['stage'] == 'phase115_supervisor_evidence_review']
    checks['manifest_schema'] = list(rows[0]) == ['result_id', 'stage', 'artifact_type', 'absolute_path', 'generated_at', 'validation_status', 'validation_evidence', 'notes']
    checks['manifest_unique_ids'] = len(ids) == len(set(ids))
    checks['phase115_rows'] = len(phase_rows) == expected_rows
    checks['phase115_paths_exist'] = all(Path(row['absolute_path']).exists() for row in phase_rows)
    checks['completion_markers'] = all(
        marker in (PROJECT / name).read_text(encoding='utf-8-sig')
        for name, marker in [
            ('PROJECT_STATUS.md', 'Phase 115 supervisor evidence review and decision brief'),
            ('DECISIONS.md', 'Phase 115 supervisor-meeting evidence decisions'),
            ('CODEX_HANDOFF.md', 'Phase 115 supervisor evidence review current breakpoint'),
        ]
    )
    checks['document_validation_pass'] = json.loads((VALIDATION / 'DOCUMENT_VALIDATION_v1.json').read_text(encoding='utf-8'))['status'] == 'PASS'
    checks['delivery_exact'] = all(row['status'] == 'PASS_EXACT_COPY' for row in csv.DictReader(delivery_manifest.open(encoding='utf-8'), delimiter='\t'))
    result = {
        'status': 'PASS' if all(checks.values()) else 'FAIL',
        'generated_at': NOW, 'checks': checks,
        'phase115_manifest_rows': len(phase_rows),
        'delivery_directory': str(DELIVERY), 'active_jobs': 'none',
    }
    out = VALIDATION / 'PHASE115_STATE_UPDATE_VALIDATION_v1.json'
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    if result['status'] != 'PASS':
        raise RuntimeError(json.dumps(result, ensure_ascii=False))
    manifest = PROJECT / 'RESULTS_MANIFEST.tsv'
    with manifest.open('r', encoding='utf-8-sig', newline='') as handle:
        ids = {row['result_id'] for row in csv.DictReader(handle, delimiter='\t')}
    if 'phase115_state_update_validation_v1' not in ids:
        with manifest.open('a', encoding='utf-8', newline='') as handle:
            csv.writer(handle, delimiter='\t', lineterminator='\n').writerow([
                'phase115_state_update_validation_v1', 'phase115_project_state_update', 'state_update_validation_json', str(out), NOW, 'PASS',
                'Phase115 mandatory four-file completion gate PASS; Word render 11/11; delivery exact',
                'State validation generated before self-registration to avoid self-hash recursion',
            ])
    return out, result


def main():
    delivery_manifest = make_delivery_manifest()
    paths = formal_paths(delivery_manifest)
    file_manifest = make_file_manifest(paths)
    backup = backup_state()
    update_state(backup)
    new_rows = append_manifest(paths, file_manifest)
    state_validation, result = validate_state(delivery_manifest, new_rows)
    print(json.dumps({'state_validation': str(state_validation), **result}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
