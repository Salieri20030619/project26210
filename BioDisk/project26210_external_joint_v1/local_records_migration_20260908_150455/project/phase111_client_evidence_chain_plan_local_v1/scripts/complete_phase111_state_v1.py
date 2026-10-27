from pathlib import Path
from datetime import datetime
import csv, hashlib, json, shutil

ROOT = Path(r'D:\project26210_external_joint_v1')
PKG = ROOT / 'phase111_client_evidence_chain_plan_local_v1'
STATES = [ROOT/'PROJECT_STATUS.md', ROOT/'DECISIONS.md', ROOT/'RESULTS_MANIFEST.tsv', ROOT/'CODEX_HANDOFF.md']
NOW = datetime.now().astimezone().isoformat(timespec='seconds')

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for b in iter(lambda: fh.read(1048576), b''):
            h.update(b)
    return h.hexdigest().upper()

def append_once(path, marker, text):
    old = path.read_text(encoding='utf-8')
    if marker in old:
        raise RuntimeError(f'marker already present: {marker}')
    path.write_text(old.rstrip() + '\n\n' + text.strip() + '\n', encoding='utf-8')

def main():
    for p in STATES:
        if not p.exists():
            raise FileNotFoundError(p)
    backup = PKG / ('preupdate_backup_' + datetime.now().strftime('%Y%m%d_%H%M%S'))
    backup.mkdir()
    for p in STATES:
        shutil.copy2(p, backup / p.name)

    append_once(ROOT/'PROJECT_STATUS.md', '## Phase 111 evidence-first dual-route client plan', f'''## Phase 111 evidence-first dual-route client plan (2026-09-01)

- Status: **PASS_CLIENT_EVIDENCE_CHAIN_PLAN_V2; WORD RENDER 10/10 PAGES; FORMAL VALIDATION 15/15 PASS; ACCESSIBILITY HIGH=0, MEDIUM=0, LOW=0**.
- Scope: client-facing planning/reporting revision only; no new biological analysis was executed. The report is rebuilt around two explicit routes: F1/TYMP -> MTA -> SPP1 -> spatial structure, and F1/TYMP -> MTA -> SPP1 -> prognosis/treatment.
- Independent biological units are inherited and named by evidence source: original single-cell patients (up to 40), Phase109 bridge patients (20 abundance;14 communication), GSE175540 sections (12 FFPE plus 12 frozen), UnitedMet tumour samples (n=255), NCI60 cell lines (n=58;renal n=8), corrected CheckMate RNA patients (n=311) and JAVELIN RNA patients (n=726). The report creates no new statistical sample.
- Principal positive evidence retained: F1-TYMP association with heterogeneity; operational SPP1 macrophage state; bidirectional SPP1-endothelial communication compatibility (14/14 and 13/14 positive; both q=0.00305); GSE175540 SPP1-vascular support in FFPE and frozen.
- Principal negative/non-estimable evidence retained: public data do not support TYMP-driven measured MTA; direct MTA-induced SPP1 is unproven; original Part1 does not support a general TYMP/F1-to-SPP1 edge; existing response/PFS/OS/treatment-interaction analyses do not support clinical translation.
- Recommended breakpoint: execute the spatial route first (immutable input audit, outcome-blind niche transfer, capped LR spatial validation, official-TLS immune gradients, then section-level niche pseudobulk). Clinical projection is conditional on a frozen reproducible spatial phenotype. Direct TYMP-MTA-SPP1 closure remains perturbational/client-owned.
- Authoritative package: `{PKG}`. Final editable report: `{PKG/'output_v2'/'项目26210_双路线证据与补充分析计划_v2.docx'}`. Final PDF: `{PKG/'output_v2'/'项目26210_双路线证据与补充分析计划_v2.pdf'}`.
- Validation: Microsoft Word fixed-format render produced 10 pages/10 PNGs; page-pixel audit found no blank page and safe page margins; required evidence language, heading hierarchy, table headers, alt text and absence of internal Phase labels/placeholders passed 15/15. Interactive image inspection remained unavailable because the Windows sandbox image helper repeatedly failed with token error 1344.''')

    append_once(ROOT/'DECISIONS.md', '## 2026-09-01 - Phase 111 evidence-first dual-route plan decision', '''## 2026-09-01 - Phase 111 evidence-first dual-route plan decision

- Supersede the V1 activity-led plan with the V2 evidence-first client plan.
- Organize the next-stage report as two explicit routes: spatial structure and prognosis/treatment; each route must state existing evidence, current decision, missing edge, possible analysis, success criterion and stop rule.
- Preserve unsupported TYMP-to-MTA, unclosed MTA-to-SPP1 and negative clinical translation results in their biological positions. Downstream SPP1-endothelial/vascular support cannot rescue upstream gaps.
- Prioritize the 24-section spatial route. Do not start a new generic clinical model until a reproducible outcome-blind spatial phenotype is frozen.
- Do not expose internal Phase identifiers in the client document. Keep datasets and statistical units separate and retain patient/section/cell-line units and FDR families.
- Direct causal closure requires TYMP perturbation, direct MTA measurement, SPP1 RNA/protein and A2B/STAT3 blockade or rescue; computation alone cannot close these edges.''')

    append_once(ROOT/'CODEX_HANDOFF.md', '## Phase 111 evidence-first dual-route plan current breakpoint', f'''## Phase 111 evidence-first dual-route plan current breakpoint (2026-09-01)

- Completed: non-overwriting V2 DOCX/PDF rebuilt around the spatial and prognosis/treatment routes; no biological analysis rerun.
- Validation: 10/10 Word-rendered pages; 15/15 formal structural/page checks PASS; accessibility 0/0/0; source and final hashes recorded in the package manifest.
- Exact resume point: obtain user review of V2. If analysis is authorized, start with the GSE175540 24-section immutable input/metadata audit and preregistered outcome-blind spatial niche block. Do not begin clinical projection before a spatial phenotype is frozen.
- Active jobs: none. Blockers: direct TYMP-MTA and MTA-SPP1 causal closure requires eligible perturbation data/client experiments; interactive page-image inspection is unavailable in the current Windows sandbox because of token error 1344.
- Key package: `{PKG}`. Claims that must not be made: validated complete chain; direct MTA-induced SPP1; general TYMP-driven SPP1 polarization; current ICI response/PFS/OS predictor.''')

    manifest = ROOT/'RESULTS_MANIFEST.tsv'
    with manifest.open('r', encoding='utf-8', newline='') as fh:
        existing = list(csv.DictReader(fh, delimiter='\t'))
    ids = {r['result_id'] for r in existing}
    artifacts = [
        ('phase111_client_evidence_plan_docx_v2','client_report_docx',PKG/'output_v2'/'项目26210_双路线证据与补充分析计划_v2.docx'),
        ('phase111_client_evidence_plan_pdf_v2','rendered_report_pdf',PKG/'output_v2'/'项目26210_双路线证据与补充分析计划_v2.pdf'),
        ('phase111_two_route_evidence_map_v2','figure_png',PKG/'figures_v2'/'two_route_evidence_map_v2.png'),
        ('phase111_report_builder_v2','formal_script',PKG/'scripts'/'build_client_evidence_chain_plan_v2.py'),
        ('phase111_report_validator_v2','formal_validation_script',PKG/'scripts'/'validate_evidence_chain_report_v2.py'),
        ('phase111_report_validation_v2','validation_json',PKG/'qc_v2'/'REPORT_VALIDATION_v2.json'),
        ('phase111_a11y_validation_v2','accessibility_json',PKG/'qc_v2'/'a11y_audit_final_v2.json'),
        ('phase111_word_render_status_v2','render_status_json',PKG/'qc_v2'/'render_status.json'),
        ('phase111_file_manifest_v1','file_manifest_tsv',PKG/'FILE_MANIFEST_SHA256_v1.tsv'),
        ('phase111_state_completion_script_v1','formal_script',PKG/'scripts'/'complete_phase111_state_v1.py'),
    ]
    for rid, _, path in artifacts:
        if rid in ids or not path.exists():
            raise RuntimeError(f'bad manifest target: {rid} {path}')
    rows = []
    for rid, kind, path in artifacts:
        rows.append([rid,'phase111_client_evidence_chain_plan',kind,str(path),NOW,'PASS',f'SHA256 {sha(path)}; Word render 10/10; formal validation 15/15; a11y 0/0/0','Evidence-first spatial and prognosis/treatment route plan'])
    with manifest.open('a', encoding='utf-8', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t', lineterminator='\n')
        writer.writerows(rows)

    report = {'status':'PASS','updated_at':NOW,'backup':str(backup),'rows_added':len(rows),'checks':[]}
    with manifest.open('r', encoding='utf-8', newline='') as fh:
        final = list(csv.DictReader(fh, delimiter='\t'))
    final_ids = [r['result_id'] for r in final]
    checks = {
        'four_state_files_exist': all(p.exists() for p in STATES),
        'manifest_schema': list(final[0].keys()) == ['result_id','stage','artifact_type','absolute_path','generated_at','validation_status','validation_evidence','notes'],
        'manifest_unique_ids': len(final_ids) == len(set(final_ids)),
        'phase111_rows_present': all(rid in final_ids for rid,_,_ in artifacts),
        'phase111_paths_exist': all(path.exists() for _,_,path in artifacts),
        'phase111_hashes_match': all(f'SHA256 {sha(path)}' in next(r['validation_evidence'] for r in final if r['result_id']==rid) for rid,_,path in artifacts),
    }
    report['checks'] = checks
    if not all(checks.values()):
        report['status'] = 'FAIL'
    out = PKG/'qc_v2'/'PHASE111_STATE_UPDATE_VALIDATION_v1.json'
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    if report['status'] != 'PASS':
        raise SystemExit(1)
    print(json.dumps(report, ensure_ascii=False, indent=2))

if __name__ == '__main__':
    main()
