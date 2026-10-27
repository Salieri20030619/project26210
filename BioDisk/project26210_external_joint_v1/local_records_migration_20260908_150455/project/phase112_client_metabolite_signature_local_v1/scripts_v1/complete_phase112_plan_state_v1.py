from pathlib import Path
from datetime import datetime
import csv, hashlib, json, shutil

ROOT=Path(r'D:\project26210_external_joint_v1')
PKG=ROOT/'phase112_client_metabolite_signature_local_v1'
STATES=[ROOT/'PROJECT_STATUS.md',ROOT/'DECISIONS.md',ROOT/'RESULTS_MANIFEST.tsv',ROOT/'CODEX_HANDOFF.md']
NOW=datetime.now().astimezone().isoformat(timespec='seconds')

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as fh:
        for b in iter(lambda:fh.read(1048576),b''): h.update(b)
    return h.hexdigest().upper()

def append_once(path,marker,text):
    old=path.read_text(encoding='utf-8')
    if marker in old: raise RuntimeError('marker exists: '+marker)
    path.write_text(old.rstrip()+'\n\n'+text.strip()+'\n',encoding='utf-8')

def main():
    if not all(p.exists() for p in STATES): raise FileNotFoundError('state file missing')
    backup=PKG/('preupdate_backup_'+datetime.now().strftime('%Y%m%d_%H%M%S')); backup.mkdir()
    for p in STATES: shutil.copy2(p,backup/p.name)

    append_once(ROOT/'PROJECT_STATUS.md','## Phase 112 client-review seven-feature analysis plan',f'''## Phase 112 client-review seven-feature analysis plan (2026-09-01)

- Status: **PASS_CLIENT_ANALYSIS_PLAN_V1; WORD RENDER 16/16 PAGES; FORMAL VALIDATION 15/15 PASS; ACCESSIBILITY HIGH=0, MEDIUM=0, LOW=0; NO SEVEN-FEATURE RESULT COMPUTED**.
- Scope: planning and feature registration only. The client figure is used to freeze an unweighted shTYMP-like seven-metabolite direction score: 5-MTA down, Thymidine up, 2'-Deoxyuridine up, Methionine down, SAM down, Thymine down and Spermidine down. No value, effect weight or replicate count was estimated from the image.
- Planned independent biological units: UnitedMet RC18 tumour n=97 and RC20 tumour n=52 are the complete 7/7 primary cohorts; CPTAC tumour n=50 and CPTAC_val tumour n=56 are 3/7 sensitivity cohorts; NCI60 uses independent cell lines with all-lineage n up to 58 and renal n=8 as a 5/7 sensitivity.
- Planned methods: within-dataset rank z-scaling; complete-case equal-weight score; separate RC18/RC20 Spearman tests and two-test BH family; seven-feature component BH families; MTAP-adjusted rank sensitivity; Fisher-z meta-analysis; leave-one-feature/sample influence; 100,000 permutations; reduced datasets cannot upgrade the full-score decision.
- Conditional modules are explicit: public perturbation search; RNA surrogate only after complete-score transfer; single-cell/spatial projection only after bidirectional RNA transfer; clinical projection only after downstream support. No new wet-lab work is planned.
- Evidence boundary: client experiment, public observational transfer and downstream projection remain separate evidence layers. Public failure cannot be replaced by the client figure, and the client figure cannot be treated as a public replication.
- Authoritative package: `{PKG}`. Client DOCX: `{PKG/'output_v1'/'项目26210_无权重七特征公共数据迁移分析方案_v1.docx'}`. PDF: `{PKG/'output_v1'/'项目26210_无权重七特征公共数据迁移分析方案_v1.pdf'}`.
- Validation: 16 Word-rendered pages, no blank pages, safe page margins, eight header-marked tables, one alt-text figure, required steps present, no internal Phase labels, 15/15 checks PASS; accessibility 0/0/0. Interactive image inspection remained unavailable because the Windows image helper failed with token error 1344.
- Current breakpoint: send the plan for client review. Do not execute steps 0-7 until approval; steps 8-11 remain conditional and are not automatically authorized.''')

    append_once(ROOT/'DECISIONS.md','## 2026-09-01 - Phase 112 client-feature transfer planning decision','''## 2026-09-01 - Phase 112 client-feature transfer planning decision

- Supersede the Phase111 assumption that additional wet-lab work or the spatial route is the immediate upstream solution. The client is not expected to add wet-lab experiments.
- Treat the provided client figure as a client-owned feature-registration source, not as a replacement for failed public validation and not as a public-data sample.
- Freeze the unweighted seven-feature shTYMP-like score before analysis. Do not estimate weights from the image and do not change signs or membership after viewing public outcomes.
- Use RC18 and RC20 as the only complete 7/7 primary cohorts. CPTAC/CPTAC_val 3/7 and NCI60 5/7 analyses are explicitly reduced sensitivities and cannot validate the full score.
- Preserve the prior TYMP-MTA null result as a separate endpoint. A composite score cannot overwrite it.
- Execute only steps 0-7 after client approval. RNA, single-cell, spatial and clinical projection require their stated upstream gates. Failure at any gate is a completed result, not a reason to add models or datasets.
- Do not use the two unprovided client figures, infer their content or include them in the feature registry.''')

    append_once(ROOT/'CODEX_HANDOFF.md','## Phase 112 seven-feature plan current breakpoint',f'''## Phase 112 seven-feature plan current breakpoint (2026-09-01)

- Completed: client-facing 16-page DOCX/PDF analysis plan, seven-feature registry and source-image identity record. No seven-feature statistic or biological result has been computed.
- Exact resume point: obtain client approval. If approved, execute plan steps 0-7 only: input/provenance audit, exact metabolite mapping, pairing, unweighted scoring, RC18/RC20 primary tests, cross-cohort synthesis, robustness and reduced sensitivities.
- Active jobs: none. Analysis directories remain empty. Conditional steps 8-11 are not authorized by plan approval alone unless their upstream gates pass.
- Key paths: `{PKG/'output_v1'}`, `{PKG/'inputs_v1'/'CLIENT_FIGURE_FEATURE_REGISTRY_v1.json'}`, `{PKG/'qc_v1'/'CLIENT_ANALYSIS_PLAN_VALIDATION_v1.json'}`.
- Claims that must not be made: public replication of the client perturbation; a validated TYMP-MTA-SPP1 chain; full seven-feature validation from 3/7 or 5/7 data; any result from the two unprovided client figures.''')

    manifest=ROOT/'RESULTS_MANIFEST.tsv'
    with manifest.open('r',encoding='utf-8',newline='') as fh: existing=list(csv.DictReader(fh,delimiter='\t'))
    ids={r['result_id'] for r in existing}
    artifacts=[
        ('phase112_client_analysis_plan_docx_v1','client_analysis_plan_docx',PKG/'output_v1'/'项目26210_无权重七特征公共数据迁移分析方案_v1.docx','PASS'),
        ('phase112_client_analysis_plan_pdf_v1','rendered_plan_pdf',PKG/'output_v1'/'项目26210_无权重七特征公共数据迁移分析方案_v1.pdf','PASS'),
        ('phase112_seven_feature_workflow_png_v1','figure_png',PKG/'assets_v1'/'seven_feature_analysis_workflow_v1.png','PASS'),
        ('phase112_client_fig4_source_v1','client_owned_source_image',PKG/'inputs_v1'/'client_fig4_metabolomics_source_20260830.png','PASS_IDENTITY_ONLY'),
        ('phase112_client_feature_registry_v1','feature_registry_json',PKG/'inputs_v1'/'CLIENT_FIGURE_FEATURE_REGISTRY_v1.json','PASS'),
        ('phase112_plan_builder_v1','formal_script',PKG/'scripts_v1'/'build_client_analysis_plan_v1.py','PASS'),
        ('phase112_plan_validator_v1','formal_validation_script',PKG/'scripts_v1'/'validate_client_analysis_plan_v1.py','PASS'),
        ('phase112_plan_validation_v1','validation_json',PKG/'qc_v1'/'CLIENT_ANALYSIS_PLAN_VALIDATION_v1.json','PASS'),
        ('phase112_plan_a11y_v1','accessibility_json',PKG/'qc_v1'/'a11y_audit_v1.json','PASS'),
        ('phase112_plan_render_status_v1','render_status_json',PKG/'qc_v1'/'render_pages'/'render_status.json','PASS'),
        ('phase112_file_manifest_v1','file_manifest_tsv',PKG/'FILE_MANIFEST_SHA256_v1.tsv','PASS'),
        ('phase112_plan_state_script_v1','formal_script',PKG/'scripts_v1'/'complete_phase112_plan_state_v1.py','PASS'),
    ]
    for rid,_,path,_ in artifacts:
        if rid in ids or not path.exists(): raise RuntimeError(f'bad artifact {rid}: {path}')
    rows=[]
    for rid,kind,path,status in artifacts:
        note='Client-owned input identity only; no reanalysis' if status=='PASS_IDENTITY_ONLY' else 'Client-review unweighted seven-feature transfer plan; no analysis executed'
        rows.append([rid,'phase112_client_metabolite_signature_plan',kind,str(path),NOW,status,f'SHA256 {sha(path)}; Word render 16/16; validation 15/15; a11y 0/0/0; no seven-feature result computed',note])
    with manifest.open('a',encoding='utf-8',newline='') as fh: csv.writer(fh,delimiter='\t',lineterminator='\n').writerows(rows)

    with manifest.open('r',encoding='utf-8',newline='') as fh: final=list(csv.DictReader(fh,delimiter='\t'))
    final_ids=[r['result_id'] for r in final]
    checks={
        'four_state_files_exist':all(p.exists() for p in STATES),
        'manifest_schema':list(final[0].keys())==['result_id','stage','artifact_type','absolute_path','generated_at','validation_status','validation_evidence','notes'],
        'manifest_unique_ids':len(final_ids)==len(set(final_ids)),
        'phase112_rows_present':all(rid in final_ids for rid,_,_,_ in artifacts),
        'phase112_paths_exist':all(path.exists() for _,_,path,_ in artifacts),
        'phase112_hashes_match':all(f'SHA256 {sha(path)}' in next(r['validation_evidence'] for r in final if r['result_id']==rid) for rid,_,path,_ in artifacts),
        'no_analysis_outputs':not any((PKG/'analysis_v1').iterdir()),
    }
    report={'status':'PASS' if all(checks.values()) else 'FAIL','updated_at':NOW,'backup':str(backup),'rows_added':len(rows),'checks':checks}
    out=PKG/'qc_v1'/'PHASE112_PLAN_STATE_VALIDATION_v1.json'; out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))
    if report['status']!='PASS': raise SystemExit(1)

if __name__=='__main__': main()
