from __future__ import annotations
import csv,hashlib,re
from datetime import datetime
from pathlib import Path
STAGE=Path(r"C:\Users\31937\project26210_phase108_state_update_v1");ROOT=Path(r"D:\project26210_external_joint_v1")
STATUS=r"""

## Phase 107-108 evidence escalation and unified client architecture (2026-08-31)

- Phase107 status: **PASS_EVIDENCE_ESCALATION_VALIDATION 26/26**. The rebuilt original Part1 supported component associations but not the full chain; later immune single-cell, spatial, measured-metabolite and cell-line datasets were targeted tests of specific unresolved edges, not additional samples in one pooled cohort.
- Phase108 status: **PASS_UNIFIED_CHAIN_VALIDATION 22/22**. The final client-facing structure retires the original Part1-4 divisions and uses one continuous biological chain: TYMP/F1 tumour state -> metabolic bridge -> macrophage subdivision/SPP1 state -> spatial vascular/immune ecology -> clinical translation -> client perturbational validation.
- Phase106 dataset separation remains the internal control layer. Dataset names, biological units, endpoints and multiplicity families remain explicit in Methods, captions and evidence tables; effects are not pooled across incompatible datasets.
- Nine chain edges were audited in Phase107. F1-TYMP RNA association is supported with heterogeneity; strict TYMP-specific function is partial/unsupported; TYMP-measured MTA is unsupported; MTA-SPP1 is not closed; general tumour-TYMP-to-SPP1 polarization is unsupported; GSE175540 provides dataset-specific SPP1-vascular support; vascular-immune exclusion and ICB/PFS/OS translation are unsupported.
- Later datasets remain narratively connected because each addresses a Part1 gap. Positive downstream edges do not rescue negative upstream edges, and the complete chain is not claimed as validated.
- The report now has seven result-first sections and six planned main figures. It does not present datasets as separate client-facing chapters and does not use original Part labels.
- Client validation remains a labelled perturbational interface. No client data or predicted outcome is included; client owns final experimental figures.
- `CLIENT_REPORT_ARCHITECTURE_v4.md` supersedes Phase105 v1, Phase106 v2 and the two-option organisation in Phase107 v3. Phase106/107 audit tables remain authoritative provenance.
- Current breakpoint: use the unified v4 architecture for content reconstruction. Each chain edge must display one of supported,context-dependent/partial,current-data-do-not-support,not-evaluable,or pending-client-validation.
"""
DECISIONS=r"""

## 2026-08-31 - Phase 107-108 unified-chain decision

- Retire Part1-4 as client-facing report divisions. Organize all final content by biological chain nodes and edges.
- Keep Phase106 dataset separation for statistical governance while using Phase107 gap-driven follow-up logic to connect datasets narratively.
- Do not pool incompatible datasets, but do not present them as disconnected project chapters. Introduce each dataset where it tests a specific unresolved edge.
- Retain negative, contradictory and non-estimable edges in their biological position. Do not search for or substitute another dataset solely to obtain a preferred direction.
- Do not claim the complete chain is validated. Report component associations and exact edge states.
- Use the seven-section structure and six-figure plan in Phase108 as the current client-report architecture.
- Reserve client-owned validation slots for direct TYMP-to-MTA and MTA-to-SPP1 perturbational evidence without importing,predicting or redrawing client data.
"""
HANDOFF=r"""

## Phase 107-108 current unified-chain state (2026-08-31)

- Phase107 validation 26/26 PASS; Phase108 validation 22/22 PASS.
- Authoritative directories: `phase107_evidence_escalation_audit_local_v1` and `phase108_unified_chain_architecture_local_v1`.
- Final client report has no Part1-4 division. Use `CLIENT_REPORT_ARCHITECTURE_v4.md` and `UNIFIED_CHAIN_NODE_MATRIX_v1.tsv`.
- Narrative chain: TYMP/F1 tumour state -> metabolic bridge -> macrophage subdivision/SPP1 state -> spatial vascular/immune ecology -> clinical translation -> client validation.
- Evidence control: Phase106 datasets remain statistically separate; Phase107 links each follow-up dataset to the Part1 gap it tested.
- Current edge states: F1-TYMP supported with heterogeneity; strict TYMP function partial; direct TYMP-MTA unsupported; MTA-SPP1 pending client; general TYMP-SPP1 unsupported; GSE175540 SPP1-vascular dataset-specific support; immune-exclusion and clinical translation unsupported.
- Do not claim a fully validated chain. The integrated diagram must colour every edge by evidence state.
- No Phase107/108 background job remains. Client validation data/figures remain outside project analysis.
"""
def read(p):
    b=p.read_bytes();return b.decode('utf-8-sig'),b.startswith(b'\xef\xbb\xbf')
def write(p,t,bom):p.write_text(t,encoding='utf-8-sig' if bom else 'utf-8',newline='')
def append(p,m,a):
    t,b=read(p);c=t.count(m)
    if c>1:raise RuntimeError(m)
    if c==0:
        n='\r\n' if '\r\n' in t else '\n';write(p,t.rstrip('\r\n')+a.replace('\n',n)+n,b)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for x in iter(lambda:f.read(1024*1024),b''):h.update(x)
    return h.hexdigest().upper()
append(STAGE/'PROJECT_STATUS.md','## Phase 107-108 evidence escalation and unified client architecture (2026-08-31)',STATUS)
append(STAGE/'DECISIONS.md','## 2026-08-31 - Phase 107-108 unified-chain decision',DECISIONS)
append(STAGE/'CODEX_HANDOFF.md','## Phase 107-108 current unified-chain state (2026-08-31)',HANDOFF)
p=STAGE/'RESULTS_MANIFEST.tsv';text,bom=read(p);r=csv.DictReader(text.splitlines(),delimiter='\t');header=r.fieldnames;rows=list(r)
if header is None:raise RuntimeError('header')
ids={x['result_id'] for x in rows};ts=datetime.now().astimezone().isoformat(timespec='seconds')
sources=[('phase107_evidence_escalation_audit',ROOT/'phase107_evidence_escalation_audit_local_v1',26),('phase108_unified_chain_architecture',ROOT/'phase108_unified_chain_architecture_local_v1',22)]
for stage,d,count in sources:
    for a in sorted(x for x in d.iterdir() if x.is_file()):
        stem=re.sub(r'[^a-z0-9]+','_',a.stem.lower()).strip('_');rid=f'{stage}_{stem}_v1'
        if rid in ids:continue
        kind={'.py':'formal_script','.tsv':'evidence_table_tsv','.json':'validation_json','.md':'audit_document_md'}.get(a.suffix,'audit_artifact')
        row={'result_id':rid,'stage':stage,'artifact_type':kind,'absolute_path':str(a),'generated_at':ts,'validation_status':'PASS','validation_evidence':f'SHA256 {sha(a)}; validation {count}/{count} PASS','notes':'Evidence-escalation and unified no-Part client architecture'}
        text=text.rstrip('\r\n')+'\r\n'+'\t'.join(row[c] for c in header)+'\r\n';ids.add(rid)
write(p,text,bom);print('PHASE107_108_STATE_STAGING_UPDATE_PASS')
