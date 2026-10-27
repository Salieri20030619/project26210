#!/usr/bin/env python3
from __future__ import annotations
import csv,json,os
from pathlib import Path

ROOT=Path('/mnt/BioDisk/wangmingyuan/project26210_phase74_virtual_knockout_enrichment_v1')
RUN=ROOT/'nominal_target_sensitivity_v2';VALID=ROOT/'nominal_target_sensitivity_validation_v2'
TARGETS={'TYMP','MTAP','SPP1','STAT3'};DATASETS=('GSE207493','GSE264586','GSE294109')

def read(path):
    with path.open('r',encoding='utf-8',newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
def genes(s):return {x for x in (s or '').split(';') if x and x!='NA'}
checks=[]
def check(name,ok,detail):checks.append({'check':name,'status':'PASS' if ok else 'FAIL','detail':str(detail)})

def main():
    VALID.mkdir(parents=True,exist_ok=True)
    check('complete_marker',(RUN/'NOMINAL_SENSITIVITY_COMPLETE_V2').exists(),RUN/'NOMINAL_SENSITIVITY_COMPLETE_V2')
    gene=read(RUN/'phase74_gene_level_target_audit_v2.tsv');gcross=read(RUN/'phase74_gene_level_target_cross_cohort_v2.tsv')
    source=read(RUN/'phase74_nominal_enrichment_target_gene_rows_v2.tsv');cross=read(RUN/'phase74_nominal_enrichment_target_gene_cross_cohort_v2.tsv');counts=read(RUN/'phase74_nominal_target_gene_term_counts_v2.tsv')
    check('gene_rows',len(gene)==48,len(gene));check('gene_cross_rows',len(gcross)==16,len(gcross));check('source_rows',len(source)==587,len(source));check('count_rows',len(counts)==16,len(counts))
    check('gene_targets',set(r['queried_target_gene'] for r in gene)==TARGETS,TARGETS)
    check('source_has_target',all(genes(r['matched_target_genes'])&TARGETS for r in source),'every row contains target')
    idx={(r['dataset'],r['knockout_gene'],r['ontology'],r['term_id']):r for r in source}
    unique_ok=len(cross)==len({(r['knockout_gene'],r['ontology'],r['term_id']) for r in cross})
    check('cross_unique',unique_ok,len(cross))
    recompute=True;common_ok=True;presence_ok=True
    for r in cross:
        rr={d:idx.get((d,r['knockout_gene'],r['ontology'],r['term_id'])) for d in DATASETS}
        p1=rr['GSE264586'] is not None;p2=rr['GSE294109'] is not None
        common=(genes(rr['GSE264586']['matched_target_genes'])&genes(rr['GSE294109']['matched_target_genes'])) if p1 and p2 else set()
        expected='TWO_EXTERNAL_NOMINAL_SAME_TARGET' if p1 and p2 and common else ('TWO_EXTERNAL_NOMINAL_DIFFERENT_TARGET' if p1 and p2 else ('ONE_EXTERNAL_NOMINAL_TARGET' if p1^p2 else 'DISCOVERY_ONLY_NOMINAL_TARGET'))
        recompute &= r['external_nominal_decision']==expected
        common_ok &= genes(r['common_external_target_genes'])==common
        presence_ok &= (r['GSE264586_present']=='TRUE')==p1 and (r['GSE294109_present']=='TRUE')==p2 and (r['GSE207493_present']=='TRUE')==(rr['GSE207493'] is not None)
    check('cross_recompute',recompute,'all decisions reproduced');check('common_targets',common_ok,'all intersections reproduced');check('presence_flags',presence_ok,'all presence flags reproduced')
    n_two=sum(r['external_nominal_decision']=='TWO_EXTERNAL_NOMINAL_SAME_TARGET' for r in cross)
    check('count_reconcile',sum(int(r['two_external_same_target_terms']) for r in counts)==n_two,f'counts={sum(int(r["two_external_same_target_terms"]) for r in counts)} cross={n_two}')
    metadata=read(RUN/'run_metadata_v2.tsv')
    check('metadata',len(metadata)==1 and metadata[0]['status']=='PASS_PHASE74_NOMINAL_TARGET_SENSITIVITY_V2_VALIDATION_PENDING' and int(metadata[0]['two_external_same_target_terms'])==n_two,metadata)
    failures=[r for r in checks if r['status']=='FAIL']
    path=VALID/'validation_checks_v2.tsv';tmp=path.with_suffix('.tsv.part')
    with tmp.open('w',encoding='utf-8',newline='') as h:
        w=csv.DictWriter(h,fieldnames=['check','status','detail'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(checks)
    os.replace(tmp,path)
    summary={'status':'PASS' if not failures else 'FAIL','checks':len(checks),'passed':len(checks)-len(failures),'failed':len(failures),'v1_cross_status':'INVALID_PRESERVED'}
    (VALID/'validation_summary_v2.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    if failures:raise SystemExit('NOMINAL_TARGET_VALIDATION_FAIL '+json.dumps(failures))
    (VALID/'INDEPENDENT_VALIDATION_PASS_V2').write_text('PHASE74_NOMINAL_TARGET_SENSITIVITY_VALIDATION_PASS_V2\n',encoding='utf-8')
    print(json.dumps(summary))

if __name__=='__main__':main()
