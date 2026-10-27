#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
from pathlib import Path

ROOT=Path('/mnt/BioDisk/wangmingyuan/project26210_phase74_virtual_knockout_enrichment_v1')
RUN=ROOT/'enrichment_run_v1';RESULTS=ROOT/'results_v1';VALID=ROOT/'independent_validation_v1'
DATASETS=('GSE207493','GSE264586','GSE294109');TARGETS=('STAT3','SPP1','TYMP','MTAP')

def read(path):
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',encoding='utf-8',newline='') as h:return list(csv.DictReader(h,delimiter='\t'))

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

checks=[]
def check(name,ok,detail):checks.append({'check':name,'status':'PASS' if ok else 'FAIL','detail':str(detail)})

def main():
    VALID.mkdir(parents=True,exist_ok=True)
    check('enrichment_complete',(RUN/'ENRICHMENT_COMPLETE').exists(),RUN/'ENRICHMENT_COMPLETE')
    counts=read(RUN/'phase74_affected_gene_counts_v1.tsv')
    check('count_rows',len(counts)==12,len(counts))
    check('count_keys',{(r['dataset'],r['knockout_gene']) for r in counts}=={(d,t) for d in DATASETS for t in TARGETS},'12 unique dataset-KO keys')
    count_numeric=all(10<=int(r['affected_fdr05'])<int(r['modeled_genes']) and 0<float(r['symbol_mapping_fraction'])<=1 for r in counts)
    check('count_numeric',count_numeric,'affected and mapping values valid')
    affected=read(RUN/'phase74_affected_genes_all_v1.tsv.gz')
    check('affected_total',len(affected)==sum(int(r['affected_fdr05']) for r in counts),len(affected))
    check('affected_fdr',all(float(r['p.adj'])<0.05 for r in affected),'all retained genes pass gene-level FDR')
    enrich=read(RUN/'phase74_enrichment_all_v1.tsv.gz')
    supported=read(RUN/'phase74_enrichment_supported_v1.tsv')
    check('enrichment_nonempty',len(enrich)>0,len(enrich))
    check('five_ontologies',set(r['ontology'] for r in enrich)=={'GO_BP','GO_CC','GO_MF','KEGG_LEGACY','KEGG_MEDICUS'},set(r['ontology'] for r in enrich))
    check('supported_rule',all(float(r['q_value'])<0.05 and int(r['affected_gene_count'])>=3 and r['supported']=='TRUE' for r in supported),len(supported))
    cross=read(RESULTS/'phase74_cross_cohort_enrichment_v1.tsv')
    check('cross_unique',len(cross)==len({(r['knockout_gene'],r['ontology'],r['term_id']) for r in cross}),len(cross))
    allowed={'TWO_EXTERNAL_ENRICHED','ONE_EXTERNAL_ENRICHED','DISCOVERY_ONLY_ENRICHED'}
    check('cross_decisions',all(r['external_decision'] in allowed for r in cross),'allowed only')
    check('cross_direction',all(r['direction_available']=='FALSE' for r in cross),'unsigned recorded')
    for r in cross:
        n=sum(bool(r[f'{d}_q_value']) for d in ('GSE264586','GSE294109'))
        expected='TWO_EXTERNAL_ENRICHED' if n==2 else ('ONE_EXTERNAL_ENRICHED' if n==1 else 'DISCOVERY_ONLY_ENRICHED')
        if r['external_decision']!=expected:
            check('cross_recompute',False,f"{r['knockout_gene']} {r['term_id']}");break
    else:check('cross_recompute',True,'all decisions reproduced')
    top=read(RESULTS/'phase74_top_two_external_terms_v1.tsv')
    check('top_only_two_external',all(r['external_decision']=='TWO_EXTERNAL_ENRICHED' for r in top),len(top))
    metadata=json.loads((RESULTS/'run_metadata_v1.json').read_text(encoding='utf-8'))
    check('metadata_status',metadata.get('status')=='PASS_PHASE74_FINALIZATION_VALIDATION_PENDING',metadata.get('status'))
    check('metadata_dimensions',metadata.get('cross_terms')==len(cross),metadata)
    manifest=read(RESULTS/'FILE_MANIFEST_SHA256_v1.tsv')
    hash_ok=all((ROOT/r['path']).exists() and (ROOT/r['path']).stat().st_size==int(r['bytes']) and sha(ROOT/r['path'])==r['sha256'] for r in manifest)
    check('manifest_hashes',hash_ok and len(manifest)>=5,len(manifest))
    failures=[r for r in checks if r['status']=='FAIL']
    path=VALID/'validation_checks_v1.tsv';tmp=path.with_suffix('.tsv.part')
    with tmp.open('w',encoding='utf-8',newline='') as h:
        w=csv.DictWriter(h,fieldnames=['check','status','detail'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(checks)
    os.replace(tmp,path)
    summary={'status':'PASS' if not failures else 'FAIL','checks':len(checks),'passed':len(checks)-len(failures),'failed':len(failures)}
    (VALID/'validation_summary_v1.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    if failures:raise SystemExit('PHASE74_VALIDATION_FAIL '+json.dumps(failures))
    (VALID/'INDEPENDENT_VALIDATION_PASS_v1').write_text('PHASE74_INDEPENDENT_VALIDATION_PASS\n',encoding='utf-8')
    print(json.dumps(summary))

if __name__=='__main__':main()
