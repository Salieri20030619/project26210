#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path('/mnt/BioDisk/wangmingyuan/project26210_phase74_virtual_knockout_enrichment_v1')
RUN = ROOT / 'enrichment_run_v1'
RESULTS = ROOT / 'results_v1'
DATASETS = ('GSE207493','GSE264586','GSE294109')
TARGETS = ('STAT3','SPP1','TYMP','MTAP')
ONTOLOGIES = ('GO_BP','GO_CC','GO_MF','KEGG_LEGACY','KEGG_MEDICUS')


def read_tsv(path: Path) -> list[dict[str,str]]:
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path,'rt',encoding='utf-8',newline='') as handle:
        return list(csv.DictReader(handle,delimiter='\t'))


def write_tsv(path: Path, rows: list[dict[str,object]], fields: list[str]) -> None:
    tmp = path.with_suffix(path.suffix+'.part')
    with tmp.open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    os.replace(tmp,path)


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''): h.update(block)
    return h.hexdigest()


def main() -> None:
    RESULTS.mkdir(parents=True,exist_ok=True)
    if not (RUN/'ENRICHMENT_COMPLETE').exists(): raise RuntimeError('Phase74 enrichment incomplete')
    supported=read_tsv(RUN/'phase74_enrichment_supported_v1.tsv')
    counts=read_tsv(RUN/'phase74_affected_gene_counts_v1.tsv')
    index={(r['dataset'],r['knockout_gene'],r['ontology'],r['term_id']):r for r in supported}
    terms=sorted({(r['knockout_gene'],r['ontology'],r['term_id'],r['term_name']) for r in supported})
    cross=[]
    for target,ontology,term_id,term_name in terms:
        ds={d:index.get((d,target,ontology,term_id)) for d in DATASETS}
        ext=[d for d in ('GSE264586','GSE294109') if ds[d] is not None]
        if len(ext)==2: decision='TWO_EXTERNAL_ENRICHED'
        elif len(ext)==1: decision='ONE_EXTERNAL_ENRICHED'
        else: decision='DISCOVERY_ONLY_ENRICHED'
        row={'knockout_gene':target,'ontology':ontology,'term_id':term_id,'term_name':term_name,
             'external_decision':decision,'external_supported_cohorts':';'.join(ext),
             'discovery_supported':'TRUE' if ds['GSE207493'] else 'FALSE','direction_available':'FALSE'}
        for d in DATASETS:
            row[f'{d}_q_value']=ds[d]['q_value'] if ds[d] else ''
            row[f'{d}_affected_gene_count']=ds[d]['affected_gene_count'] if ds[d] else ''
            row[f'{d}_affected_genes']=ds[d]['affected_genes'] if ds[d] else ''
        cross.append(row)
    cross.sort(key=lambda r:(TARGETS.index(r['knockout_gene']),ONTOLOGIES.index(r['ontology']),
                             {'TWO_EXTERNAL_ENRICHED':0,'ONE_EXTERNAL_ENRICHED':1,'DISCOVERY_ONLY_ENRICHED':2}[r['external_decision']],r['term_id']))
    fields=list(cross[0]) if cross else ['knockout_gene','ontology','term_id','term_name','external_decision']
    write_tsv(RESULTS/'phase74_cross_cohort_enrichment_v1.tsv',cross,fields)

    two=[r for r in cross if r['external_decision']=='TWO_EXTERNAL_ENRICHED']
    one=[r for r in cross if r['external_decision']=='ONE_EXTERNAL_ENRICHED']
    top=[]
    for target in TARGETS:
        pool=[r for r in two if r['knockout_gene']==target]
        pool.sort(key=lambda r:(max(float(r['GSE264586_q_value']),float(r['GSE294109_q_value'])),-int(r['GSE264586_affected_gene_count'])-int(r['GSE294109_affected_gene_count']),r['term_id']))
        for rank,r in enumerate(pool[:20],1):
            z=dict(r);z['within_knockout_rank']=rank;top.append(z)
    top_fields=['within_knockout_rank']+fields
    write_tsv(RESULTS/'phase74_top_two_external_terms_v1.tsv',top,top_fields)

    decision_counts=Counter(r['external_decision'] for r in cross)
    summary=['# Phase74 virtual-knockout differential regulation and enrichment summary v1','',
             '## Differential network regulation','']
    for r in counts:
        summary.append(f"- {r['dataset']} {r['knockout_gene']} KO: {r['affected_fdr05']}/{r['modeled_genes']} genes affected at gene-level FDR<0.05.")
    summary.extend(['','## Cross-cohort enrichment decisions','',
                    f"- TWO_EXTERNAL_ENRICHED: {decision_counts.get('TWO_EXTERNAL_ENRICHED',0)} terms.",
                    f"- ONE_EXTERNAL_ENRICHED: {decision_counts.get('ONE_EXTERNAL_ENRICHED',0)} terms.",
                    f"- DISCOVERY_ONLY_ENRICHED: {decision_counts.get('DISCOVERY_ONLY_ENRICHED',0)} terms.",'',
                    '## Top two-external terms by knockout',''])
    for target in TARGETS:
        selected=[r for r in top if r['knockout_gene']==target][:10]
        summary.append(f'### {target} KO')
        summary.append('')
        if selected:
            summary.extend([f"- {r['ontology']} | {r['term_name']} | q(ext1)={r['GSE264586_q_value']} | q(ext2)={r['GSE294109_q_value']}" for r in selected])
        else:
            summary.append('- No exact term achieved two-external enrichment.')
        summary.append('')
    summary.extend(['## Interpretation boundary','',
                    'The input is unsigned scTenifoldKnk network perturbation, not a post-KO expression matrix. Enrichment describes over-representation among network-affected genes and cannot be called upregulated or downregulated.',''])
    (RESULTS/'PHASE74_ENRICHMENT_SUMMARY_v1.md').write_text('\n'.join(summary),encoding='utf-8')

    metadata={'status':'PASS_PHASE74_FINALIZATION_VALIDATION_PENDING','completed_at':datetime.now(timezone.utc).astimezone().isoformat(),
              'supported_input_rows':len(supported),'cross_terms':len(cross),'two_external_terms':len(two),
              'one_external_terms':len(one),'direction_available':False}
    (RESULTS/'run_metadata_v1.json').write_text(json.dumps(metadata,indent=2)+'\n',encoding='utf-8')
    manifest=[]
    for path in sorted(RESULTS.iterdir()):
        if path.is_file() and path.name!='FILE_MANIFEST_SHA256_v1.tsv':
            manifest.append({'path':str(path.relative_to(ROOT)),'bytes':path.stat().st_size,'sha256':sha256(path)})
    write_tsv(RESULTS/'FILE_MANIFEST_SHA256_v1.tsv',manifest,['path','bytes','sha256'])
    print(json.dumps(metadata))


if __name__=='__main__': main()
