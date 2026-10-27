#!/usr/bin/env python3
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


def bh(values):
    p=np.asarray(values,float); out=np.full(len(p),np.nan); good=np.isfinite(p)
    if not good.any(): return out
    x=p[good]; order=np.argsort(x); ranked=x[order]*len(x)/np.arange(1,len(x)+1)
    ranked=np.minimum.accumulate(ranked[::-1])[::-1]; restored=np.empty(len(x)); restored[order]=np.minimum(ranked,1)
    out[np.where(good)[0]]=restored; return out


def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser(); p.add_argument('--training',required=True); p.add_argument('--out',required=True); a=p.parse_args()
    training=Path(a.training).resolve(); out=Path(a.out).resolve(); part=Path(str(out)+'.part')
    if out.exists() or part.exists(): raise RuntimeError('Refusing existing validation output')
    part.mkdir(parents=True)
    direct=pd.read_csv(training/'GSE249082_direct_MTA_edgeR_v1.tsv',sep='\t')
    ko=pd.read_csv(training/'GSE276571_MtapKO_edgeR_v1.tsv',sep='\t')
    mtapi=pd.read_csv(training/'GSE298469_MTAPi_kidney_edgeR_v1.tsv',sep='\t')
    atac=pd.read_csv(training/'GSE166213_MTA_CD8_ATAC_promoter_v1.tsv',sep='\t')
    cand=pd.read_csv(training/'phase77_direct_candidates_cross_context_audit_v1.tsv',sep='\t')
    sig=pd.read_csv(training/'phase77_weighted_signatures_v1.tsv',sep='\t')
    counts=pd.read_csv(training/'phase77_signature_counts_v1.tsv',sep='\t')
    qual=pd.read_csv(training/'phase77_consensus_qualification_v1.tsv',sep='\t')
    anchors=pd.read_csv(training/'phase77_macrophage_anchor_audit_v1.tsv',sep='\t')
    checks=[]
    def ck(name,condition,observed,expected): checks.append({'check':name,'status':'PASS' if condition else 'FAIL','observed':observed,'expected':expected})
    ck('direct_replicates',set(direct[['n_control','n_perturbed']].drop_duplicates().itertuples(index=False,name=None))=={(4,4)},str(direct[['n_control','n_perturbed']].drop_duplicates().values.tolist()),'[[4,4]]')
    ck('ko_replicates',set(ko[['n_control','n_perturbed']].drop_duplicates().itertuples(index=False,name=None))=={(3,3)},str(ko[['n_control','n_perturbed']].drop_duplicates().values.tolist()),'[[3,3]]')
    ck('mtapi_replicates',set(mtapi[['n_control','n_perturbed']].drop_duplicates().itertuples(index=False,name=None))=={(6,5)},str(mtapi[['n_control','n_perturbed']].drop_duplicates().values.tolist()),'[[6,5]]')
    qerr=max(float(np.nanmax(np.abs(bh(x.PValue)-x.FDR))) for x in [direct,ko,mtapi,atac])
    ck('bh_replay',qerr<1e-10,qerr,'<1e-10')
    expected=direct[(direct.FDR<.05)&(direct.logFC.abs()>=.5)].copy()
    ck('candidate_rule',set(expected.human_symbol)==set(cand.human_symbol),len(cand),len(expected))
    replay=(cand.validation_testable.astype(bool)&(cand.ko_same.astype(bool)|cand.mtapi_same.astype(bool))&~(cand.ko_opposite.astype(bool)|cand.mtapi_opposite.astype(bool)))
    ck('consensus_rule_replay',np.array_equal(replay.to_numpy(),cand.consensus_pass.astype(bool).to_numpy()),int(np.sum(replay.to_numpy()!=cand.consensus_pass.astype(bool).to_numpy())),0)
    direct_sig=sig[sig.signature_id=='DIRECT_MTA_GSE249082']; consensus=sig[sig.signature_id=='CONSENSUS_MTA_CROSS_CONTEXT']
    ck('direct_signature_subset',set(direct_sig.human_symbol).issubset(set(cand.human_symbol)),len(direct_sig),'<='+str(len(cand)))
    ck('consensus_signature_subset',set(consensus.human_symbol).issubset(set(cand.loc[cand.consensus_pass.astype(bool),'human_symbol'])),len(consensus),'<='+str(int(cand.consensus_pass.astype(bool).sum())))
    merged=sig.merge(direct[['human_symbol','logFC']],on='human_symbol',suffixes=('','_source'))
    weight_error=float(np.max(np.abs(merged.weight-np.minimum(np.abs(merged.logFC_source),3)))) if len(merged) else 0.0
    direction_error=int(np.sum(merged.direction!=np.where(merged.logFC_source>0,'UP','DOWN'))) if len(merged) else 0
    ck('weight_replay',weight_error<1e-12,weight_error,'<1e-12'); ck('direction_replay',direction_error==0,direction_error,0)
    max_per_direction=int(sig.groupby(['signature_id','direction']).size().max()) if len(sig) else 0
    ck('top100_cap',max_per_direction<=100,max_per_direction,'<=100')
    n_up=int((consensus.direction=='UP').sum()); n_down=int((consensus.direction=='DOWN').sum())
    expected_qual='QUALIFIED_PRIMARY_TRANSFER' if n_up>=10 and n_down>=10 else 'CONSENSUS_SIGNATURE_NOT_QUALIFIED'
    ck('qualification_replay',qual.qualification.iloc[0]==expected_qual,qual.qualification.iloc[0],expected_qual)
    ck('anchor_set',set(anchors.human_symbol)=={'ARG1','VEGFA','IL10','TIMP1'},';'.join(sorted(anchors.human_symbol)),'ARG1;IL10;TIMP1;VEGFA')
    ck('no_forced_anchor_insertion',set(sig.human_symbol).issubset(set(cand.human_symbol)),len(set(sig.human_symbol)-set(cand.human_symbol)),0)
    ck('unique_signature_gene_direction',not sig.duplicated(['signature_id','human_symbol']).any(),int(sig.duplicated(['signature_id','human_symbol']).sum()),0)
    frame=pd.DataFrame(checks); frame.to_csv(part/'validation_checks_v1.tsv',sep='\t',index=False)
    failed=int((frame.status=='FAIL').sum()); summary={'status':'PASS' if failed==0 else 'FAIL','checks':len(frame),'passed':int((frame.status=='PASS').sum()),'failed':failed,'qualification':expected_qual}
    (part/'validation_summary_v1.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    core=sorted(x for x in training.iterdir() if x.is_file())
    pd.DataFrame([{'file':x.name,'bytes':x.stat().st_size,'sha256':sha256(x)} for x in core]).to_csv(part/'CORE_OUTPUT_SHA256_v1.tsv',sep='\t',index=False)
    if failed: raise RuntimeError(summary)
    (part/'INDEPENDENT_VALIDATION_PASS_v1').write_text('status=PASS\n',encoding='utf-8')
    os.replace(part,out); print(json.dumps(summary))


if __name__=='__main__': main()
