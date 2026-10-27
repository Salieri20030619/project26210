#!/usr/bin/env python3
import argparse, json, math, os
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats

OUTCOMES={
 'R01_DIRECT_MTA_to_SPP1_STATE':('SPP1_STATE',['GPNMB','LGALS3','FABP5','CD9','LPL'],['Epithelial','Myeloid'],1),
 'R02_DIRECT_MTA_to_VASCULAR_NO_VEGFA':('VASCULAR_NO_VEGFA',['MMP9','MMP14','ANGPTL4','IL1B','LGALS3'],['Epithelial','Myeloid','Vascular'],1),
}
ABUND={'Epithelial':['Epithelial cell'],'Myeloid':['Macrophage','Monocyte'],'Vascular':['Endothelial','Pericyte/SMC']}

def bh(p):
 p=np.asarray(p,float);o=np.argsort(p);x=p[o]*len(p)/np.arange(1,len(p)+1);x=np.minimum.accumulate(x[::-1])[::-1];r=np.empty(len(p));r[o]=np.minimum(x,1);return r
def prep(s):return 'FFPE' if int(s.replace('GSM',''))<=5924041 else 'FROZEN'
def dense_gene(a,g):
 names=np.asarray([str(x).upper() for x in a.var_names]);idx=np.where(names==g)[0]
 if len(idx)!=1: raise RuntimeError(f'{g}: expected one row, found {len(idx)}')
 x=a.X[:,int(idx[0])];return np.asarray(x.toarray() if sparse.issparse(x) else x,float).reshape(-1)
def weighted_score(a,sig,order,exclude):
 names=set(str(x).upper() for x in a.var_names);s=sig[~sig.human_symbol.isin(exclude)].copy();s=s[s.human_symbol.isin(names)]
 up=int((s.direction=='UP').sum());down=int((s.direction=='DOWN').sum());coverage=len(s)/max(1,len(sig[~sig.human_symbol.isin(exclude)]))
 if coverage<.70 or up<8 or down<8: raise RuntimeError(f'signature coverage failed {coverage} {up} {down}')
 m=np.column_stack([dense_gene(a,g)[order] for g in s.human_symbol]);mean=m.mean(axis=0);sd=m.std(axis=0);sd[sd==0]=1;z=(m-mean)/sd
 sign=np.where(s.direction.to_numpy()=='UP',1.,-1.);w=s.weight.to_numpy(float);return (z@(w*sign))/w.sum(),coverage,up,down,';'.join(s.human_symbol)
def program(a,genes,order):
 m=np.column_stack([dense_gene(a,g)[order] for g in genes]);sd=m.std(axis=0);sd[sd==0]=1;return ((m-m.mean(axis=0))/sd).mean(axis=1)
def abundance(d,name):
 cols=[f'mean__{x}' for x in ABUND[name]]
 if any(x not in d for x in cols):raise RuntimeError(f'missing abundance {name}')
 return d[cols].sum(axis=1).to_numpy(float)
def residual(v,covs):
 y=stats.rankdata(v);x=np.column_stack([np.ones(len(v))]+[stats.rankdata(c) for c in covs]);return y-x@np.linalg.lstsq(x,y,rcond=None)[0]
def partial(x,y,covs):
 m=np.column_stack([x,y]+covs);m=m[np.all(np.isfinite(m),axis=1)];rx=residual(m[:,0],[m[:,i] for i in range(2,m.shape[1])]);ry=residual(m[:,1],[m[:,i] for i in range(2,m.shape[1])]);return float(np.corrcoef(rx,ry)[0,1]),len(m)

def main():
 p=argparse.ArgumentParser();p.add_argument('--phase54',required=True);p.add_argument('--phase55',required=True);p.add_argument('--signature',required=True);p.add_argument('--out',required=True);a=p.parse_args()
 p54=Path(a.phase54);p55=Path(a.phase55);out=Path(a.out);part=Path(str(out)+'.part')
 if out.exists() or part.exists():raise RuntimeError('refusing existing output')
 part.mkdir(parents=True);sig=pd.read_csv(a.signature,sep='\t');sig=sig[sig.signature_id=='DIRECT_MTA_GSE249082'].copy()
 if len(sig)!=89:raise RuntimeError('unexpected signature')
 ids={x.name for x in (p54/'run_v1/sections_v1').iterdir() if x.is_dir()};ids55={x.name for x in (p55/'run_v1/sections_v1').iterdir() if x.is_dir()}
 if ids!=ids55 or len(ids)!=24:raise RuntimeError('section mismatch')
 rows=[];qc=[]
 for sid in sorted(ids):
  asec=ad.read_h5ad(p54/'run_v1/sections_v1'/sid/'processed_section_v1.h5ad');score=pd.read_csv(p54/'run_v1/sections_v1'/sid/'spot_metadata_scores_v1.tsv.gz',sep='\t');abu=pd.read_csv(p55/'run_v1/sections_v1'/sid/'spot_cell_abundance_v1.tsv.gz',sep='\t');score=score[score.qc_pass.astype(str).str.lower().isin(['true','1'])]
  d=score.merge(abu,on='barcode',validate='one_to_one');
  if set(d.barcode)!=set(asec.obs_names) or set(d.barcode)!=set(abu.barcode):raise RuntimeError(f'{sid}: barcode mismatch')
  order=pd.Index(asec.obs_names).get_indexer(d.barcode);vals={'SPP1_STATE':d['score_primary__spp1_state_no_spp1'].to_numpy(float),'VASCULAR_NO_VEGFA':program(asec,OUTCOMES['R02_DIRECT_MTA_to_VASCULAR_NO_VEGFA'][1],order)}
  logumi=np.log1p(d.total_counts.to_numpy(float));ab={k:abundance(d,k) for k in ABUND};prep_id=prep(sid)
  for rid,(outcome,exclude,covnames,direction) in OUTCOMES.items():
   x,cov,up,down,present=weighted_score(asec,sig,order,exclude);rho,n=partial(x,vals[outcome],[ab[c] for c in covnames]+[logumi]);rows.append({'sample_id':sid,'preparation':prep_id,'relation_id':rid,'outcome':outcome,'expected_direction':direction,'n_spots':n,'rho':rho,'coverage':cov,'present_up':up,'present_down':down,'evidence_use':'EXPLORATORY_DIRECT_SIGNATURE_ONLY'})
  qc.append({'sample_id':sid,'preparation':prep_id,'n_spots':len(d),'barcode_match':'PASS'});del asec
 rel=pd.DataFrame(rows);meta=[]
 for keys,g in rel.groupby(['preparation','relation_id'],sort=True):
  values=g.rho.to_numpy(float);test=stats.wilcoxon(values,alternative='two-sided');meta.append({'preparation':keys[0],'relation_id':keys[1],'n_sections':len(values),'median_rho':float(np.median(values)),'min_rho':float(values.min()),'max_rho':float(values.max()),'expected_direction_fraction':float(np.mean(values>0)),'p_value':float(test.pvalue)})
 meta=pd.DataFrame(meta);meta['q_value']=np.nan
 for _,idx in meta.groupby('preparation').groups.items():meta.loc[list(idx),'q_value']=bh(meta.loc[list(idx),'p_value'])
 meta['statistical_status']=np.where((meta.q_value<.05)&(meta.median_rho>0)&(meta.expected_direction_fraction>=8/12),'SUPPORTED','NOT_SUPPORTED');meta['evidence_use']='EXPLORATORY_DIRECT_SIGNATURE_ONLY'
 decisions=[]
 for rid,g in meta.groupby('relation_id'):
  s=dict(zip(g.preparation,g.statistical_status));n=sum(s.get(x)=='SUPPORTED' for x in ['FFPE','FROZEN']);decisions.append({'relation_id':rid,'ffpe_status':s.get('FFPE','MISSING'),'frozen_status':s.get('FROZEN','MISSING'),'preparations_supported':n,'statistical_decision':'TWO_PREPARATION_SUPPORTED' if n==2 else ('ONE_PREPARATION_SUPPORTED' if n==1 else 'NOT_SUPPORTED'),'evidence_use':'EXPLORATORY_ONLY_CONSENSUS_NOT_QUALIFIED'})
 pd.DataFrame(qc).to_csv(part/'section_qc_v1.tsv',sep='\t',index=False);rel.to_csv(part/'relations_by_section_v1.tsv',sep='\t',index=False);meta.to_csv(part/'relation_meta_v1.tsv',sep='\t',index=False);pd.DataFrame(decisions).to_csv(part/'decisions_v1.tsv',sep='\t',index=False)
 (part/'run_metadata_v1.json').write_text(json.dumps({'status':'PASS','sections':24,'spots':int(pd.DataFrame(qc).n_spots.sum()),'signature':'DIRECT_MTA_GSE249082','evidence_use':'EXPLORATORY_ONLY_CONSENSUS_NOT_QUALIFIED'},indent=2)+'\n')
 os.replace(part,out);print(json.dumps({'status':'PASS_PHASE77_SPATIAL_TRANSFER','sections':24}))
if __name__=='__main__':main()
