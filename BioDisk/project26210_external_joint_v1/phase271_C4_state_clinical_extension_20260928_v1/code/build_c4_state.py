from pathlib import Path
import numpy as np,pandas as pd,json,hashlib,datetime,warnings
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
p=Path(__file__).resolve().parents[1];r=p.parent
src=r/'phase200_union_C4_C10_C15_reassessment_20260918_v1/core/contrasts/tumor_C4.npz'
z=np.load(src,allow_pickle=True);genes=z['genes'];studies=z['study'];patients=z['patients']
exclude={'SPP1','TYMP','FN1','CLEC5A','OLR1','IL1RN'}
def select(train,tag):
 e=z['effects'][train];a=z['target_log2CPM'][train];b=z['reference_log2CPM'][train];det=z['target_detection'][train];st=studies[train];n=len(e)
 evaluable=((a>=1)|(b>=1)).sum(0)>=max(3,int(np.ceil(.2*n)))
 evaluable &= np.any(e!=0,axis=0)
 pv=np.full(len(genes),np.nan)
 with warnings.catch_warnings():
  warnings.simplefilter('ignore')
  pv[evaluable]=wilcoxon(e[:,evaluable],axis=0,zero_method='wilcox',correction=True,method='approx').pvalue
 q=np.full(len(genes),np.nan);q[evaluable]=multipletests(pv[evaluable],method='fdr_bh')[1]
 se=np.array([np.median(e[st==s],axis=0) for s in sorted(set(st))])
 t=pd.DataFrame({'gene':genes,'patients':n,'studies':len(set(st)),'median_delta':np.median(e,axis=0),'positive_fraction':(e>0).mean(0),'positive_study_fraction':(se>0).mean(0),'median_study_delta':np.median(se,axis=0),'median_detection':np.median(det,axis=0),'P':pv,'FDR':q,'tested':evaluable})
 t['excluded']=t.gene.isin(exclude)|t.gene.str.startswith(('MT-','RPL','RPS'))
 t['eligible']=(~t.excluded)&(t.P<.05)&(t.median_delta>=.5)&(t.positive_fraction>=.7)&(t.positive_study_fraction>=.8)&(t.median_detection>=.1)
 t['selected']=False
 ordered=t[t.eligible].sort_values(['median_study_delta','median_delta','gene'],ascending=[False,False,True]).head(20)
 t.loc[ordered.index,'selected']=True
 t.to_csv(p/f'tables/MARKERS_{tag}.tsv.gz',sep='\t',index=False)
 return ordered.gene.tolist(),t
panel,tab=select(np.ones(len(patients),bool),'full')
assert len(panel)>=16,('insufficient fixed markers',len(panel))
tab[tab.selected].sort_values(['median_study_delta','median_delta','gene'],ascending=[False,False,True]).to_csv(p/'tables/FROZEN_C4_SIGNATURE.tsv',sep='\t',index=False)
ref=r/'phase219_spatial_reference_identifiability_20260920_v1/reference'
rm=pd.read_csv(ref/'PATIENT_STATE_METADATA.tsv',sep='\t');rx=np.load(ref/'PATIENT_STATE_COUNTS.npz')
assert np.array_equal(genes,rx['genes'])
rc=rx['counts'];rl=np.log2(rc/rc.sum(1)[:,None]*1e6+1)
gmap={g:i for i,g in enumerate(genes)};pred=[];foldpan=[];lin=[]
for st in sorted(set(studies)):
 tr=studies!=st;test=~tr
 gs,t=select(tr,st);assert len(gs)>=16,(st,len(gs))
 ids=[gmap[g] for g in gs]
 for g in gs:foldpan.append({'heldout_study':st,'gene':g})
 for i in np.flatnonzero(test):
  v=z['target_log2CPM'][i,ids].mean();w=z['reference_log2CPM'][i,ids].mean()
  pred.append({'patient':patients[i],'study':st,'gene_n':len(gs),'C4_score':v,'other_macro_score':w,'delta':v-w,'positive':v>w})
 for i in rm.index[rm.study==st]:
  lin.append({'patient':rm.at[i,'patient'],'study':st,'label':rm.at[i,'label'],'cells':int(rm.at[i,'cells']),'score':float(rl[i,ids].mean())})
pr=pd.DataFrame(pred);ls=pd.DataFrame(lin);fp=pd.DataFrame(foldpan)
pr.to_csv(p/'tables/LOSO_PAIRED_STATE_VALIDATION.tsv',sep='\t',index=False)
ls.to_csv(p/'tables/LOSO_LINEAGE_SPECIFICITY.tsv',sep='\t',index=False)
fp.to_csv(p/'tables/LOSO_TRAINING_GENES.tsv',sep='\t',index=False)
ss=pr.groupby('study').agg(n=('patient','nunique'),positive=('positive','sum'),median_delta=('delta','median')).reset_index()
ss.to_csv(p/'tables/LOSO_STUDY_SUMMARY.tsv',sep='\t',index=False)
summ=ls.groupby(['study','label']).agg(patients=('patient','nunique'),median_score=('score','median')).reset_index()
summ.to_csv(p/'tables/LOSO_LINEAGE_SUMMARY.tsv',sep='\t',index=False)
status={'created':datetime.datetime.now().astimezone().isoformat(),'genes':panel,'gene_n':len(panel),'source_patients':len(patients),'source_studies':len(set(studies)),'LOSO_positive':int(pr.positive.sum()),'LOSO_n':len(pr),'LOSO_median_delta':float(pr.delta.median()),'all_studies_positive_median':bool((ss.median_delta>0).all()),'state_contrast_gate':bool(pr.positive.mean()>=.8 and (ss.median_delta>0).all()),'no_clinical_outcomes_loaded':True,'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'fixed_members_sha256':hashlib.sha256((p/'tables/FROZEN_C4_SIGNATURE.tsv').read_bytes()).hexdigest(),'interpretation':'C4-associated transcriptomic program, not cell fraction. Cross-lineage specificity described separately.'}
(p/'qc/SIGNATURE_FREEZE_BEFORE_CLINICAL.json').write_text(json.dumps(status,indent=2))
print(json.dumps(status,indent=2));print(ss.to_string(index=False));print(summ.pivot(index='label',columns='study',values='median_score').round(2).to_string())
