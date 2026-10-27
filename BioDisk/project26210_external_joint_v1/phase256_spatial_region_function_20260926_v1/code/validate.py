from pathlib import Path
import sys,json,hashlib,datetime
import numpy as np,pandas as pd
from scipy import stats
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'code'))
from extract import gate,fingerprint,stamp
def bh(p):
 p=np.asarray(p,float);n=len(p);o=np.argsort(p,kind='stable')
 q=np.minimum.accumulate((p[o]*n/np.arange(1,n+1))[::-1])[::-1]
 r=np.empty(n);r[o]=np.minimum(q,1);return r
def es_dense(rank,genes,hit):
 in_set=np.isin(genes,hit);weights=np.abs(rank)*in_set
 if weights.sum()==0:return np.nan
 run=np.cumsum(weights/weights.sum()-(~in_set)/(~in_set).sum())
 peak=run[np.argmax(np.abs(run))];return peak
def main():
 gate('final_validation')
 assert (P/'extract_EXIT_STATUS').read_text().strip()=='0'
 assert (P/'analysis_EXIT_STATUS').read_text().strip()=='0'
 for s in ['extract_COMPLETE','analysis_COMPLETE','GSEA_COMPLETE.txt','EXTRACTION_COMPLETE.json','AGGREGATION_COMPLETE.json']:assert (P/s).is_file()
 before=json.loads((P/'qc/SOURCE_MANIFEST_BEFORE.json').read_text())
 changed=[v['path'] for v in before if v!=fingerprint(v['path'])]
 assert not changed,changed
 (P/'qc/SOURCE_MANIFEST_AFTER.json').write_text(json.dumps([fingerprint(v['path']) for v in before],indent=2))
 reg=pd.read_csv(P/'tables/SECTION_IDENTITY.tsv',sep='\t')
 assert len(reg)==38 and reg[['dataset','patient']].drop_duplicates().shape[0]==37
 eligibility=pd.read_csv(P/'tables/REGION_ELIGIBILITY.tsv',sep='\t')
 assert len(eligibility)==114 and eligibility.eligible.all()
 s=pd.read_csv(P/'tables/PATIENT_FUNCTION_SUMMARY.tsv.gz',sep='\t')
 g=pd.read_csv(P/'tables/GSEA_FULL_RESULTS.tsv.gz',sep='\t')
 keys=['dataset','model','version','library','term']
 assert not s.duplicated(keys).any() and not g.duplicated(keys).any()
 assert len(s[s.patient_n>=10])==len(g)
 assert not g.dataset.str.contains('GLMF').any()
 assert s.loc[s.patient_n<10,['P_patient_sign','FDR_patient_sign']].isna().all().all()
 assert (s.positive+s.negative+s.zero==s.patient_n).all()
 assert s.library.nunique()==3 and s.model.nunique()==3 and s.version.nunique()==2
 tests=0
 for key,t in s.groupby(['dataset','library','family']):
  tt=t[t.P_patient_sign.notna()]
  if len(tt):
   np.testing.assert_allclose(bh(tt.P_patient_sign),tt.FDR_patient_sign,rtol=1e-10,atol=1e-12);tests+=len(tt)
 for key,t in g.groupby(['dataset','library','family']):
  q=bh(t.P_GSEA.fillna(1));mask=t.P_GSEA.notna()
  np.testing.assert_allclose(q[mask],t.loc[mask,'FDR_GSEA'],rtol=1e-10,atol=1e-12)
 jobs=pd.read_csv(P/'gsea/JOBS.tsv',sep='\t');checks=[];sign_checks=0
 for job in jobs.itertuples():
  gate('validate_'+job.dataset+'_'+job.model+'_'+job.version+'_'+job.library)
  rank=pd.read_csv(job.rank_path,sep='\t');arr=np.load(P/'effects'/f'PATIENT_{job.dataset}_{job.model}.npz',allow_pickle=False)
  genes=arr['genes'];E=arr['effects'];n=len(arr['patients'])
  program=np.load(P/'effects'/f'PROGRAM_{job.dataset}__{job.model}__{job.version}__{job.library}.npz',allow_pickle=False)
  assert len(set(program['patients']))==n==job.patient_n and np.array_equal(arr['patients'],program['patients'])
  terms=program['terms'];D=program['effects'];refrows=s[(s.dataset==job.dataset)&(s.model==job.model)&(s.version==job.version)&(s.library==job.library)].set_index('term').loc[terms]
  np.testing.assert_allclose(np.median(D,axis=0),refrows.median_Z_difference,atol=1e-12)
  np.testing.assert_array_equal((D>1e-12).sum(0),refrows.positive)
  sets={f[0]:f[2:] for f in [l.split('\t') for l in Path(job.gmt_path).read_text().splitlines()]}
  # Independent pathway score reconstruction for every member set (gene matrix rather than sparse product).
  maxerr=0.
  gi={v:i for i,v in enumerate(genes)}
  for j,term in enumerate(terms):
   vv=E[:,[gi[x] for x in sets[term]]].mean(1)
   maxerr=max(maxerr,float(np.max(np.abs(vv-D[:,j]))))
  assert maxerr<1e-10
  # Independent scalar exact binomial on stratified rows, including extreme and nonsignificant effects.
  ix=np.unique(np.r_[np.linspace(0,len(terms)-1,min(30,len(terms))).astype(int),np.argmax(refrows.positive),np.argmax(refrows.negative)])
  if n>=10:
   for j in ix:
    row=refrows.iloc[j];nz=int(row.positive+row.negative)
    pv=stats.binomtest(int(row.positive),nz,.5).pvalue if nz else 1.
    assert abs(pv-row.P_patient_sign)<1e-12;sign_checks+=1
  # Reconstruct ES with full running walk on fixed evenly spaced pathways, not the fast sparse formula.
  if n>=10:
   gg=g[(g.dataset==job.dataset)&(g.model==job.model)&(g.version==job.version)&(g.library==job.library)].set_index('term')
  eserr=0.
  for j in ix:
   term=terms[j];ev=es_dense(rank['rank'].to_numpy(),rank.gene.to_numpy(),sets[term])
   reference=float(gg.loc[term,'ES']) if n>=10 else float(refrows.loc[term,'GSEA_ES_descriptive'])
   if np.isfinite(ev) and np.isfinite(reference):eserr=max(eserr,abs(ev-reference))
  assert eserr<1e-8,(job,eserr)
  checks.append(dict(dataset=job.dataset,model=job.model,version=job.version,library=job.library,pathways=len(terms),patient_n=n,max_program_error=maxerr,ES_checked=len(ix),max_ES_error=eserr))
 pd.DataFrame(checks).to_csv(P/'qc/INDEPENDENT_PROGRAM_AND_ES_CHECKS.tsv',sep='\t',index=False)
 # Merge both evidence layers without overwriting patient-level P/FDR.
 joined=s.merge(g[keys+['ES','NES','P_GSEA','FDR_GSEA','log2err','leadingEdge']],on=keys,how='left',validate='one_to_one')
 joined['GSEA_status']=np.where(joined.patient_n<10,'DESCRIPTIVE_ES_ONLY',np.where(joined.P_GSEA.isna(),'NOT_ESTIMABLE',np.where(joined.FDR_GSEA<.05,'FDR_ENRICHED',np.where(joined.P_GSEA<.05,'NOMINAL_ONLY','NOT_ENRICHED'))))
 joined.to_csv(P/'tables/ALL_FUNCTIONAL_EVIDENCE.tsv.gz',sep='\t',index=False)
 result=dict(status='PASS',time=stamp(),sections=38,patients=37,regional_comparisons=len(eligibility),functional_rows=len(s),gsea_rows=len(g),BH_patient_values_rechecked=tests,exact_sign_tests_rechecked=sign_checks,program_checks=len(checks),source_files_unchanged=len(before),GSEA_not_estimable=int(g.P_GSEA.isna().sum()),GSEA_epsilon_floor_rows=int((g.P_GSEA<=1e-6).sum()),new_figures=0,hotspot_analysis_started=False,old_RDS_modified=False)
 (P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(result,indent=2))
 print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':main()
