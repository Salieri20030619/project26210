from pathlib import Path
import json,hashlib,subprocess
import numpy as np,pandas as pd
from scipy.stats import rankdata,t,wilcoxon
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');R=P/'part2_rebuild_v1';O=P/'phase157_ordered_bridge_reassessment_20260915_v1';V=O/'05_validation';V.mkdir(exist_ok=False)
with (V/'memory.log').open('w') as f:subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],stdout=f,stderr=subprocess.STDOUT,check=True)
checks=[]
def ck(name,value):checks.append({'check':name,'passed':bool(value)})
def bh(v):
 v=np.asarray(v,float);ok=np.isfinite(v);out=np.full(len(v),np.nan);a=v[ok];idx=np.argsort(a);q=np.minimum.accumulate((a[idx]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1];back=np.empty(len(a));back[idx]=np.minimum(q,1);out[ok]=back;return out
def rz(x):
 x=rankdata(x);return (x-x.mean())/x.std(ddof=1)
def partial(x,y,C):
 Q=C@np.linalg.pinv(C);a=x-Q@x;b=y-Q@y;r=np.corrcoef(a,b)[0,1];df=len(x)-np.linalg.matrix_rank(C)-1
 return r,2*t.sf(abs(r)*np.sqrt(df/max(1-r*r,1e-15)),df)
base=pd.read_csv(P/'phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv',sep='\t')
ck('16_unique_patients',len(base)==16 and base.patient.is_unique and base.sample_id.is_unique)
identity=pd.read_csv(O/'02_inputs/CELL_IDENTITY.tsv.gz',sep='\t',dtype={'role':str})
frozen=set(pd.read_csv(R/'09_downstream_crosswalk/C4_CANDIDATE_CELL_IDS.tsv.gz',sep='\t').cell_id)
ck('frozen_C4_exact',set(identity.loc[identity.role=='C4','cell_id'])==frozen)
ck('no_cell_duplicates',identity.cell_id.is_unique)
mal=pd.read_csv(O/'02_inputs/MATCHED_MALIGNANT_LOG2CPM.tsv.gz',sep='\t',index_col=0).loc[base.patient]
ck('TYMP_equals_previous_direct_test',np.allclose(mal.TYMP,base.TYMP_pseudobulk))
sig=pd.read_csv(O/'02_inputs/MTA_SIGNATURE_GENE_REGISTRY.tsv',sep='\t')
ck('signature_coverage',len(sig)==84 and (sig.direction=='UP').sum()==68 and (sig.direction=='DOWN').sum()==16)
ck('signature_no_endpoint_leakage',not set(sig.human_symbol)&{'TYMP','SPP1'})
cs=pd.read_csv(O/'02_inputs/CELL_MTA_RESPONSE_SCORES.tsv.gz',sep='\t')
cs['compartment']=np.where(cs.role=='Malignant','Malignant','Macrophage')
for key,g in cs.groupby(['study','compartment']):ck('signature_centering_'+str(key),abs(g.MTA_response_transcript_score.mean())<1e-10)
pairs=pd.read_csv(O/'03_MTA_to_new_C4/C4_MTA_RESPONSE_PATIENT_CONTRASTS.tsv',sep='\t');ss=pd.read_csv(O/'03_MTA_to_new_C4/C4_MTA_RESPONSE_SUMMARY.tsv',sep='\t')
for _,row in ss.iterrows():
 d=cs[cs.role!='Malignant']
 if row.scope=='tumor_only':d=d[d.source_tissue.str.lower()=='tumor']
 agg=d.groupby(['curated_patient_id','role']).agg(n=('cell_id','size'),score=('MTA_response_transcript_score','mean'))
 dif=[]
 for patient in d.curated_patient_id.unique():
  if (patient,'C4') not in agg.index or (patient,'OtherMacrophage') not in agg.index:continue
  a=agg.loc[(patient,'C4')];b=agg.loc[(patient,'OtherMacrophage')]
  if min(a.n,b.n)>=20:dif.append(a.score-b.score)
 dif=np.asarray(dif)
 ck('MTA_state_'+row.scope,len(dif)==row.n_patients and np.isclose(np.median(dif),row.median_effect) and (dif>0).sum()==row.positive and np.isclose(wilcoxon(dif,method='approx',correction=True).pvalue,row.P))
ck('MTA_state_BH',np.allclose(bh(ss.P),ss.FDR_two_scopes))
other=pd.read_csv(O/'04_other_X_screen/ALL_CANDIDATE_BRIDGE_RESULTS.tsv',sep='\t')
X=pd.read_csv(O/'04_other_X_screen/PATIENT_CANDIDATE_VALUES.tsv',sep='\t',index_col=0).loc[base.patient]
mta=pd.read_csv(O/'03_MTA_to_new_C4/MTA_PROXY_C4_BRIDGE_RESULTS.tsv',sep='\t')
mtax=pd.read_csv(O/'03_MTA_to_new_C4/MATCHED_MTA_PROXY_PATIENTS.tsv',sep='\t').set_index('patient').loc[base.patient]
for label,df,values,col in [('MTA',mta,mtax,'metric'),('other',other,X,'candidate')]:
 for i,row in df.iterrows():
  A=rz(base.TYMP_pseudobulk);M=rz(values[row[col]]);Y=rz(base.C4_fraction)
  C=np.ones((16,1)) if row.adjustment=='raw' else pd.get_dummies(base.study,dtype=float).to_numpy()
  ar,ap=partial(A,M,C);yr,yp=partial(M,Y,C);br,bp=partial(M,Y,np.column_stack([C,A]))
  designA=np.column_stack([C,A]);designB=np.column_stack([C,A,M])
  ca=(np.linalg.pinv(designA)@M)[-1];cb=np.linalg.pinv(designB)@Y;ct=(np.linalg.pinv(designA)@Y)[-1]
  ck(label+str(i)+'_design_rank',np.linalg.matrix_rank(designB)==designB.shape[1])
  ck(label+str(i)+'_correlations',np.allclose([ar,yr,br],[row.a_rho,row.XY_rho,row.b_partial_rho],atol=1e-10))
  ck(label+str(i)+'_P',np.allclose([ap,yp,bp],[row.a_P,row.XY_P,row.b_P],atol=1e-10))
  ck(label+str(i)+'_indirect',np.allclose([ca,cb[-1],ca*cb[-1],cb[-2],ct],[row.a_beta,row.b_beta,row.indirect_ab,row.direct_cprime,row.total_c],atol=1e-10))
  ck(label+str(i)+'_decomposition',abs(ct-cb[-2]-ca*cb[-1])<1e-10)
 if label=='MTA':
  for key,g in df.groupby('adjustment'):
   for p in ['a_P','XY_P','b_P']:ck(label+str(key)+p+'_BH',np.allclose(bh(g[p]),g[p+'_FDR_MTA_four']))
 else:
  for key,g in df.groupby(['family','adjustment']):
   for p in ['a_P','XY_P','b_P']:ck(label+str(key)+p+'_BH',np.allclose(bh(g[p]),g[p+'_FDR_family']))
reg=pd.read_csv(O/'04_other_X_screen/CANDIDATE_REGISTRY.tsv',sep='\t')
ck('116_candidates',len(reg)==116 and reg.candidate.is_unique)
ck('all_X_exclude_A_Y',all(not set(x.split(';'))&{'TYMP','SPP1'} for x in reg.genes))
ck('no_new_figures',not any(p.suffix.lower() in ['.png','.jpg','.svg','.pdf'] for p in O.rglob('*')))
assert all(x['passed'] for x in checks),[x for x in checks if not x['passed']]
(V/'CHECKS.json').write_text(json.dumps(checks,indent=2))
status={'status':'PASS','checks':len(checks),'passed':len(checks),'limitations':'Bootstrap intervals are conditional fixed-rank exploratory intervals; no causal or multiplicity-adjusted mediation declaration.'}
(V/'STATUS.json').write_text(json.dumps(status,indent=2))
sources=[
 P/'phase35_part1_c01_c04_cell_communication_v1/resource_lock_v1/cellchat_interactions_v1.tsv',
 R/'06_functional_programs/run_v1/resources/MSigDB_Hallmark_2020.gmt',
 R/'08_regulons/input_v1/CELL_METADATA.tsv.gz',
 Path('/mnt/BioDisk/wangmingyuan/project26210_phase77_mta_perturbation_signature_v1/training_v1/phase77_weighted_signatures_v1.tsv'),
 Path('/mnt/BioDisk/wangmingyuan/project26210_phase98_public_paired_metabolomics_v1/analysis_v2/cohort_associations_v2.tsv'),
 Path('/mnt/BioDisk/wangmingyuan/project26210_phase98_public_paired_metabolomics_v1/analysis_v2/cross_cohort_meta_analysis_v2.tsv'),
 P/'phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv']
manifest=[dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sources]
(O/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(status))
