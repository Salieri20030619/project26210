from pathlib import Path
import pandas as pd,numpy as np,json
from scipy.stats import rankdata,wilcoxon,t as td,norm
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase187_IL6_project_receiver_validation_20260917_v1';old=r/'phase182_IL10_LGALS3_receiver_audit_20260916_v1'
audit=json.loads((d/'STAT3_CIRCULARITY_AUDIT.json').read_text());assert not any(v['IL6_in_targets'] for v in audit.values())
(d/'STAT3_SUPPLEMENT_PLAN.md').write_text('IL6不在冻结两套STAT3靶集中，补做IL6—STAT3代理；collectri zmean为主，dorothea zmean敏感性。沿同19患者资格及深度主调整、nGenes/C4敏感性。此为独立预设STAT3联系终点，不改变已计算MTA_IL6/IL6_SPP1家族。公共255例调整MTA+MTAP+髓系，缺氧敏感性。IL6ST在两靶集中，禁止用其与评分相关作为受体证据。')
meta=pd.read_csv(old/'SOURCE_CELL_METADATA.tsv.gz',sep='\t',index_col=0);z=np.load(old/'SOURCE_PANEL_CACHE.npz');counts=z['counts'];genes=list(z['genes'])
m=pd.read_csv(r/'part2_rebuild_v1/08_regulons/input_v1/CELL_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str});m=m[(~m.candidate.isin(['1','Monocyte']))&(m.source_tissue.str.lower()=='tumor')]
eligible=pd.read_csv(d/'MACRO_ELIGIBILITY.tsv',sep='\t');pats=set(eligible[eligible.eligible].patient);cache=np.load(old/'SCORE_EDGE_CACHE.npz')
def pc(x,y,cv):
 C=np.column_stack([np.ones(len(x))]+[rankdata(np.round(v,12)) for v in cv]);a=rankdata(np.round(x,12));b=rankdata(np.round(y,12))
 a-=C@np.linalg.lstsq(C,a,rcond=None)[0];b-=C@np.linalg.lstsq(C,b,rcond=None)[0]
 return np.corrcoef(a,b)[0,1]
out=[];checks=[]
for pat,g in m.groupby('curated_patient_id'):
 if pat not in pats:continue
 idx=meta.index.get_indexer(g.index);raw=counts[idx];dep=g.total_counts_raw.to_numpy()
 il=np.log1p(raw[:,genes.index('IL6')]/dep*1e4);sp=np.log1p(raw[:,genes.index('SPP1')]/dep*1e4)
 for res in ['collectri','dorothea']:
  a=cache[res+'__zmean__'+pat];assert a.shape[0]==len(g) and np.allclose(a[:,2],sp);checks.append(True)
  for adj,cv in [('depth',[a[:,3]]),('depth_genes',[a[:,3],a[:,4]]),('depth_C4',[a[:,3],a[:,5]])]:
   out.append({'patient':pat,'study':g.study.iloc[0],'resource':res,'adjustment':adj,'rho':pc(il,a[:,0],cv)})
out=pd.DataFrame(out);out.to_csv(d/'IL6_STAT3_MACRO_PATIENTS.tsv',sep='\t',index=False);summ=[]
for (res,adj),g in out.groupby(['resource','adjustment']):
 v=g.rho.to_numpy();summ.append({'resource':res,'adjustment':adj,'n':len(v),'positive':sum(v>0),'median_rho':np.median(v),'P':wilcoxon(v,method='approx',correction=True).pvalue})
pd.DataFrame(summ).to_csv(d/'IL6_STAT3_MACRO_SUMMARY.tsv',sep='\t',index=False)
es=[]
for cohort in ['RC18','RC20','CPTAC','CPTAC_val']:
 f=pd.read_csv(d/f'PUBLIC_PATIENTS_{cohort}.tsv',sep='\t',index_col=0)
 for res in ['collectri','dorothea']:
  frozen=pd.read_csv(r/'phase179_STAT3_external_activity_validation_20260916_v1/run_v2'/f'BULK_PATIENTS_{cohort}_{res}_zmean.tsv',sep='\t',index_col=0)
  for adj,cv in [('main',['MTAP','MYELOID_PROXY','MET:5-methylthioadenosine (MTA)']),('hypoxia',['MTAP','MYELOID_PROXY','MET:5-methylthioadenosine (MTA)','HYPOXIA_PROXY'])]:
   rho=pc(f.IL6.to_numpy(),frozen.loc[f.index,'activity'].to_numpy(),[f[c].to_numpy() for c in cv]);n=len(f);k=len(cv)
   es.append({'cohort':cohort,'resource':res,'adjustment':adj,'n':n,'k':k,'rho':rho,'P':2*td.sf(abs(rho*np.sqrt((n-k-2)/(1-rho*rho))),n-k-2)})
es=pd.DataFrame(es);es.to_csv(d/'IL6_STAT3_PUBLIC_COHORT.tsv',sep='\t',index=False);po=[]
for (res,adj),g in es.groupby(['resource','adjustment']):
 w=g.n-g.k-3;z=np.sum(w*np.arctanh(g.rho))/sum(w);se=1/np.sqrt(sum(w))
 po.append({'resource':res,'adjustment':adj,'n':int(g.n.sum()),'rho':np.tanh(z),'CI_low':np.tanh(z-1.96*se),'CI_high':np.tanh(z+1.96*se),'P':2*norm.sf(abs(z/se)),'positive_cohorts':sum(g.rho>0)})
pd.DataFrame(po).to_csv(d/'IL6_STAT3_PUBLIC_META.tsv',sep='\t',index=False)
(d/'STAT3_ALIGNMENT_VALIDATION.json').write_text(json.dumps({'patient_score_alignment_checks':len(checks),'all_pass':all(checks),'IL6_excluded_from_both_targets':True},indent=2))
print('MACRO',pd.DataFrame(summ).to_string(index=False));print('BULK',pd.DataFrame(po).to_string(index=False))
