from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
from scipy.stats import rankdata,wilcoxon,t as tdist,norm
from statsmodels.stats.multitest import multipletests
root=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=root/'phase187_IL6_project_receiver_validation_20260917_v1';d.mkdir(exist_ok=True)
plan="""# 计算前固定方案
继承Phase186候选；用户授权无分歧连续推进。
冻结75库/细胞注释，不扩展成纤维标签。成纤维来源关联因只有1肿瘤患者且IL6未检出，登记不可估计，不重分群或改用周细胞冒充。
整体肿瘤成熟巨噬：每患者>=50细胞，IL6与SPP1各>=10检出细胞；主检验细胞内IL6—SPP1的深度调整秩相关，以患者效应推断，n>=10才Wilcoxon；单个预设主要终点无FDR。控制检出基因数/C4为敏感性，保留全部结果。原STAT3靶集中若含IL6，不使用原评分验证IL6—STAT3，不换评分追阳性。
公共实测MTA四队列255患者：A MTA—IL6调整MTAP+髓系；B IL6—SPP1同时调整MTA+MTAP+髓系。FDR在A/B两项家族内，分别报告队列与Fisher固定效应合并。额外缺氧调整为独立敏感性家族。无中介ab、无新签名/图。
"""
(d/'PLAN_BEFORE_ANALYSIS_CN.md').write_text(plan)
old=root/'phase182_IL10_LGALS3_receiver_audit_20260916_v1'
meta=pd.read_csv(old/'SOURCE_CELL_METADATA.tsv.gz',sep='\t',index_col=0);z=np.load(old/'SOURCE_PANEL_CACHE.npz')
assert np.array_equal(meta.index.to_numpy(str),z['cell_ids']);assert len(meta)==258120 and meta.sample_id.nunique()==75
raw=z['counts'];genes=list(z['genes'])
targets=pd.read_csv(old/'STAT3_TARGETS_WITHOUT_LGALS3.tsv',sep='\t')
circ={res:{'IL6_in_targets':bool((g.gene=='IL6').any()),'IL6ST_in_targets':bool((g.gene=='IL6ST').any())} for res,g in targets.groupby('resource')}
(d/'STAT3_CIRCULARITY_AUDIT.json').write_text(json.dumps(circ,indent=2))
ss=pd.read_csv(old/'SOURCE_PANEL_SUMMARY.tsv',sep='\t');ss[ss.gene=='IL6'].to_csv(d/'IL6_SOURCE_SUMMARY.tsv',sep='\t',index=False)
rc=pd.read_csv(old/'RECEIVER_COVERAGE_SUMMARY.tsv',sep='\t');rc[rc.feature.str.contains('IL6')].to_csv(d/'IL6_RECEIVER_COVERAGE.tsv',sep='\t',index=False)
mm=pd.read_csv(root/'part2_rebuild_v1/08_regulons/input_v1/CELL_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str})
mac=mm[(~mm.candidate.isin(['1','Monocyte']))&(mm.source_tissue.str.lower()=='tumor')]
def partial(a,b,cv):
 a=rankdata(np.round(a,12));b=rankdata(np.round(b,12))
 C=np.column_stack([np.ones(len(a))]+[rankdata(np.round(x,12)) for x in cv])
 ar=a-C@np.linalg.lstsq(C,a,rcond=None)[0];br=b-C@np.linalg.lstsq(C,b,rcond=None)[0]
 if np.std(ar)<1e-10 or np.std(br)<1e-10:return np.nan
 return float(np.corrcoef(ar,br)[0,1])
effects=[];elig=[]
for pat,m in mac.groupby('curated_patient_id'):
 idx=meta.index.get_indexer(m.index);assert (idx>=0).all()
 v=raw[idx];dep=m.total_counts_raw.to_numpy(float)
 il=v[:,genes.index('IL6')];sp=v[:,genes.index('SPP1')]
 ok=len(m)>=50 and sum(il>0)>=10 and sum(sp>0)>=10
 elig.append({'patient':pat,'study':m.study.iloc[0],'cells':len(m),'IL6_detected':sum(il>0),'SPP1_detected':sum(sp>0),'eligible':ok})
 if not ok:continue
 il=np.log1p(il/dep*1e4);sp=np.log1p(sp/dep*1e4)
 for adj,cv in [('depth',[np.log1p(dep)]),('depth_genes',[np.log1p(dep),m.n_genes_raw.to_numpy()]),('depth_C4',[np.log1p(dep),(m.candidate=='4').to_numpy(int)])]:
  effects.append({'patient':pat,'study':m.study.iloc[0],'adjustment':adj,'rho':partial(il,sp,cv)})
ef=pd.DataFrame(effects);ef.to_csv(d/'MACRO_PATIENT_EFFECTS.tsv',sep='\t',index=False);pd.DataFrame(elig).to_csv(d/'MACRO_ELIGIBILITY.tsv',sep='\t',index=False)
summ=[];leave=[]
for adj,g in ef.groupby('adjustment'):
 v=g.rho.dropna().to_numpy();rng=np.random.default_rng(187);bs=np.median(v[rng.integers(0,len(v),(5000,len(v)))],axis=1)
 summ.append({'adjustment':adj,'n_patients':len(v),'positive':sum(v>0),'median_rho':np.median(v),'CI_low':np.quantile(bs,.025),'CI_high':np.quantile(bs,.975),'P':wilcoxon(v,method='approx',correction=True).pvalue if len(v)>=10 else np.nan})
 for st in g.study.unique():
  a=g[g.study!=st].rho.dropna();leave.append({'adjustment':adj,'omit_study':st,'n':len(a),'median_rho':a.median(),'positive':sum(a>0)})
pd.DataFrame(summ).to_csv(d/'MACRO_SUMMARY.tsv',sep='\t',index=False);pd.DataFrame(leave).to_csv(d/'MACRO_LEAVE_STUDY.tsv',sep='\t',index=False)
pub=Path('/mnt/BioDisk/wangmingyuan/project26210_phase98_public_paired_metabolomics_v1');mapping=pd.read_csv(pub/'MasterMapping_ccRCC.csv',dtype=str).fillna('')
edges=[];mapchecks=[]
for cohort in ['RC18','RC20','CPTAC','CPTAC_val']:
 f=pd.read_csv(root/'phase179_STAT3_external_activity_validation_20260916_v1/run_v2'/f'BULK_PATIENTS_{cohort}_collectri_zmean.tsv',sep='\t',index_col=0)
 rp=pub/'downloads_v1'/f'matched_tpm_{cohort}.csv'
 first=pd.read_csv(rp,nrows=0).columns[0]
 rr=pd.read_csv(rp,usecols=[first,'IL6','SPP1','TYMP']).set_index(first)
 mp=mapping[(mapping.Dataset==cohort)&(mapping.TN=='Tumor')].copy()
 mp['patient']=mp.CommonID.where(mp.CommonID!='',mp.RNAID.str.replace(r'-T$','',regex=True)).str.replace(r'^RC(?:18|20):','',regex=True)
 e=np.log2(rr.loc[mp.RNAID]+1);e.index=mp.patient;e=e.loc[f.index]
 assert e.index.is_unique and f.index.is_unique
 assert np.allclose(e.SPP1,f.SPP1) and np.allclose(e.TYMP,f.TYMP)
 f['IL6']=e.IL6
 f.to_csv(d/f'PUBLIC_PATIENTS_{cohort}.tsv',sep='\t')
 mapchecks.append({'cohort':cohort,'n':len(f),'SPP1_TYMP_mapping_reproduced':True})
 for adj,covs in [('main',['MTAP','MYELOID_PROXY']),('hypoxia',['MTAP','MYELOID_PROXY','HYPOXIA_PROXY'])]:
  for edge,x,y,cc in [('MTA_IL6','MET:5-methylthioadenosine (MTA)','IL6',covs),('IL6_SPP1','IL6','SPP1',covs+['MET:5-methylthioadenosine (MTA)'])]:
   t=f[[x,y]+cc].dropna();rho=partial(t[x].to_numpy(),t[y].to_numpy(),[t[c].to_numpy() for c in cc]);n=len(t);k=len(cc);df=n-k-2
   p=2*tdist.sf(abs(rho*np.sqrt(df/(1-rho*rho))),df)
   edges.append({'cohort':cohort,'adjustment':adj,'edge':edge,'n':n,'k':k,'rho':rho,'P':p})
ed=pd.DataFrame(edges);ed['FDR']=ed.groupby(['cohort','adjustment']).P.transform(lambda p:multipletests(p,method='fdr_bh')[1]);ed.to_csv(d/'PUBLIC_COHORT_EDGES.tsv',sep='\t',index=False)
def pool(g):
 w=g.n-g.k-3;zz=np.arctanh(g.rho);zzm=np.sum(w*zz)/sum(w);se=1/np.sqrt(sum(w))
 return {'rho':np.tanh(zzm),'CI_low':np.tanh(zzm-1.96*se),'CI_high':np.tanh(zzm+1.96*se),'P':2*norm.sf(abs(zzm/se)),'n':int(g.n.sum()),'positive_cohorts':int(sum(g.rho>0)),'cohorts':len(g)}
po=[];lo=[]
for (adj,edge),g in ed.groupby(['adjustment','edge']):
 po.append({'adjustment':adj,'edge':edge,**pool(g)})
 for cohort in g.cohort:lo.append({'adjustment':adj,'edge':edge,'omit_cohort':cohort,**pool(g[g.cohort!=cohort])})
po=pd.DataFrame(po);po['FDR']=po.groupby('adjustment').P.transform(lambda p:multipletests(p,method='fdr_bh')[1]);po.to_csv(d/'PUBLIC_META.tsv',sep='\t',index=False);pd.DataFrame(lo).to_csv(d/'PUBLIC_LEAVE_COHORT.tsv',sep='\t',index=False)
pd.DataFrame(mapchecks).to_csv(d/'PUBLIC_MAPPING_VALIDATION.tsv',sep='\t',index=False)
assert sum(x['n'] for x in mapchecks)==255
print('CIRC',circ);print('MACRO',pd.DataFrame(summ).to_string(index=False));print('PUBLIC',po.to_string(index=False))
(d/'COMPUTE_STATUS.json').write_text(json.dumps({'status':'computed','patients_bulk':255,'macro_patients':len(elig),'source_fibroblast_tumor_patients':1,'cell_alignment_pass':True,'old_models_modified':False},indent=2))
