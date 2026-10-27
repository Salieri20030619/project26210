from pathlib import Path
import ast,json,hashlib,datetime,sys
import numpy as np,pandas as pd,scipy
from scipy.stats import rankdata,wilcoxon
from statsmodels.stats.multitest import multipletests
O=Path(__file__).resolve().parent;R=O.parent;S=R/'phase162_MTA_polyamine_part1_reassessment_20260915_v1';P174=R/'phase174_SMS_targeted_short_path_20260916_v1'
manifest=[]
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def track(p):
 p=Path(p);manifest.append(dict(path=str(p),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns,sha256=sha(p)));return p
def read(p,**kw):return pd.read_csv(track(p),sep='\t',**kw)
original=track(P174/'run_analysis.py')
tree=ast.parse(original.read_text());selected=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['rz','partial','fit','boot','RE']]
assert len(selected)==5
helper='import numpy as np\nfrom scipy.stats import rankdata,t,norm\n\n'+'\n\n'.join(ast.unparse(n) for n in selected)
(O/'stats_helpers.py').write_text(helper);sys.path.insert(0,str(O));from stats_helpers import rz,partial,fit,boot,RE
track(O/'PLAN_BEFORE_ANALYSIS_CN.md')
m=read(S/'PART1_TARGET_CELL_METADATA.tsv.gz')
cov=read(P174/'MALIGNANT_COVARIATES.tsv.gz').set_index('cell_id')
oldSMS=read(P174/'MALIGNANT_PATIENT_SMS_EFFECTS.tsv')
programs={'production_core':['ODC1','SRM','SMS'],'SAM_supply':['AMD1','MAT2A']}
models={'raw':[],'depth':['depth'],'depth_cycle_hypoxia':['depth','cycle','hypoxia'],'depth_genes_cycle_hypoxia':['depth','genes','cycle','hypoxia']}
elig=[];effects=[];drop=[];cache={};matched=[]
def z(v):
 v=np.asarray(v,float);return (v-v.mean(0))/v.std(0,ddof=1)
for patient,dm in m.groupby('curated_patient_id',sort=True):
 for module,genes in programs.items():
  detections={g:int((dm['count:'+g]>0).sum()) for g in ['TYMP']+genes}
  ok=len(dm)>=50 and min(detections.values())>=10
  elig.append(dict(patient=patient,study=dm.study.iloc[0],module=module,cells=len(dm),detections=json.dumps(detections),eligible=ok))
  if not ok:continue
  C0=cov.loc[dm.exact_cell_id]
  assert (C0.patient==patient).all()
  E=dm[['expr:'+g for g in genes]].to_numpy()
  Z=z(E);RK=np.column_stack([rz(E[:,j]) for j in range(E.shape[1])])
  scores={'gene_z_mean':Z.mean(1),'gene_rank_mean':RK.mean(1)}
  X=rz(dm['expr:TYMP'])
  for score,sv in scores.items():
   Y=rz(sv)
   for model,cn in models.items():
    C=np.column_stack([np.ones(len(dm))]+[rz(C0[x]) for x in cn]);DX=np.column_stack([C,X]);rank=np.linalg.matrix_rank(DX);cond=np.linalg.cond(DX);passed=rank==DX.shape[1] and cond<30
    rho=partial(X,Y,C) if passed else np.nan
    key='fit'+str(len(effects));cache[key+'_C']=C;cache[key+'_x']=X;cache[key+'_y']=Y
    effects.append(dict(patient=patient,study=dm.study.iloc[0],module=module,score=score,model=model,cells=len(dm),rho=rho,eligible=passed,rank=rank,columns=DX.shape[1],condition_number=cond,key=key))
    if module=='production_core' and score=='gene_z_mean':
     sm=oldSMS[(oldSMS.patient==patient)&(oldSMS.model==model)]
     assert len(sm)==1
     matched.append(dict(patient=patient,model=model,module_rho=rho,SMS_rho=float(sm.iloc[0].rho),difference=rho-float(sm.iloc[0].rho)))
     for j,g in enumerate(genes):
      y=rz(np.delete(Z,j,axis=1).mean(1));rr=partial(X,y,C)
      drop.append(dict(patient=patient,model=model,omitted_gene=g,rho=rr,full_module_rho=rho,note='descriptive_only_no_selection'))
D=pd.DataFrame(effects);D.to_csv(O/'MALIGNANT_PROGRAM_EFFECTS.tsv',sep='\t',index=False);pd.DataFrame(elig).to_csv(O/'MALIGNANT_ELIGIBILITY.tsv',sep='\t',index=False)
pd.DataFrame(drop).to_csv(O/'GENE_LEAVE_ONE_DIAGNOSTIC.tsv',sep='\t',index=False);pd.DataFrame(matched).to_csv(O/'MATCHED_CORE_VS_SMS.tsv',sep='\t',index=False);np.savez_compressed(O/'MALIGNANT_CACHE.npz',**cache)
rng=np.random.default_rng(262175);summary=[];subgroups=[]
for (module,score,model),d in D.groupby(['module','score','model']):
 v=d.loc[d.eligible,'rho'].to_numpy();n=len(v)
 bs=np.median(v[rng.integers(0,n,(5000,n))],axis=1)
 summary.append(dict(module=module,score=score,model=model,n=n,positive=int((v>0).sum()),negative=int((v<0).sum()),median_rho=float(np.median(v)),P=float(wilcoxon(v).pvalue),low=float(np.quantile(bs,.025)),high=float(np.quantile(bs,.975))))
 for study in sorted(d.study.unique()):
  for mode,dd in [('within',d[d.study==study]),('omit',d[d.study!=study])]:
   vv=dd.loc[dd.eligible,'rho'].to_numpy();subgroups.append(dict(module=module,score=score,model=model,study=study,subset=mode,n=len(vv),positive=int((vv>0).sum()),median_rho=float(np.median(vv)) if len(vv) else np.nan))
SS=pd.DataFrame(summary)
for key,ix in SS.groupby(['score','model']).groups.items():SS.loc[ix,'FDR2']=multipletests(SS.loc[ix,'P'],method='fdr_bh')[1]
SS.to_csv(O/'MALIGNANT_PROGRAM_SUMMARY.tsv',sep='\t',index=False);pd.DataFrame(subgroups).to_csv(O/'MALIGNANT_STUDY_SENSITIVITY.tsv',sep='\t',index=False)
print('MALIGNANT COMPLETE',flush=True)
adjs={'raw':[],'MTAP':['MTAP'],'MTAP_myeloid':['MTAP','MYELOID_PROXY'],'MTAP_myeloid_hypoxia':['MTAP','MYELOID_PROXY','HYPOXIA_PROXY']}
rows=[];draws={};values={}
for ci,cohort in enumerate(['RC18','RC20','CPTAC','CPTAC_val']):
 d=read(S/('PUBLIC_PATIENTS_'+cohort+'.tsv'),index_col=0);assert d.index.is_unique
 for mi,(module,genes) in enumerate(programs.items()):
  assert np.allclose(z(d[genes]).mean(1),d[module],atol=1e-12)
  for ai,(adj,cc) in enumerate(adjs.items()):
   dd=d[['TYMP',module,'MET:5-methylthioadenosine (MTA)']+cc].dropna();arr=dd.to_numpy();assert len(dd)==len(d)
   rr=fit(arr);bb,ii=boot(arr,262175+ci*100+mi*10+ai)
   key=cohort+'__'+module+'__'+adj;draws[key]=bb;values[key+'_values']=arr;values[key+'_indices']=ii[:20]
   rows.append(dict(cohort=cohort,module=module,adjustment=adj,**rr,ab_low=np.quantile(bb,.025),ab_high=np.quantile(bb,.975),ab_var=np.var(bb,ddof=1),bootstrap_positive_fraction=(bb>0).mean()))
 print('PUBLIC COMPLETE',cohort,flush=True)
F=pd.DataFrame(rows)
for key,ix in F.groupby(['cohort','adjustment']).groups.items():
 for edge in ['a','b']:F.loc[ix,edge+'_FDR2']=multipletests(F.loc[ix,edge+'_P'],method='fdr_bh')[1]
F.to_csv(O/'PUBLIC_PROGRAM_SHORT_COHORT.tsv',sep='\t',index=False)
meta=[];edges=[];loo=[]
for (module,adj),g in F.groupby(['module','adjustment']):
 rr,w=RE(g.ab,g.ab_var);bd=np.stack([draws[x.cohort+'__'+module+'__'+adj] for x in g.itertuples(index=False)]);pooled=w@bd;draws['POOLED__'+module+'__'+adj]=pooled
 meta.append(dict(module=module,adjustment=adj,n=int(g.n.sum()),positive_cohorts=int((g.ab>0).sum()),**rr,bootstrap_low=np.quantile(pooled,.025),bootstrap_high=np.quantile(pooled,.975),bootstrap_positive_fraction=(pooled>0).mean(),weights=json.dumps(dict(zip(g.cohort,w)))))
 for edge in ['a','b']:
  oo,ww=RE(np.arctanh(g[edge+'_rho']),1/(g.n-g[edge+'_k']-3))
  edges.append(dict(module=module,adjustment=adj,edge=edge,n=int(g.n.sum()),positive_cohorts=int((g[edge+'_rho']>0).sum()),rho=np.tanh(oo['estimate']),P=oo['P'],low=np.tanh(oo['low']),high=np.tanh(oo['high']),I2=oo['I2']))
 for omit in g.cohort:
  h=g[g.cohort!=omit];oo,ww=RE(h.ab,h.ab_var);loo.append(dict(module=module,adjustment=adj,omitted=omit,n=int(h.n.sum()),**oo))
MM=pd.DataFrame(meta);EE=pd.DataFrame(edges)
for key,ix in MM.groupby('adjustment').groups.items():MM.loc[ix,'FDR2']=multipletests(MM.loc[ix,'P'],method='fdr_bh')[1]
for key,ix in EE.groupby(['adjustment','edge']).groups.items():EE.loc[ix,'FDR2']=multipletests(EE.loc[ix,'P'],method='fdr_bh')[1]
legacy=read(S/'PUBLIC_META_AXIS_ASSOCIATIONS.tsv');legacy=legacy[legacy.gene.isin(programs)]
legacy.to_csv(O/'LEGACY_PROGRAM_ASSOCIATIONS.tsv',sep='\t',index=False)
H=read(R/'phase163_polyamine_MTA_serial_receiver_20260915_v1/META_CONDITIONAL_EDGES.tsv')
for i,r in EE.iterrows():
 if r.module=='production_core' and r.adjustment!='MTAP_myeloid_hypoxia':
  old=H[(H.candidate==r.module)&(H.adjustment==r.adjustment)&(H.edge==('a' if r.edge=='a' else 'd'))];assert len(old)==1
  assert abs(float(old.iloc[0].random_rho)-r.rho)<1e-10
  EE.loc[i,'historical_FDR3']=float(old.iloc[0].random_FDR)
 if r.edge=='b' and r.adjustment in ['raw','MTAP','MTAP_myeloid']:
  aj='myeloid' if r.adjustment=='MTAP_myeloid' else 'raw';ed='b' if r.adjustment=='raw' else 'b_adjust_MTAP'
  old=legacy[(legacy.gene==r.module)&(legacy.adjustment==aj)&(legacy.edge==ed)&legacy.outcome.str.contains('MTA')]
 elif r.edge=='a' and r.adjustment=='raw':
  old=legacy[(legacy.gene==r.module)&(legacy.adjustment=='raw')&(legacy.edge=='a')&legacy.outcome.str.contains('MTA')]
 else:old=pd.DataFrame()
 if len(old):
  assert len(old)==1 and abs(float(old.iloc[0].random_rho)-r.rho)<1e-10;EE.loc[i,'historical_FDR2']=float(old.iloc[0].random_FDR)
MM.to_csv(O/'PUBLIC_PROGRAM_SHORT_META.tsv',sep='\t',index=False);EE.to_csv(O/'PUBLIC_PROGRAM_EDGES_META.tsv',sep='\t',index=False);pd.DataFrame(loo).to_csv(O/'PUBLIC_LEAVE_ONE_COHORT.tsv',sep='\t',index=False)
np.savez_compressed(O/'PUBLIC_BOOTSTRAPS.npz',**draws);np.savez_compressed(O/'PUBLIC_INPUT_CACHE.npz',**values)
pd.DataFrame(manifest).to_csv(O/'INPUT_MANIFEST.tsv',sep='\t',index=False)
status=dict(stage='Phase175',status='COMPUTATION_COMPLETE_PENDING_VALIDATION',modules=programs,malignant_n={'production_core':14,'SAM_supply':16},public_n=255,public_models=len(F),bootstrap_per_model=5000,SPP1_in_model=False,MTAP_in_score=False,new_figures=0,python=sys.version,numpy=np.__version__,scipy=scipy.__version__,pandas=pd.__version__,completed_at=datetime.datetime.now().astimezone().isoformat())
(O/'COMPUTE_STATUS.json').write_text(json.dumps(status,indent=2))
print(SS.to_string(index=False));print(EE.to_string(index=False));print(MM.to_string(index=False))
