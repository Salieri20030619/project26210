from pathlib import Path
import json,gzip,hashlib,datetime,sys,subprocess
import numpy as np,pandas as pd,scipy
from scipy.io import mmread
from scipy.stats import rankdata,wilcoxon,t,norm
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits
O=Path(__file__).resolve().parent;R=O.parent;S=R/'phase162_MTA_polyamine_part1_reassessment_20260915_v1';Q=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1')
manifest=[]
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def track(p):
 p=Path(p);manifest.append(dict(path=str(p),size=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns,sha256=sha(p)));return p
def read(p,**kw):return pd.read_csv(track(p),sep='\t',**kw)
def rz(x):
 x=rankdata(np.asarray(x,float));sd=x.std(ddof=1)
 return (x-x.mean())/sd if sd>0 else np.zeros_like(x)
def partial(x,y,C):
 rx=x-C@np.linalg.lstsq(C,x,rcond=None)[0];ry=y-C@np.linalg.lstsq(C,y,rcond=None)[0]
 return float(np.corrcoef(rx,ry)[0,1])
gate=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
(O/'MEMORY_GATE.log').write_text(gate.stdout+gate.stderr);assert gate.returncode==0
track(O/'PLAN_BEFORE_ANALYSIS_CN.md')
legacy=read(S/'PART1_PATIENT_TYMP_GENE_CORRELATIONS.tsv');legacy=legacy[legacy.gene=='SMS']
legacy_summary=read(S/'PART1_TYMP_GENE_CORRELATION_SUMMARY.tsv');legacy_summary=legacy_summary[legacy_summary.gene=='SMS']
meta=read(S/'PART1_TARGET_CELL_METADATA.tsv.gz')
inp=Q/'07_pyscenic_v1/input'
genes=read(inp/'genes_v1.tsv').gene.tolist();ids=read(inp/'cells_v1.tsv').exact_cell_id.tolist()
with gzip.open(track(inp/'raw_counts_genes_by_cells_v1.mtx.gz'),'rb') as f:
 with threadpool_limits(limits=1):mat=mmread(f).tocsr()
assert mat.shape==(len(genes),11950)
meta=meta.set_index('exact_cell_id').loc[ids].reset_index()
tot=np.asarray(mat.sum(axis=0)).ravel();ng=np.asarray((mat>0).sum(axis=0)).ravel()
raw=mat[[genes.index('TYMP'),genes.index('SMS')]].toarray().T
assert np.array_equal(raw[:,0],meta['count:TYMP']) and np.array_equal(raw[:,1],meta['count:SMS']) and np.array_equal(tot,meta.total_counts)
gmt=track(R/'part2_rebuild_v1/06_functional_programs/run_v1/resources/MSigDB_Hallmark_2020.gmt')
exclude=set('ODC1 SRM SMS MTAP AMD1 MAT2A TYMP SLC16A3'.split());gs={}
for line in gmt.read_text().splitlines():
 fields=line.split('\t')
 if fields[0] in ['G2-M Checkpoint','Hypoxia']:gs[fields[0]]=[g for g in fields[2:] if g in genes and g not in exclude]
assert set(gs)=={'G2-M Checkpoint','Hypoxia'} and all(len(x)>10 for x in gs.values())
(O/'COVARIATE_GENESETS.json').write_text(json.dumps(dict(gene_sets=gs,excluded=sorted(exclude),source=str(gmt)),indent=2))
expr={name:np.log1p(mat[[genes.index(g) for g in names]].toarray().T/tot[:,None]*1e4) for name,names in gs.items()}
del mat
patids=sorted(legacy.patient.unique());assert len(patids)==16
models={'raw':[],'depth':['depth'],'depth_cycle':['depth','cycle'],'depth_genes':['depth','genes'],'depth_hypoxia':['depth','hypoxia'],'depth_cycle_hypoxia':['depth','cycle','hypoxia'],'depth_genes_cycle_hypoxia':['depth','genes','cycle','hypoxia']}
rows=[];cache={};covrows=[]
for pi,pid in enumerate(patids):
 ix=np.flatnonzero(meta.curated_patient_id.to_numpy()==pid);dm=meta.iloc[ix]
 assert len(ix)>=50 and (raw[ix,0]>0).sum()>=10 and (raw[ix,1]>0).sum()>=10
 vals={'depth':np.log1p(tot[ix]),'genes':ng[ix]}
 for label,name in [('cycle','G2-M Checkpoint'),('hypoxia','Hypoxia')]:
  arr=expr[name][ix];sd=arr.std(0,ddof=1);keep=sd>0
  vals[label]=((arr[:,keep]-arr[:,keep].mean(0))/sd[keep]).mean(1)
 x=rz(dm['expr:TYMP']);y=rz(dm['expr:SMS'])
 for mi,(model,names) in enumerate(models.items()):
  C=np.column_stack([np.ones(len(ix))]+[rz(vals[n]) for n in names]);design=np.column_stack([C,x])
  rank=np.linalg.matrix_rank(design);cond=np.linalg.cond(design)
  ok=rank==design.shape[1] and cond<30
  rho=partial(x,y,C) if ok else np.nan
  key=f'p{pi}_m{mi}';cache[key+'_C']=C;cache[key+'_x']=x;cache[key+'_y']=y
  rows.append(dict(patient=pid,study=dm.study.iloc[0],model=model,cells=len(ix),TYMP_detected=int((raw[ix,0]>0).sum()),SMS_detected=int((raw[ix,1]>0).sum()),rank=rank,columns=design.shape[1],condition_number=cond,eligible=ok,rho=rho,cache_key=key))
 covrows.extend(dict(patient=pid,cell_id=ids[j],depth=vals['depth'][k],genes=vals['genes'][k],cycle=vals['cycle'][k],hypoxia=vals['hypoxia'][k]) for k,j in enumerate(ix))
D=pd.DataFrame(rows);D.to_csv(O/'MALIGNANT_PATIENT_SMS_EFFECTS.tsv',sep='\t',index=False)
pd.DataFrame(covrows).to_csv(O/'MALIGNANT_COVARIATES.tsv.gz',sep='\t',index=False)
np.savez_compressed(O/'MALIGNANT_DESIGN_CACHE.npz',**cache)
rng=np.random.default_rng(262174);summ=[];loo=[];groups=[]
for model,df in D.groupby('model'):
 v=df.loc[df.eligible,'rho'].to_numpy();n=len(v)
 P=float(wilcoxon(v).pvalue) if n>=10 else np.nan
 bs=np.median(v[rng.integers(0,n,(5000,n))],axis=1) if n>=10 else np.array([np.nan])
 old=legacy_summary[legacy_summary.adjustment==model]
 summ.append(dict(model=model,n_patients=n,positive=int((v>0).sum()),negative=int((v<0).sum()),median_rho=float(np.median(v)),P=P,median_low=float(np.nanquantile(bs,.025)),median_high=float(np.nanquantile(bs,.975)),historical_FDR_axis6=float(old.iloc[0].FDR_axis6) if len(old) else np.nan))
 for row in df[df.eligible].itertuples(index=False):
  vv=df.loc[df.eligible&(df.patient!=row.patient),'rho'].to_numpy();loo.append(dict(model=model,omitted_patient=row.patient,n=len(vv),median_rho=np.median(vv),positive=int((vv>0).sum())))
 for study in sorted(df.study.unique()):
  for kind,dd in [('within',df[df.study==study]),('omit',df[df.study!=study])]:
   vv=dd.loc[dd.eligible,'rho'].to_numpy()
   groups.append(dict(model=model,study=study,subset=kind,n=len(vv),positive=int((vv>0).sum()),median_rho=np.median(vv) if len(vv) else np.nan))
SS=pd.DataFrame(summ);SS['supplement_FDR_4']=np.nan
mask=~SS.model.isin(['raw','depth','depth_cycle']);SS.loc[mask,'supplement_FDR_4']=multipletests(SS.loc[mask,'P'].fillna(1),method='fdr_bh')[1]
SS.to_csv(O/'MALIGNANT_SMS_SUMMARY.tsv',sep='\t',index=False);pd.DataFrame(loo).to_csv(O/'MALIGNANT_LEAVE_ONE_PATIENT.tsv',sep='\t',index=False);pd.DataFrame(groups).to_csv(O/'MALIGNANT_STUDY_SENSITIVITY.tsv',sep='\t',index=False)
D0=D[D.model.isin(['raw','depth','depth_cycle'])].merge(legacy,left_on=['patient','model'],right_on=['patient','adjustment'],suffixes=('_new','_old'))
assert len(D0)==48 and np.allclose(D0.rho_new,D0.rho_old,atol=1e-9)
D0.to_csv(O/'LEGACY_SMS_REPRODUCTION.tsv',sep='\t',index=False)
print('MALIGNANT COMPLETE',flush=True)
# The shorter A-X-Y path excludes SPP1 altogether.
adjs={'raw':[],'MTAP':['MTAP'],'MTAP_myeloid':['MTAP','MYELOID_PROXY'],'MTAP_myeloid_hypoxia':['MTAP','MYELOID_PROXY','HYPOXIA_PROXY']}
def fit(arr):
 n=len(arr);A=rz(arr[:,0]);X=rz(arr[:,1]);Y=rz(arr[:,2]);C=np.column_stack([np.ones(n)]+[rz(arr[:,j]) for j in range(3,arr.shape[1])])
 DX=np.column_stack([C,A]);DY=np.column_stack([C,A,X])
 a=float(np.linalg.lstsq(DX,X,rcond=None)[0][-1]);yc=np.linalg.lstsq(DY,Y,rcond=None)[0];b=float(yc[-1]);cp=float(yc[-2]);total=float(np.linalg.lstsq(DX,Y,rcond=None)[0][-1])
 ar=partial(A,X,C);br=partial(X,Y,DX);k=C.shape[1]-1
 def pv(r,df):return float(2*t.sf(abs(r)*np.sqrt(df/max(1-r*r,1e-15)),df))
 return dict(n=n,a=a,b=b,ab=a*b,cprime=cp,total=total,identity_error=total-cp-a*b,a_rho=ar,a_P=pv(ar,n-k-2),a_k=k,b_rho=br,b_P=pv(br,n-k-3),b_k=k+1,condition_X=float(np.linalg.cond(DX)),condition_Y=float(np.linalg.cond(DY)))
def boot(arr,seed,B=5000):
 rng=np.random.default_rng(seed);ii=rng.integers(0,len(arr),(B,len(arr)))
 vr=[rankdata(arr[:,j][ii],axis=1) for j in range(arr.shape[1])]
 vr=[(v-v.mean(1)[:,None])/v.std(1,ddof=1)[:,None] for v in vr]
 A,X,Y=vr[:3];C=[np.ones_like(A)]+vr[3:]
 def reg(y,preds):
  d=np.stack(preds,axis=2);g=np.einsum('bni,bnj->bij',d,d);rhs=np.einsum('bni,bn->bi',d,y)
  return np.einsum('bij,bj->bi',np.linalg.pinv(g),rhs)
 aa=reg(X,C+[A])[:,-1];bb=reg(Y,C+[A,X])[:,-1]
 return aa*bb,ii
RR=[];draws={};data={};jk=[]
for ci,cohort in enumerate(['RC18','RC20','CPTAC','CPTAC_val']):
 f=read(S/('PUBLIC_PATIENTS_'+cohort+'.tsv'),index_col=0);assert f.index.is_unique
 for ai,(adj,covs) in enumerate(adjs.items()):
  dat=f[['TYMP','SMS','MET:5-methylthioadenosine (MTA)']+covs].dropna();arr=dat.to_numpy()
  assert len(arr)==len(f)
  fitr=fit(arr);bb,ii=boot(arr,262174+ci*100+ai)
  key=cohort+'__'+adj;draws[key]=bb;data[key+'_values']=arr;data[key+'_indices']=ii[:25];data[key+'_ids']=dat.index.to_numpy(str)
  assert np.isfinite(bb).all()
  RR.append(dict(cohort=cohort,adjustment=adj,**fitr,ab_low=float(np.quantile(bb,.025)),ab_high=float(np.quantile(bb,.975)),ab_bootstrap_variance=float(np.var(bb,ddof=1)),bootstrap_positive_fraction=float((bb>0).mean()),bootstrap_n=len(bb)))
  if adj=='MTAP_myeloid':
   for i,pid in enumerate(dat.index):
    v=fit(np.delete(arr,i,axis=0));jk.append(dict(cohort=cohort,omitted_patient=pid,n=v['n'],a=v['a'],b=v['b'],ab=v['ab']))
 print('PUBLIC COMPLETE',cohort,flush=True)
F=pd.DataFrame(RR);F['a_supplement_FDR2']=np.nan;F['b_supplement_FDR2']=np.nan
for i,row in F.iterrows():
 q=multipletests([row.a_P,row.b_P],method='fdr_bh')[1];F.loc[i,['a_supplement_FDR2','b_supplement_FDR2']]=q
# Preserve historical three-candidate FDR where the identical edges exist.
H=read(R/'phase163_polyamine_MTA_serial_receiver_20260915_v1/COHORT_CONDITIONAL_EDGES.tsv')
for i,row in F.iterrows():
 for edge,new in [('a','a'),('d','b')]:
  h=H[(H.cohort==row.cohort)&(H.candidate=='SMS')&(H.adjustment==row.adjustment)&(H.edge==edge)]
  if len(h):
   assert abs(float(h.iloc[0].rho)-row[new+'_rho'])<1e-9
   F.loc[i,new+'_historical_FDR3']=float(h.iloc[0].FDR)
F.to_csv(O/'PUBLIC_SHORT_PATH_COHORT.tsv',sep='\t',index=False)
np.savez_compressed(O/'PUBLIC_AB_BOOTSTRAPS.npz',**draws);np.savez_compressed(O/'PUBLIC_INPUT_CACHE.npz',**data)
pd.DataFrame(jk).to_csv(O/'PUBLIC_LEAVE_ONE_PATIENT.tsv',sep='\t',index=False)
def RE(est,var):
 w=1/np.asarray(var);est=np.asarray(est);mu=np.average(est,weights=w);Q=float(np.sum(w*(est-mu)**2));den=w.sum()-(w*w).sum()/w.sum();tau=max(0,(Q-len(est)+1)/den)
 rw=1/(np.asarray(var)+tau);m=float(np.average(est,weights=rw));se=float(np.sqrt(1/rw.sum()))
 return dict(estimate=m,SE=se,low=m-1.96*se,high=m+1.96*se,P=float(2*norm.sf(abs(m/se))),tau2=tau,I2=max(0,(Q-len(est)+1)/Q)*100 if Q>0 else 0),rw/rw.sum()
MM=[];EE=[];LO=[]
hm=read(R/'phase163_polyamine_MTA_serial_receiver_20260915_v1/META_CONDITIONAL_EDGES.tsv')
for adj,g in F.groupby('adjustment'):
 out,weights=RE(g.ab,g.ab_bootstrap_variance)
 bd=np.stack([draws[x.cohort+'__'+adj] for x in g.itertuples(index=False)])
 pooled=weights@bd;draws['POOLED__'+adj]=pooled
 MM.append(dict(adjustment=adj,cohorts=4,n=int(g.n.sum()),positive_cohorts=int((g.ab>0).sum()),**out,bootstrap_low=float(np.quantile(pooled,.025)),bootstrap_high=float(np.quantile(pooled,.975)),bootstrap_positive_fraction=float((pooled>0).mean()),weight_rule='fixed_point_estimate_RE_weights_for_bootstrap',weights=json.dumps(dict(zip(g.cohort,weights)))))
 for edge in ['a','b']:
  rr,ww=RE(np.arctanh(g[edge+'_rho']),1/(g.n-g[edge+'_k']-3))
  old=hm[(hm.candidate=='SMS')&(hm.adjustment==adj)&(hm.edge==('a' if edge=='a' else 'd'))]
  EE.append(dict(adjustment=adj,edge=edge,n=int(g.n.sum()),positive_cohorts=int((g[edge+'_rho']>0).sum()),rho=np.tanh(rr['estimate']),low=np.tanh(rr['low']),high=np.tanh(rr['high']),P=rr['P'],I2=rr['I2'],historical_FDR3=float(old.iloc[0].random_FDR) if len(old) else np.nan))
 for omit in g.cohort:
  gg=g[g.cohort!=omit];rr,ww=RE(gg.ab,gg.ab_bootstrap_variance)
  LO.append(dict(adjustment=adj,omitted_cohort=omit,n=int(gg.n.sum()),positive_cohorts=int((gg.ab>0).sum()),**rr))
M=pd.DataFrame(MM);M['sensitivity_FDR3']=np.nan;ix=M.adjustment!='MTAP_myeloid';M.loc[ix,'sensitivity_FDR3']=multipletests(M.loc[ix,'P'],method='fdr_bh')[1]
E=pd.DataFrame(EE)
for adj,ix in E.groupby('adjustment').groups.items():E.loc[ix,'supplement_FDR2']=multipletests(E.loc[ix,'P'],method='fdr_bh')[1]
M.to_csv(O/'PUBLIC_SHORT_PATH_META.tsv',sep='\t',index=False);E.to_csv(O/'PUBLIC_CONDITIONAL_EDGES_META.tsv',sep='\t',index=False);pd.DataFrame(LO).to_csv(O/'PUBLIC_LEAVE_ONE_COHORT.tsv',sep='\t',index=False)
np.savez_compressed(O/'PUBLIC_AB_BOOTSTRAPS.npz',**draws)
pd.DataFrame(manifest).to_csv(O/'INPUT_MANIFEST.tsv',sep='\t',index=False)
status=dict(stage='Phase174',status='COMPUTATION_COMPLETE_PENDING_VALIDATION',malignant_patients=16,malignant_models=len(models),public_cohorts=4,public_patients=255,short_models=len(F),bootstrap_per_model=5000,only_candidate='SMS',SPP1_in_models=False,new_figures=0,python=sys.version,numpy=np.__version__,scipy=scipy.__version__,pandas=pd.__version__,created_at=datetime.datetime.now().astimezone().isoformat())
(O/'COMPUTE_STATUS.json').write_text(json.dumps(status,indent=2))
print(SS.to_string(index=False));print(M.to_string(index=False));print(E.to_string(index=False));print(F[F.adjustment=='MTAP_myeloid'][['cohort','a','b','ab','ab_low','ab_high','a_rho','a_P','b_rho','b_P']].to_string(index=False));print(json.dumps(status))
