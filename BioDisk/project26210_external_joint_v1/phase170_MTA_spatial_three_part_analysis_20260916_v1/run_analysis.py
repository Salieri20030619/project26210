from pathlib import Path
import sys,json,hashlib,datetime
import numpy as np,pandas as pd,h5py,scipy,anndata
from anndata.io import read_elem
from scipy.sparse import load_npz,csr_matrix
from scipy.stats import rankdata,binomtest
O=Path(__file__).resolve().parent;R=O.parent;Q=R/'phase167_A_24section_endpoint_check_20260916_v1'
(O/'design_cache').mkdir(exist_ok=False)
inputs=pd.read_csv(R/'phase165_A_input_confirmation_20260915_v1/SECTION_INPUT_QC.tsv',sep='\t')
meta=pd.read_csv(R/'phase168_reference_source_audit_20260916_v1/SPATIAL_26_SECTION_PATIENT_MAP.tsv',sep='\t').set_index('section')
eligible=pd.read_csv(Q/'SECTION_ENDPOINT_SUMMARY.tsv',sep='\t').set_index('section')
genes='TYMP SRM SMS SLC16A3 ODC1 MTAP AMD1 MAT2A SPP1 COL1A1 COL1A2 DCN LUM CA9 EGLN3 VEGFA BNIP3 PDK1'.split()
manifest=[];results=[];checks=[];summaries=[];patients=[]
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def record(p):
 p=Path(p);manifest.append(dict(path=str(p),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns,sha256=sha(p)))
def check(g,k,v):
 checks.append(dict(section=g,check=k,passed=bool(v)))
 if not v:raise RuntimeError(str((g,k)))
def z(a):
 a=np.asarray(a,dtype=float)
 if not len(a):return a
 return (a-a.mean(axis=0))/np.where(a.std(axis=0)>0,a.std(axis=0),1)
for p in [R/'phase165_A_input_confirmation_20260915_v1/SECTION_INPUT_QC.tsv',R/'phase168_reference_source_audit_20260916_v1/SPATIAL_26_SECTION_PATIENT_MAP.tsv',Q/'SECTION_ENDPOINT_SUMMARY.tsv',Q/'FROZEN_CONFIG.json',O/'FROZEN_ANALYSIS_PLAN.md']:record(p)
for s in inputs.itertuples(index=False):
 gsm=s.section
 for p in [s.input_path,s.abundance_path,Q/(gsm+'_SPOT_TECHNICAL_QC.tsv.gz'),Q/(gsm+'_adjacency.npz')]:record(p)
 with h5py.File(s.input_path,'r') as f:
  obs=read_elem(f['obs']);var=read_elem(f['var']);c=read_elem(f['layers/counts']);xy=read_elem(f['obsm/spatial'])
 raw=c[:,[list(var.index).index(g) for g in genes]].toarray().astype(float)
 tot=np.asarray(c.sum(axis=1)).ravel();e=pd.DataFrame(np.log1p(raw/tot[:,None]*1e4),columns=genes)
 sq=pd.read_csv(Q/(gsm+'_SPOT_TECHNICAL_QC.tsv.gz'),sep='\t')
 check(gsm,'frozen_spot_order',list(obs.index)==sq.spot_id.tolist())
 ab=pd.read_csv(s.abundance_path,sep='\t',index_col=0)
 check(gsm,'abundance_spot_order',list(obs.index)==list(ab.index))
 q=ab.filter(regex='^q05__');fra=q.div(q.sum(axis=1),axis=0)
 mal=fra['q05__Malignant'].to_numpy()
 mc=['q05__C1QC-state Macrophage','q05__Other Macrophage','q05__SPP1-state Macrophage']
 mac=fra[mc].sum(axis=1).to_numpy()
 other=fra.drop(columns=['q05__Malignant']+mc).max(axis=1).to_numpy()
 sh=(mal>=.2)&(mal>=np.maximum(mac,other));mh=(mac>=.1)&(mac>=np.median(mac))
 check(gsm,'frozen_region_membership',np.array_equal(sh,sq.source_region) and np.array_equal(mh,sq.macrophage_region))
 A=load_npz(Q/(gsm+'_adjacency.npz'));rec=A.multiply(mh[None,:]).tocsr();rn=np.asarray(rec.sum(axis=1)).ravel()
 W=csr_matrix(rec.multiply(1/np.maximum(rn,1)[:,None]))
 check(gsm,'frozen_neighbor_counts',np.array_equal(rn,sq.eligible_macrophage_neighbors))
 masks={'source':sh,'receiver_link':sh&(rn>=3)}
 check(gsm,'frozen_receiver_membership',np.array_equal(masks['receiver_link'],sq.eligible_receiver_link))
 st=z(e[['COL1A1','COL1A2','DCN','LUM']]).mean(axis=1);hy=z(e[['CA9','EGLN3','VEGFA','BNIP3','PDK1']]).mean(axis=1)
 base=np.column_stack([np.log1p(tot),mal,mac]);ext=np.column_stack([base,st,hy])
 nbrbase=W@base;nbrext=W@ext;ne=W@e.to_numpy();coords=z(xy)
 coordpoly=np.column_stack([coords,coords[:,0]**2,coords[:,0]*coords[:,1],coords[:,1]**2])
 cache={};cachekey=0
 for role,mask0 in masks.items():
  Cbase=base if role=='source' else np.column_stack([base,nbrbase,e.TYMP])
  Cext=ext if role=='source' else np.column_stack([ext,nbrext,e.TYMP])
  for candidate in ['SRM','SMS','SLC16A3','MTAP']:
   modelnames=['extended_primary'] if candidate=='MTAP' else ['extended_primary','composition','mtap_adjusted','coordinate_adjusted','largest_component']
   if gsm=='GSM5924048' and candidate!='MTAP':
    modelnames += ['component_'+str(k) for k in sorted(sq.component.unique()) if int((mask0&(sq.component.to_numpy()==k)).sum())>=50]
   for model in modelnames:
    mask=mask0.copy();C=Cbase.copy() if model=='composition' else Cext.copy()
    if model=='mtap_adjusted':
     C=np.column_stack([C,e.MTAP]) if role=='source' else np.column_stack([C,e.MTAP,ne[:,genes.index('MTAP')]])
    if model=='coordinate_adjusted':C=np.column_stack([C,coordpoly])
    if model=='largest_component':mask &= sq.component.to_numpy()==sq.component.value_counts().idxmax()
    if model.startswith('component_'):mask &= sq.component.to_numpy()==int(model.split('_')[1])
    x=e.TYMP.to_numpy() if role=='source' else e[candidate].to_numpy()
    y=e[candidate].to_numpy() if role=='source' else ne[:,genes.index('SPP1')]
    n=int(mask.sum())
    xr=rankdata(x[mask]);yr=rankdata(y[mask]);cr=rankdata(C[mask],axis=0) if n else C[mask]
    X=np.column_stack([np.ones(n),z(xr),z(cr)])
    rank=int(np.linalg.matrix_rank(X)) if n else 0;cond=float(np.linalg.cond(X)) if n else np.inf
    reasons=[]
    if n<50:reasons.append('fewer_than_50_spots')
    if rank!=X.shape[1]:reasons.append('rank_deficient')
    if not np.isfinite(cond) or cond>=30:reasons.append('condition_number_ge30')
    if candidate!='MTAP' and not bool(eligible.loc[gsm,'source_all_designs_eligible' if role=='source' else 'receiver_all_designs_eligible']):reasons.append('frozen_endpoint_ineligible')
    rho=np.nan
    if not reasons:
     Cmat=X[:,np.r_[0,2:X.shape[1]]]
     rx=X[:,1]-Cmat@np.linalg.lstsq(Cmat,X[:,1],rcond=None)[0]
     ry=yr-Cmat@np.linalg.lstsq(Cmat,yr,rcond=None)[0]
     if np.std(rx)<1e-10 or np.std(ry)<1e-10:reasons.append('constant_residual')
     else:rho=float(np.dot(rx,ry)/np.sqrt(np.dot(rx,rx)*np.dot(ry,ry)))
    key='fit'+str(cachekey);cachekey+=1
    cache[key+'_X']=X;cache[key+'_y']=yr;cache[key+'_spots']=np.flatnonzero(mask)
    results.append(dict(section=gsm,patient_id=meta.loc[gsm,'patient_id'],cohort=meta.loc[gsm,'cohort'],preparation=s.preparation,role=role,candidate=candidate,model=model,spots=n,columns=X.shape[1],rank=rank,condition_number=cond,rho=rho,eligible=not reasons,reason=';'.join(reasons),cache_key=key))
 np.savez_compressed(O/'design_cache'/(gsm+'.npz'),**cache)
 print('section complete',gsm,flush=True)
d=pd.DataFrame(results);d.to_csv(O/'SECTION_EFFECTS.tsv',sep='\t',index=False)
pd.DataFrame(manifest).to_csv(O/'INPUT_MANIFEST.tsv',sep='\t',index=False)
pd.DataFrame(checks).to_csv(O/'INPUT_VALIDATION.tsv',sep='\t',index=False)
rng=np.random.default_rng(20260916)
def summarize(sub,role,candidate,model,stratum,policy,infer=True):
 grouped=[]
 for pid,a in sub[sub.eligible].groupby('patient_id',sort=True):
  rr=np.clip(a.rho.to_numpy(),-1+1e-12,1-1e-12)
  value=float(np.tanh(np.arctanh(rr).mean()))
  grouped.append(value)
  patients.append(dict(patient_id=pid,cohort=a.cohort.iloc[0],role=role,candidate=candidate,model=model,stratum=stratum,policy=policy,sections=len(a),section_ids=';'.join(a.section),rho=value))
 values=np.array(grouped);n=len(values);pos=int((values>0).sum());neg=int((values<0).sum());neff=pos+neg
 p=float(binomtest(pos,neff,p=.5,alternative='two-sided').pvalue) if infer and n>=10 and neff>=10 else np.nan
 lo=hi=np.nan
 if infer and n>=10:
  boot=np.median(values[rng.integers(0,n,size=(10000,n))],axis=1);lo,hi=np.quantile(boot,[.025,.975])
 fam='MTAP_background' if candidate=='MTAP' else role+'_candidates'
 summaries.append(dict(role=role,candidate=candidate,model=model,stratum=stratum,policy=policy,family=fam,n_patients=n,n_sections=int(sub.eligible.sum()),positive=pos,negative=neg,zero=n-neff,test_n=neff,median_rho=float(np.median(values)) if n else np.nan,q25=float(np.quantile(values,.25)) if n else np.nan,q75=float(np.quantile(values,.75)) if n else np.nan,median_ci_low=lo,median_ci_high=hi,P=p))
for (role,candidate,model),sub in d.groupby(['role','candidate','model'],sort=True):
 summarize(sub,role,candidate,model,'all','equal_fisher_mean',infer=not model.startswith('component_'))
 if model=='extended_primary':
  for prep in ['FFPE','frozen']:summarize(sub[sub.preparation==prep],role,candidate,model,'preparation_'+prep,'equal_fisher_mean')
  for cohort in ['ExhauCRF','IMM','BIONIKK']:
   summarize(sub[sub.cohort==cohort],role,candidate,model,'cohort_'+cohort,'equal_fisher_mean',False)
   summarize(sub[sub.cohort!=cohort],role,candidate,model,'exclude_cohort_'+cohort,'equal_fisher_mean',False)
  for policy,drop in [('c2_FFPE_only',['GSM5924050']),('c2_frozen_only',['GSM5924030']),('exclude_c2',['GSM5924030','GSM5924050'])]:
   summarize(sub[~sub.section.isin(drop)],role,candidate,model,'all',policy)
ss=pd.DataFrame(summaries);ss['FDR']=np.nan
for key,a in ss.groupby(['model','stratum','policy','family']):
 idx=a.index;p=a.P.fillna(1).to_numpy();N=2 if key[-1]=='MTAP_background' else 3
 assert len(p)<=N
 pp=np.r_[p,np.ones(N-len(p))];order=np.argsort(pp);q=np.minimum.accumulate((pp[order]*N/np.arange(1,N+1))[::-1])[::-1];out=np.empty(N);out[order]=np.minimum(q,1)
 ss.loc[idx,'FDR']=np.where(a.P.notna(),out[:len(p)],np.nan)
ss.to_csv(O/'PATIENT_SUMMARY.tsv',sep='\t',index=False)
pat=pd.DataFrame(patients);pat.to_csv(O/'PATIENT_EFFECTS.tsv',sep='\t',index=False)
loo=[]
for key,a in pat[(pat.model=='extended_primary')&(pat.stratum=='all')&(pat.policy=='equal_fisher_mean')].groupby(['role','candidate']):
 for x in a.itertuples(index=False):
  vals=a[a.patient_id!=x.patient_id].rho.to_numpy()
  loo.append(dict(role=key[0],candidate=key[1],omitted_patient=x.patient_id,n=len(vals),median_rho=float(np.median(vals)),positive=int((vals>0).sum()),negative=int((vals<0).sum())))
pd.DataFrame(loo).to_csv(O/'LEAVE_ONE_PATIENT.tsv',sep='\t',index=False)
p0=pat[(pat.model=='extended_primary')&(pat.stratum=='all')&(pat.policy=='equal_fisher_mean')]
p1=pat[(pat.model=='mtap_adjusted')&(pat.stratum=='all')&(pat.policy=='equal_fisher_mean')]
paired=p0.merge(p1,on=['patient_id','role','candidate'],suffixes=('_primary','_mtap'))
paired['rho_change']=paired.rho_mtap-paired.rho_primary;paired.to_csv(O/'MTAP_ADJUSTMENT_PAIRED.tsv',sep='\t',index=False)
status=dict(stage='Phase170',status='COMPUTATION_COMPLETE_PENDING_VALIDATION',section_effect_rows=len(d),patient_effect_rows=len(pat),summary_rows=len(ss),spot_level_tests=0,new_figures=0,python=sys.version,numpy=np.__version__,scipy=scipy.__version__,pandas=pd.__version__,anndata=anndata.__version__,h5py=h5py.__version__,completed_at=datetime.datetime.now().astimezone().isoformat())
(O/'COMPUTE_STATUS.json').write_text(json.dumps(status,indent=2))
print(ss[(ss.model=='extended_primary')&(ss.stratum=='all')&(ss.policy=='equal_fisher_mean')].to_string(index=False),flush=True)
print(json.dumps(status),flush=True)
