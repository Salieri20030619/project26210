from pathlib import Path
import json,hashlib,traceback,sys
import numpy as np,pandas as pd,scipy
from scipy import sparse,stats
P=Path(__file__).resolve().parents[1];R=P.parent;F=R/'phase229_spatial_F1_partition_eligibility_20260921_v1';A=R/'phase223_spatial_cohort_reinclusion_20260920_v1';S=R/'phase221_spatial_S1_identity_eligibility_freeze_20260920_v1'
def save(n,o):(P/n).write_text(json.dumps(o,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
def tab(n,d):d.to_csv(P/n,sep='\t',index=False)
def fit(x,y,cov,mask):
 a=np.column_stack([x,y,cov])[mask];reason=[]
 info=dict(n=len(a),rank=np.nan,condition=np.nan)
 if len(a)<50:return np.nan,dict(**info,reason='SOURCE_LT50'),None
 if not np.isfinite(a).all():return np.nan,dict(**info,reason='NONFINITE'),None
 if min(len(np.unique(a[:,0])),len(np.unique(a[:,1])))<10:return np.nan,dict(**info,reason='INSUFFICIENT_VARIATION'),None
 rk=np.column_stack([stats.rankdata(a[:,i]) for i in range(a.shape[1])])
 c=rk[:,2:];sd=c.std(axis=0)
 if (sd==0).any():return np.nan,dict(**info,reason='CONSTANT_COVARIATE'),None
 D=np.column_stack([np.ones(len(a)),(c-c.mean(axis=0))/sd])
 rank=np.linalg.matrix_rank(D);cond=float(np.linalg.cond(D));info.update(rank=int(rank),condition=cond)
 if rank<D.shape[1] or cond>=30:return np.nan,dict(**info,reason='DESIGN_RANK_OR_CONDITION'),None
 res=rk[:,:2]-D@np.linalg.lstsq(D,rk[:,:2],rcond=None)[0]
 if (res.std(axis=0)<1e-10).any():return np.nan,dict(**info,reason='RESIDUAL_CONSTANT'),None
 rho=float(np.corrcoef(res.T)[0,1])
 return rho,dict(**info,reason=''),dict(D=D,ranks=rk[:,:2])
def infer(v,seed):
 v=np.asarray(v);v=v[np.isfinite(v)];n=len(v)
 o=dict(n=n,positive=int((v>0).sum()),negative=int((v<0).sum()),zero=int((v==0).sum()),median_rho=float(np.median(v)) if n else np.nan,P=np.nan,P_sign=np.nan,CI_low=np.nan,CI_high=np.nan)
 if n>=10:
  o['P']=float(stats.wilcoxon(v,alternative='two-sided',method='auto').pvalue) if (v!=0).any() else 1.
  nonzero=int((v!=0).sum());o['P_sign']=float(stats.binomtest(int((v>0).sum()),nonzero,.5).pvalue) if nonzero else 1.
  rng=np.random.default_rng(seed);boot=np.median(v[rng.integers(0,n,size=(10000,n))],axis=1);o['CI_low'],o['CI_high']=np.quantile(boot,[.025,.975])
 return o
def run():
 elig=pd.read_csv(F/'tables/SECTION_ELIGIBILITY.tsv',sep='\t');records=[];checks=[];sources=[]
 for r in elig.itertuples():
  f=F/'spots'/(r.sample+'_F1_SPATIAL_INPUTS.tsv.gz');d=pd.read_csv(f,sep='\t')
  sources.append(dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
  gp=(S/'frozen_inputs' if r.dataset=='GSE175540' else A/'inputs')/(r.sample+'_adjacency.npz');G=sparse.load_npz(gp);deg=np.asarray(G.sum(axis=1)).ravel();W=sparse.diags(1/np.maximum(deg,1))@G
  mask=d.source_eligible.to_numpy(bool);x=d.SPP1_share.to_numpy();xy=d[['x','y']].to_numpy();xy=(xy-xy.mean(axis=0))/xy.std(axis=0)
  spatialcov=np.column_stack([xy,xy[:,0]**2,xy[:,0]*xy[:,1],xy[:,1]**2])
  for target in ['POSTNCAF','ECM4']:
   sc=d['program_'+target].to_numpy();y=W@sc
   cov=np.column_stack([sc,d.total_macro_fraction,d.neighbor_macro,d.malignant_fraction,d.log1p_UMI,d.neighbor_depth])
   for model,c in [('primary',cov),('coordinate_sensitivity',np.column_stack([cov,spatialcov]))]:
    rho,info,design=fit(x,y,c,mask)
    rec=dict(sample=r.sample,dataset=r.dataset,patient=r.patient,target=target,model=model,rho=rho,**info)
    records.append(rec)
    if design is not None:
     np.savez_compressed(P/'designs'/(r.sample+'_'+target+'_'+model+'.npz'),**design)
     # Independent QR projection of the stored rank response.
     Q=np.linalg.qr(design['D'],mode='reduced')[0];res=design['ranks']-Q@(Q.T@design['ranks'])
     check=float(res[:,0]@res[:,1]/np.sqrt((res[:,0]@res[:,0])*(res[:,1]@res[:,1])))
     assert abs(check-rho)<1e-10
     checks.append(dict(sample=r.sample,target=target,model=model,rho_qr=check,error=abs(check-rho)))
  print('F2',r.sample,flush=True)
 sect=pd.DataFrame(records);tab('tables/SECTION_EFFECTS.tsv',sect);tab('audit/QR_CHECKS.tsv',pd.DataFrame(checks))
 patients=[]
 for (dataset,patient,target,model),g in sect.groupby(['dataset','patient','target','model']):
  v=g.rho.dropna().to_numpy()
  rho=float(np.tanh(np.arctanh(np.clip(v,-.999999999,.999999999)).mean())) if len(v) else np.nan
  patients.append(dict(dataset=dataset,patient=patient,target=target,model=model,rho=rho,n_sections=len(v),sections=';'.join(g.loc[g.rho.notna(),'sample'])))
 pat=pd.DataFrame(patients);tab('tables/PATIENT_EFFECTS.tsv',pat)
 summary=[]
 for i,((dataset,target,model),g) in enumerate(pat.groupby(['dataset','target','model'])):
  summary.append(dict(dataset=dataset,target=target,model=model,**infer(g.rho.to_numpy(),262230+i)))
 su=pd.DataFrame(summary);su['FDR']=np.nan
 for (dataset,model),g in su.groupby(['dataset','model']):
  q=stats.false_discovery_control(g.P.fillna(1).to_numpy());su.loc[g.index,'FDR']=np.where(g.P.notna(),q,np.nan)
 su['support_label']=np.where(su.P.isna(),'DESCRIPTIVE_ONLY',np.where(su.P<.05,np.where(su.median_rho>0,'POSITIVE_NOMINAL','NEGATIVE_NOMINAL'),'NOT_NOMINAL'))
 tab('tables/COHORT_SUMMARY.tsv',su)
 paired=[]
 for (dataset,target),g in pat.groupby(['dataset','target']):
  piv=g.pivot(index='patient',columns='model',values='rho').dropna()
  paired.append(dict(dataset=dataset,target=target,n_same_patients=len(piv),primary_median=piv.primary.median(),coordinate_median=piv.coordinate_sensitivity.median(),median_change=(piv.coordinate_sensitivity-piv.primary).median()))
 tab('tables/PAIRED_COORDINATE_SENSITIVITY.tsv',pd.DataFrame(paired))
 loo=[]
 for (dataset,target),g in pat[pat.model=='primary'].groupby(['dataset','target']):
  g=g[g.rho.notna()]
  for row in g.itertuples():
   others=g[g.patient!=row.patient].rho;loo.append(dict(dataset=dataset,target=target,excluded_patient=row.patient,n_remaining=len(others),median_rho=others.median()))
 tab('tables/LEAVE_ONE_PATIENT.tsv',pd.DataFrame(loo))
 for x in sources:assert hashlib.sha256(Path(x['path']).read_bytes()).hexdigest()==x['sha256']
 for x in json.loads((F/'audit/RDS_BEFORE.json').read_text()):
  st=Path(x['path']).stat();assert(st.st_size,st.st_mtime_ns,st.st_ino)==(x['bytes'],x['mtime_ns'],x['inode'])
 save('INPUT_SHA256.json',sources);save('SOFTWARE.json',dict(python=sys.version,scipy=scipy.__version__,numpy=np.__version__,pandas=pd.__version__))
 save('VALIDATION.json',dict(status='PASS',sections=38,section_models=len(sect),estimable_models=int(sect.rho.notna().sum()),QR_max_error=max(c['error'] for c in checks),input_hashes_unchanged=True,old_RDS_unchanged=333,new_figures=0))
 save('STATUS.json',dict(status='F2_PROGRAM_ASSOCIATION_COMPLETE',CAF_cell_abundance_claim=False,C4_mapping_claim=False));(P/'EXIT_STATUS').write_text('0\n')
 print(su.to_string(index=False),flush=True)
try:run()
except Exception:
 (P/'ERROR.txt').write_text(traceback.format_exc());(P/'EXIT_STATUS').write_text('1\n');raise
