from pathlib import Path
import json,datetime,hashlib,ast,gc,traceback,subprocess
import numpy as np,pandas as pd,h5py
from scipy import sparse,stats
from anndata.io import read_elem
P=Path(__file__).resolve().parents[1];R=P.parent
OLD=R/'phase234_part4_replan_execution_20260921_v1'
MARK=R/'phase246_spatial_marker_structure_20260924_v1'
def now():return datetime.datetime.now().astimezone().isoformat()
def dump(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
def tab(d,f):pd.DataFrame(d).to_csv(P/'tables'/f,sep='\t',index=False)
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def gate(stage):
 q=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
 with open(P/'logs/memory_gate.log','a') as f:f.write(stage+'\n'+q.stdout+q.stderr)
 if q.returncode:raise RuntimeError('MEMORY_GATE')
def status(s,**kw):dump(P/'STATUS.json',dict(status=s,at=now(),**kw))
tree=ast.parse((OLD/'scripts/common.py').read_text())
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['fit','infer']],type_ignores=[]),'inherited','exec'),globals())
KEYS=['dataset','family','endpoint','model','effect_type','width_um']
def summarize(rows):
 d=pd.DataFrame(rows);tab(d,'SECTION_EFFECTS.tsv');out=[]
 for k,g in d.groupby(KEYS+['patient'],sort=True,dropna=False):
  v=g.effect.dropna().to_numpy()
  e=(float(np.tanh(np.arctanh(np.clip(v,-.999999999,.999999999)).mean())) if k[4]=='rho' else float(v.mean())) if len(v) else np.nan
  out.append(dict(zip(KEYS+['patient'],k),effect=e,n_sections=len(v),sections=';'.join(g.loc[g.effect.notna(),'sample'])))
 pat=pd.DataFrame(out);tab(pat,'PATIENT_EFFECTS.tsv');out=[]
 for i,(k,g) in enumerate(pat.groupby(KEYS,sort=True,dropna=False)):
  out.append(dict(zip(KEYS,k),**infer(g.effect,26210275+i)))
 su=pd.DataFrame(out)
 for _,g in su.groupby(['dataset','family']):
  for pc,qc in [('P','FDR'),('P_wilcoxon','FDR_wilcoxon')]:
   su.loc[g.index,qc]=np.where(g[pc].notna(),stats.false_discovery_control(g[pc].fillna(1).to_numpy()),np.nan)
  su.loc[g.index,'family_size_registered']=len(g);su.loc[g.index,'family_size_estimable']=g.P.notna().sum()
 su['evidence']=np.where(su.P.isna(),'DESCRIPTIVE',np.where(su.P<.05,np.where(su.FDR<.05,'FDR_SIGNIFICANT','NOMINAL_SIGNIFICANT'),'NOT_SUPPORTED'))
 tab(su,'COHORT_SUMMARY.tsv')
 loo=[];common=[]
 for k,g in pat[pat.width_um==400].groupby(KEYS,sort=True):
  v=g.loc[g.effect.notna()]
  if len(v)>=10:
   for rr in v.itertuples():
    z=v.loc[v.patient!=rr.patient,'effect'].to_numpy();nz=np.count_nonzero(z)
    loo.append(dict(zip(KEYS,k),omitted_patient=rr.patient,n_remaining=len(z),median=float(np.median(z)),positive=int((z>0).sum()),negative=int((z<0).sum()),P_sign=float(stats.binomtest(int((z>0).sum()),nz,.5).pvalue) if nz and len(z)>=10 else np.nan))
 for k,g in pat.groupby(['dataset','endpoint','model','effect_type']):
  q=g.pivot(index='patient',columns='width_um',values='effect').dropna()
  for w in [200,400,800]:
   vals=q[w].to_numpy() if w in q else np.array([])
   common.append(dict(zip(['dataset','endpoint','model','effect_type'],k),width_um=w,n_common=len(vals),median=float(np.median(vals)) if len(vals) else np.nan,positive=int((vals>0).sum()),negative=int((vals<0).sum())))
 tab(loo,'LEAVE_ONE_PATIENT_OUT.tsv');tab(common,'COMMON_PATIENT_WIDTH_SENSITIVITY.tsv')
 return pat,su
def ols_side(y,inside,mask,cov=None):
 good=mask&np.isfinite(y)
 if cov is not None:good &= np.isfinite(cov).all(1)
 z=inside[good].astype(float);v=y[good];ni=int(z.sum());no=len(z)-ni
 if min(ni,no)<20 or len(z)<50:return np.nan,'INSUFFICIENT_SIDE_SPOTS',np.nan
 C=np.empty((len(z),0)) if cov is None else cov[good]
 if C.shape[1]:
  sd=C.std(0)
  if (sd==0).any():return np.nan,'CONSTANT_COVARIATE',np.nan
  C=(C-C.mean(0))/sd
 D=np.column_stack([np.ones(len(z)),z,C]);cond=np.linalg.cond(D)
 if np.linalg.matrix_rank(D)<D.shape[1] or cond>=30:return np.nan,'RANK_OR_CONDITION',cond
 b=np.linalg.lstsq(D,v,rcond=None)[0][1]
 return float(b),'',float(cond)
def means(A,mask,v):
 count=np.asarray(A@mask.astype(float)).ravel()
 return np.asarray(A@np.where(mask,v,0.)).ravel()/np.maximum(count,1),count
def main():
 gate('start');freeze=json.loads((P/'FREEZE.json').read_text())
 assert sha(P/'PROTOCOL_BEFORE_TESTING.md')==freeze['protocol_sha256']
 regs=pd.read_csv(OLD/'tables/INPUT_SECTION_REGISTER.tsv',sep='\t');assert len(regs)==38 and regs.patient.nunique()==37
 regs.to_csv(P/'tables/INPUT_SECTION_REGISTER.tsv',sep='\t',index=False)
 elig=pd.read_csv(MARK/'tables/MARKER_ELIGIBILITY.tsv',sep='\t')
 effects=[];qc=[];readmeta=[];sources=[];side_stats=[];numeric=[]
 for i,r in enumerate(regs.itertuples()):
  gate(r.sample);status('RUNNING',completed=i,total=38,current=r.sample)
  paths=[OLD/'inputs'/f'{r.sample}_INPUTS.tsv.gz',OLD/'graphs'/f'{r.sample}_distance800.npz',MARK/'inputs'/f'{r.sample}_MARKERS.tsv.gz']
  sources += [dict(path=str(f),sha256=sha(f)) for f in paths]
  d=pd.read_csv(paths[0],sep='\t');m=pd.read_csv(paths[2],sep='\t');assert np.array_equal(d.spot_id,m.spot_id)
  assert d.spot_id.is_unique and r.physical_eligible
  inside=d.tumor_expression_region.to_numpy(bool);outside=~inside
  assert np.array_equal(inside,d.malignant_fraction.to_numpy()>=.20)
  boundary=d.expression_interface.to_numpy(bool);dist=d.expression_distance_um.to_numpy()
  assert np.all(dist[boundary]==0)
  D=sparse.load_npz(paths[1]).tocsr()
  assert D.shape==(len(d),len(d)) and not D.diagonal().any() and (D-D.T).nnz==0
  dr,dc=D.nonzero();assert np.all(d.component.to_numpy()[dr]==d.component.to_numpy()[dc])
  assert D.data.min()>0 and D.data.max()<=800.0001
  f=Path(r.input_path);st=f.stat();readmeta.append(dict(path=str(f),bytes=st.st_size,mtime_ns=st.st_mtime_ns))
  genes=['CD3D','CD3E','TRAC','CD8A','CD8B'];panels={'Tcell':genes[:3],'CD8lineage':genes}
  with h5py.File(f,'r') as h:
   obs=read_elem(h['obs']);var=read_elem(h['var']);assert np.array_equal(obs.index,d.spot_id) and var.index.is_unique
   X=read_elem(h['layers/counts'] if 'counts' in h.get('layers',{}) else h['X']).astype(float)
   assert sparse.issparse(X) and np.isfinite(X.data).all() and (X.data>=0).all()
   total=np.asarray(X.sum(1)).ravel();assert (total>0).all()
   used=[g for g in genes if g in var.index];raw=X[:,var.index.get_indexer(used)].toarray()
  log=np.log1p(raw/total[:,None]*10000);lookup={g:j for j,g in enumerate(used)}
  cd8detect=np.zeros(len(d),bool);maxerr=0
  for name,panel in panels.items():
   ids=[lookup[g] for g in panel if g in lookup];v=log[:,ids];sd=v.std(0);s=((v-v.mean(0))/np.where(sd>0,sd,1)).mean(1)
   detect=(raw[:,ids]>0).sum(1)>=2
   if name=='CD8lineage':
    detect &= raw[:,lookup['CD8B']]>0 if 'CD8B' in lookup else False
    detect &= (raw[:,[lookup[g] for g in genes[:3] if g in lookup]]>0).any(1)
    cd8detect=detect
   ok=(len(ids)==len(panel)) and detect.sum()>=20 and np.std(s)>0
   er=elig[(elig['sample']==r.sample)&(elig.marker==name)].iloc[0]
   assert ok==er.eligible and int(detect.sum())==er.detected_spots
   if ok:
    err=float(np.max(np.abs(s-m[name])));assert err<1e-10;maxerr=max(maxerr,err)
   else:assert m[name].isna().all()
  cd8ok=bool(elig[(elig['sample']==r.sample)&(elig.marker=='CD8lineage')].iloc[0].eligible)
  tq=d.Tcell_q05.to_numpy();yT=np.log2(tq) if (tq>0).all() else np.full(len(d),np.nan)
  ys={'T_abundance_log2':yT,'CD8_RNA':m.CD8lineage.to_numpy(),'T_RNA':m.Tcell.to_numpy(),'CD8_codetection':cd8detect.astype(float) if cd8ok else np.full(len(d),np.nan)}
  targets={'T_abundance_log2':'log2_difference','CD8_RNA':'zscore_difference','T_RNA':'zscore_difference','CD8_codetection':'proportion_difference'}
  depth=d.log1p_UMI.to_numpy();macro=d.total_macro_q05.to_numpy();mal=d.malignant_fraction.to_numpy()
  covband=np.column_stack([depth,np.log1p(macro)])
  def add(module,target,model,w,e,n,reason='',typ=None,condition=np.nan):
   effects.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,preparation=r.preparation,family=module+('_primary' if w==400 else '_sensitivity'),endpoint=module+'__'+target,model=model,effect_type=typ or targets[target],width_um=w,effect=e,eligible_spots=n,reason=reason,condition=condition))
  def adjacency(w):
   A=D.copy();A.data=(A.data<=w+1e-6).astype(float);A.eliminate_zeros();return A+sparse.eye(len(d),format='csr')
  A200=adjacency(200);den=np.asarray(A200@macro).ravel();exposure=np.asarray(A200@d.V1_SPP1_q05.to_numpy()).ravel()/den
  macro200=np.asarray(A200@macro).ravel()/np.asarray(A200.sum(1)).ravel()
  centers=boundary&(d.total_macro_fraction.to_numpy()>=.05)
  for w in [400,200,800]:
   mask=np.isfinite(dist)&(dist<=w+1e-6)
   for target,y in ys.items():
    module='band_primary_readouts' if target in ['T_abundance_log2','CD8_RNA'] else 'band_auxiliary_readouts'
    vi=mask&inside&np.isfinite(y);vo=mask&outside&np.isfinite(y)
    side_stats.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,width_um=w,target=target,n_inside=int(vi.sum()),n_outside=int(vo.sum()),mean_inside=float(y[vi].mean()) if vi.any() else np.nan,mean_outside=float(y[vo].mean()) if vo.any() else np.nan,raw_T_mean_inside=float(tq[mask&inside].mean()) if (mask&inside).any() else np.nan,raw_T_mean_outside=float(tq[mask&outside].mean()) if (mask&outside).any() else np.nan))
    for model in ['unadjusted','depth_macro_adjusted']:
     e,why,cond=ols_side(y,inside,mask,None if model=='unadjusted' else covband)
     add(module,target,model,w,e,int(vi.sum()+vo.sum()),why,condition=cond)
     if model=='unadjusted' and np.isfinite(e):assert abs(e-(y[vi].mean()-y[vo].mean()))<1e-9
   A=adjacency(w)
   di,ni=means(A,inside,depth);do,no=means(A,outside,depth)
   mi,_=means(A,inside,mal);mo,_=means(A,outside,mal)
   nn=ni+no;meanD=(di*ni+do*no)/nn;meanM=(mi*ni+mo*no)/nn
   C=np.column_stack([meanD,di-do,np.log1p(macro200),meanM,mi-mo,np.log(np.maximum(ni,1)/np.maximum(no,1))])
   cmask=centers&(ni>=5)&(no>=5)
   cd=pd.DataFrame(dict(spot_id=d.spot_id,center_eligible=cmask,SPP1_share200=exposure,n_inside=ni,n_outside=no))
   for j,name in enumerate(['depth_mean','depth_inside_minus_outside','log1p_macro200','malignant_mean','malignant_inside_minus_outside','log_side_count_ratio']):cd[name]=C[:,j]
   for target in ['T_abundance_log2','CD8_RNA']:
    yi,_=means(A,inside,ys[target]);yo,_=means(A,outside,ys[target]);diff=yi-yo
    cd[target+'_inside_minus_outside']=diff
    for model in ['unadjusted','depth_macro_composition_adjusted']:
     e,n,why=fit(exposure,diff,None if model=='unadjusted' else C,mask=cmask)
     add('SPP1_coupling',target,model,w,e,n,why,'rho')
   cd.loc[centers].to_csv(P/'inputs'/f'{r.sample}_centers_{w}.tsv.gz',sep='\t',index=False)
   # Save a deterministic local geometry calculation for independent verification.
   candidates=np.flatnonzero(cmask)
   if len(candidates):
    j=int(candidates[0]);nb=A[j].indices
    numeric.append(dict(sample=r.sample,width_um=w,center_index=j,center_id=str(d.spot_id.iloc[j]),n_inside=int(inside[nb].sum()),n_outside=int(outside[nb].sum()),T_delta=float(yT[nb][inside[nb]].mean()-yT[nb][outside[nb]].mean()),saved_T_delta=float(cd['T_abundance_log2_inside_minus_outside'].iloc[j])))
   qc.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,width_um=w,spots=len(d),tumor_spots=int(inside.sum()),outside_spots=int(outside.sum()),interface_spots=int(boundary.sum()),reachable_spots=int(np.isfinite(dist).sum()),band_inside=int((mask&inside).sum()),band_outside=int((mask&outside).sum()),macro_interface_centers=int(centers.sum()),paired_eligible_centers=int(cmask.sum()),CD8_eligible=cd8ok,CD8_detected_spots=int(cd8detect.sum()),marker_reconstruction_error=maxerr,components=d.component.nunique()))
  pd.DataFrame(dict(spot_id=d.spot_id,CD8_codetection=cd8detect)).to_csv(P/'inputs'/f'{r.sample}_CD8_detection.tsv.gz',sep='\t',index=False)
  print('SECTION_COMPLETE',i+1,r.sample,flush=True)
  del X,raw,log,D,A,A200,cd;gc.collect()
 tab(qc,'SECTION_GEOMETRY_ELIGIBILITY.tsv');tab(side_stats,'BAND_SIDE_DESCRIPTIVES.tsv');tab(numeric,'LOCAL_GEOMETRY_NUMERIC_CHECK.tsv')
 pat,su=summarize(effects)
 for row in readmeta:
  st=Path(row['path']).stat();assert st.st_size==row['bytes'] and st.st_mtime_ns==row['mtime_ns']
 for row in sources:assert sha(row['path'])==row['sha256']
 dump(P/'qc/READONLY_INPUTS.json',dict(raw_h5ad_metadata=readmeta,derived_sha256=sources))
 assert np.nanmax(np.abs(pd.DataFrame(numeric).T_delta-pd.DataFrame(numeric).saved_T_delta))<1e-10
 dump(P/'qc/COMPUTATION_COMPLETED.json',dict(status='COMPUTED_AWAITING_INDEPENDENT_VALIDATION',at=now(),sections=38,patients=37,section_effects=len(effects),cohort_endpoints=len(su),figures=0,protocol_sha256=sha(P/'PROTOCOL_BEFORE_TESTING.md')))
 status('COMPUTED_AWAITING_INDEPENDENT_VALIDATION',completed=38,total=38)
 print('ALL_COMPUTATION_COMPLETED',len(su),flush=True)
if __name__=='__main__':
 try:main()
 except Exception:
  status('FAILED',traceback=traceback.format_exc());raise
