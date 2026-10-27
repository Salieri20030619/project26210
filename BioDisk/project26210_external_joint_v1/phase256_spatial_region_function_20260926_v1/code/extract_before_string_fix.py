from pathlib import Path
import os,sys,json,hashlib,datetime,subprocess,importlib.metadata as im
import numpy as np,pandas as pd,h5py
from scipy import sparse
from anndata.io import read_elem
P=Path(__file__).resolve().parents[1];ROOT=P.parent
A=ROOT/'phase234_part4_replan_execution_20260921_v1'
GATE='/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'
MODELS=['absolute','adjusted','macro_share']
def stamp():return datetime.datetime.now().astimezone().isoformat()
def gate(label):
 r=subprocess.run([GATE],capture_output=True,text=True)
 with (P/'logs/memory_gate.log').open('a') as f:f.write(label+'\t'+r.stdout+r.stderr)
 if r.returncode:raise RuntimeError(r.stdout+r.stderr)
def fingerprint(path):
 p=Path(path);st=p.stat();out=dict(path=str(p.resolve()),bytes=st.st_size,mtime_ns=st.st_mtime_ns,inode=st.st_ino)
 if st.st_size<64000000:out['sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
 else:
  with p.open('rb') as f:
   h=hashlib.sha256(f.read(1048576));f.seek(-1048576,2);h.update(f.read(1048576))
  out['head_tail_1MiB_sha256']=h.hexdigest()
 return out
def writej(name,obj):(P/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
def main():
 gate('extract_start')
 reg=pd.read_csv(P/'inputs/INPUT_SECTION_REGISTER.tsv',sep='\t').fillna('')
 assert len(reg)==38 and reg[['dataset','patient']].drop_duplicates().shape[0]==37
 assert reg['sample'].is_unique
 vers={s:im.version(s) for s in ['numpy','pandas','scipy','h5py','anndata']}
 vers.update(python=sys.executable,time=stamp(),thread_environment={s:os.environ.get(s) for s in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']})
 writej('qc/PYTHON_ENVIRONMENT.json',vers)
 manifest=[];regions=[];checks=[];sample_index=[]
 for row in reg.itertuples():
  sid=row.sample;gate('extract_'+sid)
  sp=A/'inputs'/f'{sid}_INPUTS.tsv.gz';paths=[Path(row.input_path),sp,Path(row.abundance_path)]
  before=[fingerprint(f) for f in paths];manifest.extend(before)
  d=pd.read_csv(sp,sep='\t')
  assert d.spot_id.is_unique and len(d)==int(row.spots)
  with h5py.File(row.input_path,'r') as h:
   obs=read_elem(h['obs']);var=read_elem(h['var'])
   assert np.array_equal(obs.index.astype(str),d.spot_id.astype(str))
   source='layers/counts' if 'counts' in h.get('layers',{}) else 'X'
   X=read_elem(h[source]).tocsr().astype(np.float64)
  genes=np.asarray(var.index.astype(str));assert len(set(genes))==len(genes)
  assert X.shape==(len(d),len(genes)) and np.isfinite(X.data).all() and (X.data>=0).all()
  assert np.max(np.abs(X.data-np.round(X.data)))<1e-6
  totals=np.asarray(X.sum(1)).ravel();assert (totals>0).all()
  np.testing.assert_allclose(np.log1p(totals),d.log1p_UMI,rtol=1e-10,atol=1e-10)
  L=X.multiply((10000/totals)[:,None]).tocsr();L.data=np.log1p(L.data)
  effects=np.full((3,len(genes)),np.nan);pb=np.full_like(effects,np.nan)
  detected=np.zeros((3,len(genes)),int);region_matrix=np.full((3,len(d)),-1,np.int8)
  for mi,model in enumerate(MODELS):
   gate('model_'+sid+'_'+model)
   eligible=(d.total_macro_fraction.to_numpy()>=.05) if model=='macro_share' else np.ones(len(d),bool)
   v=d['V3_SPP1_macro_share' if model=='macro_share' else 'V1_SPP1_q05'].to_numpy(float)
   eligible &= np.isfinite(v)
   values=v[eligible];threshold=float(np.quantile(values,.8)) if len(values) else np.nan
   high=eligible&(v>=threshold);low=eligible&~high
   nh=int(high.sum());nl=int(low.sum());fraction=nh/(nh+nl) if nh+nl else np.nan
   region_matrix[mi,low]=0;region_matrix[mi,high]=1
   rr=dict(sample=sid,dataset=row.dataset,patient=row.patient,model=model,eligible_spots=int(eligible.sum()),high_spots=nh,background_spots=nl,threshold=threshold,high_fraction=fraction,eligible=False,reason='',condition_number=np.nan,group_information_fraction=np.nan)
   for c in ['log1p_UMI','total_macro_fraction','malignant_fraction','endothelial_fraction','pericyte_fraction']:
    rr[c+'_high']=float(d.loc[high,c].mean()) if nh else np.nan
    rr[c+'_background']=float(d.loc[low,c].mean()) if nl else np.nan
   if nh<50 or nl<50 or not(.1<=fraction<=.35) or np.ptp(values)<=1e-12:
    rr['reason']='REGION_SIZE_OR_DEGENERATE_THRESHOLD';regions.append(rr);continue
   Z=L[eligible];m=np.asarray(Z.mean(0)).ravel()
   sd=np.sqrt(np.maximum(0,np.asarray(Z.power(2).mean(0)).ravel()-m*m))
   variable=sd>1e-10
   detect=np.asarray((X[eligible]>0).sum(0)).ravel();detected[mi]=detect
   if model=='adjusted':
    cols=['log1p_UMI','total_macro_fraction','malignant_fraction','endothelial_fraction','pericyte_fraction']
    C=d.loc[eligible,cols].to_numpy(float);assert np.isfinite(C).all()
    cs=C.std(0);C=(C[:,cs>1e-10]-C[:,cs>1e-10].mean(0))/cs[cs>1e-10]
    C=np.column_stack([np.ones(len(C)),C]);cond=float(np.linalg.cond(C));rr['condition_number']=cond
    high_e=high[eligible].astype(float);res=high_e-C@np.linalg.lstsq(C,high_e,rcond=None)[0]
    info=float(np.dot(res,res)/np.sum((high_e-high_e.mean())**2));rr['group_information_fraction']=info
    if np.linalg.matrix_rank(C)!=C.shape[1] or cond>=30 or info<.05:
     rr['reason']='COVARIATE_COLLINEARITY_OR_NO_RESIDUAL_GROUP_INFORMATION';regions.append(rr);continue
    w=res/np.dot(res,res)
   else:
    w=np.where(high[eligible],1/nh,-1/nl)
    sumh=np.asarray(X[high].sum(0)).ravel();suml=np.asarray(X[low].sum(0)).ravel()
    pb[mi]=np.log2(sumh/sumh.sum()*1e6+.5)-np.log2(suml/suml.sum()*1e6+.5)
   diff=np.asarray(Z.T@w).ravel()
   effects[mi]=np.divide(diff,sd,out=np.zeros(len(sd)),where=variable)
   rr['eligible']=True;rr['variable_genes']=int(variable.sum())
   # Independent dense recomputation on fixed lexical expressed genes.
   candidates=np.flatnonzero((detect>=10)&variable)
   fixed=candidates[np.linspace(0,len(candidates)-1,min(12,len(candidates))).astype(int)]
   dense=Z[:,fixed].toarray();dense=(dense-dense.mean(0))/dense.std(0)
   if model=='adjusted':
    ref=np.linalg.lstsq(np.column_stack([C,high_e]),dense,rcond=None)[0][-1]
   else:ref=dense[high[eligible]].mean(0)-dense[low[eligible]].mean(0)
   err=float(np.max(np.abs(ref-effects[mi,fixed])));assert err<1e-9,(sid,model,err)
   checks.append(dict(sample=sid,model=model,genes_checked=len(fixed),max_error=err))
   regions.append(rr)
  np.savez_compressed(P/'effects'/f'{sid}.npz',genes=genes,effects=effects,pseudobulk_log2CPM_difference=pb,detected_spots=detected,models=np.asarray(MODELS))
  pd.DataFrame({'spot_id':d.spot_id,**{m:region_matrix[i] for i,m in enumerate(MODELS)}}).to_csv(P/'inputs'/f'{sid}_REGIONS.tsv.gz',sep='\t',index=False)
  after=[fingerprint(f) for f in paths];assert before==after,(sid,'source changed')
  sample_index.append(dict(sample=sid,dataset=row.dataset,patient=row.patient,preparation=row.preparation,cohort=row.cohort,spots=len(d),genes=len(genes),counts_source=source,counts_sum=float(totals.sum())))
  print('EXTRACTED',sid,len(d),len(genes),flush=True)
  del X,L,Z
 pd.DataFrame(regions).to_csv(P/'tables/REGION_ELIGIBILITY.tsv',sep='\t',index=False)
 pd.DataFrame(checks).to_csv(P/'qc/INDEPENDENT_GENE_EFFECT_CHECKS.tsv',sep='\t',index=False)
 pd.DataFrame(sample_index).to_csv(P/'tables/SECTION_IDENTITY.tsv',sep='\t',index=False)
 writej('qc/SOURCE_MANIFEST_BEFORE.json',manifest)
 writej('EXTRACTION_COMPLETE.json',dict(status='PASS',time=stamp(),sections=len(sample_index),models=MODELS,checks=len(checks)))
 print('EXTRACTION_COMPLETE',flush=True)
if __name__=='__main__':main()
