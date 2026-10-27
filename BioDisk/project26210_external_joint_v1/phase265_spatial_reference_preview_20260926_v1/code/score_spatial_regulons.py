from pathlib import Path
import os,json,pickle,hashlib,sys,gc
os.environ['OMP_NUM_THREADS']='1'
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
import numpy as np,pandas as pd,h5py
from scipy.sparse import csr_matrix
from ctxcore.genesig import Regulon
from pyscenic.aucell import aucell,create_rankings,aucell4r
P=Path(__file__).resolve().parents[1];R=P.parent
for d in ['code','tables','logs','qc','delivery','plot_data']:(P/d).mkdir(exist_ok=True)
tumor=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/07_pyscenic_v1/intermediate/regulons_grnboost2_v1.dat')
myeloid=R/'part2_rebuild_v1/08_regulons/network_v2_20260915/regulons.dat'
regs=[];sources=[]
for f,targets in [(tumor,['HNF4A','HOXD3','ILF2']),(myeloid,['JUNB','CEBPD'])]:
 with f.open('rb') as h:allregs=pickle.load(h)
 for tf in targets:
  found=[v for v in allregs if v.name==tf+'(+)'];assert len(found)==1,(tf,len(found))
  regs.append(found[0]);sources.append(dict(tf=tf,source=str(f),n_targets=len(found[0].gene2weight),sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
(P/'qc/TF_NETWORK_SOURCES.json').write_text(json.dumps(sources,indent=2))
reg=pd.read_csv(R/'phase234_part4_replan_execution_20260921_v1/tables/INPUT_SECTION_REGISTER.tsv',sep='\t')
coverage=[];valid=[];inputs=[]
for sample in ['GLMF1','GLMF2']:
 row=reg[reg['sample']==sample].iloc[0];path=Path(row.input_path)
 with h5py.File(path,'r') as f:
  genes=f['var/_index'].asstr()[:];ids=f['obs/_index'].asstr()[:];g=f['X']
  x=csr_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=tuple(g.attrs['shape']))
 assert x.shape==(len(ids),len(genes)) and pd.Index(genes).is_unique and pd.Index(ids).is_unique
 assert np.isfinite(x.data).all() and np.min(x.data)>=0
 order=np.argsort(genes);genes=genes[order];x=x[:,order];dense=x.toarray().astype(np.float32)
 frame=pd.DataFrame(dense,index=ids,columns=genes)
 ranks=create_rankings(frame,seed=26210)
 for v in regs:
  targets=set(v.gene2weight);present=targets&set(genes)
  coverage.append(dict(sample=sample,regulon=v.name,all_targets=len(targets),present_targets=len(present),fraction=len(present)/len(targets)))
  coverage[-1]['eligible'] = len(present)>=20 and len(present)/len(targets)>=.8
 eligible=[v for v in regs if next(z for z in coverage if z['sample']==sample and z['regulon']==v.name)['eligible']]
 score=aucell4r(ranks,eligible,auc_threshold=.05,noweights=False,normalize=False,num_workers=1)
 assert np.isfinite(score.values).all() and np.min(score.values)>=0
 score.index.name='spot_id';score.to_csv(P/f'plot_data/{sample}_spatial_RNA_regulon_AUC.tsv',sep='\t')
 # Independent direct weighted recovery-area calculation from the exact frozen ranks.
 rankcut=int(round(.05*len(genes)))-1
 for v in eligible:
  pairs=[(k,w) for k,w in v.gene2weight.items() if k in ranks.columns]
  targetr=ranks[[k for k,w in pairs]].to_numpy()
  weights=np.array([w for k,w in pairs])
  direct=(np.maximum(rankcut-targetr.astype(float),0)*weights[None,:]).sum(axis=1)/((rankcut+1)*weights.sum())
  err=float(np.max(np.abs(direct-score[v.name].to_numpy())))
  valid.append(dict(sample=sample,regulon=v.name,max_abs_error=err,rankcut=rankcut))
 inputs.append(dict(path=str(path),bytes=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns,spots=len(ids),genes=len(genes)))
 print('SCORED',sample,len(ids),score.columns.tolist(),flush=True)
 del frame,dense,ranks,x;gc.collect()
pd.DataFrame(coverage).to_csv(P/'tables/TF_TARGET_COVERAGE.tsv',sep='\t',index=False)
(P/'qc/AUC_DIRECT_CHECK.json').write_text(json.dumps(valid,indent=2))
(P/'qc/SPATIAL_INPUTS.json').write_text(json.dumps(inputs,indent=2))
print('AUC_ERRORS',valid,flush=True)

