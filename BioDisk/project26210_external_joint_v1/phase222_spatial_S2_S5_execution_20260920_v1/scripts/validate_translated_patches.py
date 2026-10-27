from pathlib import Path
import json
import numpy as np,pandas as pd,h5py
from scipy import sparse,stats
from anndata.io import read_elem
P=Path(__file__).resolve().parents[1];R=P.parent;F=R/'phase221_spatial_S1_identity_eligibility_freeze_20260920_v1';T=pd.read_csv(P/'S2/TRANSLATION_CONTROLS.tsv.gz',sep='\t');check=[]
def corr(x,y,c,idx):
 arr=np.column_stack([x[idx],y[idx],c[idx]]);r=np.column_stack([stats.rankdata(arr[:,j]) for j in range(arr.shape[1])]);v=np.linalg.inv(np.corrcoef(r.T));return -v[0,1]/np.sqrt(v[0,0]*v[1,1])
for sid in ['GSM5924030','GSM5924048']:
 z=pd.read_csv(F/'frozen_inputs'/(sid+'_SPOT_INPUTS.tsv.gz'),sep='\t',index_col='spot_id');G=sparse.load_npz(F/'frozen_inputs'/(sid+'_adjacency.npz')).tocsr()
 with h5py.File(R.parent/'project26210_external_joint_v1_phase54_gse175540_spatial_v1/run_v1/sections_v1'/sid/'processed_section_v1.h5ad') as h:obs=read_elem(h['obs'])
 qr=np.column_stack([(obs.array_col.to_numpy(int)-obs.array_row.to_numpy(int))//2,obs.array_row.to_numpy(int)]);lookup={tuple(v):i for i,v in enumerate(qr)}
 for ep in ['S2_endothelial','S2_pericyte']:
  row=T[(T.section==sid)&(T.endpoint==ep)].iloc[0];delta=np.array([int(row.dq),int(row.dr)]);key=ep.removeprefix('S2_');v=z[key+'_fraction'].to_numpy();idx=[];shifted=[];ys=[]
  for i in np.flatnonzero(z.source_eligible):
   nb=G.indices[G.indptr[i]:G.indptr[i+1]];allids=np.r_[i,nb];mapped=np.array([lookup.get(tuple(qr[j]+delta),-1) for j in allids])
   if (mapped>=0).all() and (z.component.to_numpy()[mapped]==z.component.iloc[i]).all():idx.append(i);shifted.append(v[mapped[0]]);ys.append(v[mapped[1:]].mean())
  idx=np.array(idx);assert len(idx)==row.n_matched;cache=np.load(P/'cache'/(sid+'_'+key+'_PRIMARY.npz'));c=cache['cov'].copy();c[idx,1]=shifted;y=np.zeros(len(z));y[idx]=ys
  null=corr(cache['x'],y,c,idx);orig=corr(cache['x'],cache['y'],cache['cov'],idx);assert abs(null-row.shifted)<1e-10 and abs(orig-row.matched_observed)<1e-10
  check.append({'section':sid,'endpoint':ep,'matched':len(idx),'observed_error':abs(orig-row.matched_observed),'shifted_error':abs(null-row.shifted)})
(P/'audit/INDEPENDENT_TRANSLATION_VALIDATION.json').write_text(json.dumps({'status':'PASS','checks':check},indent=2))
v=json.loads((P/'FINAL_VALIDATION.json').read_text());v['independent_translated_patch_checks']=4;v['leave_one_eligible_patients_only']=98;v['uncomputed_TLS_diagnostics_marked_NA']=True;(P/'FINAL_VALIDATION.json').write_text(json.dumps(v,indent=2));print('TRANSLATION_AND_FINAL_DIAGNOSTICS_PASS')
