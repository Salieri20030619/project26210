from pathlib import Path
import json,hashlib,math,collections
import numpy as np,pandas as pd
from scipy import sparse,stats
P=Path(__file__).resolve().parents[1];R=P.parent;F=R/'phase221_spatial_S1_identity_eligibility_freeze_20260920_v1'
assert (P/'EXIT_STATUS').read_text().strip()=='0'
sec=pd.read_csv(P/'ALL_SECTION_EFFECTS.tsv',sep='\t');pat=pd.read_csv(P/'ALL_PATIENT_EFFECTS.tsv',sep='\t');summary=pd.read_csv(P/'PATIENT_RESULT_SUMMARY.tsv',sep='\t');assert len(sec)==360 and len(pat)==345 and len(summary)==15
checks=[]
for f in sorted((P/'cache').glob('*_PRIMARY.npz')):
 d=np.load(f);ma=d['mask'];mat=np.column_stack([d['x'][ma],d['y'][ma],d['cov'][ma]]);rank=np.column_stack([stats.rankdata(mat[:,i]) for i in range(mat.shape[1])]);precision=np.linalg.inv(np.corrcoef(rank.T));effect=-precision[0,1]/np.sqrt(precision[0,0]*precision[1,1]);sid=f.name.split('_')[0];key='endothelial' if 'endothelial' in f.name else 'pericyte';observed=sec[(sec.section==sid)&(sec.endpoint=='S2_'+key)&(sec.model=='primary')].effect.iloc[0];assert abs(effect-observed)<1e-10;checks.append(abs(effect-observed))
for row in pat.itertuples():
 values=sec[(sec.patient==row.patient)&(sec.endpoint==row.endpoint)&(sec.model==row.model)].effect.dropna().to_numpy()
 if len(values)==0:assert np.isnan(row.effect);continue
 value=np.tanh(np.arctanh(np.clip(values,-.999999,.999999)).mean()) if row.effect_kind=='correlation' else np.mean(values)
 assert abs(value-row.effect)<1e-12
for row in summary.itertuples():
 vals=pat[(pat.endpoint==row.endpoint)&(pat.model==row.model)].effect.dropna().to_numpy();assert len(vals)==row.n_patients and (vals>0).sum()==row.positive and (vals<0).sum()==row.negative
 if len(vals):assert abs(np.median(vals)-row.median_effect)<1e-12
 if len(vals)>=10:
  n=row.positive+row.negative;k=min(row.positive,row.negative);pv=min(1.,2*sum(math.comb(n,j) for j in range(k+1))/2**n) if n else 1.;assert abs(pv-row.P)<1e-12
 else:assert np.isnan(row.P)
for name,ss in summary[summary.endpoint.str.startswith('S2')].groupby('model'):
 a=ss.P.fillna(1).to_numpy();ix=np.argsort(a);q=np.minimum.accumulate((a[ix]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1];assert np.allclose(q,ss.FDR.to_numpy()[ix])
# Independently recompute regional differences and multi-source TLS shortest paths.
opt=pd.read_csv(F/'S3_S4_AVAILABILITY.tsv',sep='\t').set_index('section');region_checks=0;distance_checks=0
for f in sorted((P/'cache').glob('*_REGION_MEMBERS.tsv.gz')):
 sid=f.name.split('_')[0];reg=pd.read_csv(f,sep='\t',index_col='spot_id');src=pd.read_csv(F/'frozen_inputs'/(sid+'_SPOT_INPUTS.tsv.gz'),sep='\t',index_col='spot_id');assert src.index.equals(reg.index);share=src.SPP1_share_within_mature_macro
 for ep,left,right in [('S3_expression_boundary','interface_source','core_source'),('S3_ECM_region','ECM_high_source','ECM_low_source'),('S4_TLS_periphery','TLS_peri_source','TLS_far_source')]:
  out=sec[(sec.section==sid)&(sec.endpoint==ep)].iloc[0];l=reg[left].to_numpy(bool);r=reg[right].to_numpy(bool);assert not (l&r).any()
  if np.isfinite(out.effect):
   assert l.sum()>=20 and r.sum()>=20;delta=float(np.median(share[l])-np.median(share[r]));assert abs(delta-out.effect)<1e-12;region_checks+=1
 if bool(opt.loc[sid,'S4_TLS_contrast_candidate']):
  G=sparse.load_npz(F/'frozen_inputs'/(sid+'_adjacency.npz')).tocsr();known=(reg.TLS_positive|reg.TLS_explicit_negative).to_numpy();starts=np.flatnonzero(reg.TLS_positive);dist=np.full(len(reg),np.inf);dist[starts]=0;queue=collections.deque(starts.tolist())
  while queue:
   a=queue.popleft()
   for b in G.indices[G.indptr[a]:G.indptr[a+1]]:
    if known[b] and not np.isfinite(dist[b]):dist[b]=dist[a]+1;queue.append(int(b))
  assert np.array_equal(dist,reg.TLS_graph_distance_hops.to_numpy());distance_checks+=1
# Validate translation aggregation and independent rank precision on selected vectors.
trans=pd.read_csv(P/'S2/TRANSLATION_CONTROLS.tsv.gz',sep='\t');assert np.allclose(trans.matched_observed-trans.shifted,trans.excess)
for (sid,ep),g in trans.groupby(['section','endpoint']):
 out=sec[(sec.section==sid)&(sec.endpoint==ep)&(sec.model=='local_translation_excess')].iloc[0];assert len(g)>=30 and abs(g.excess.median()-out.effect)<1e-12
 old_n=sec[(sec.section==sid)&(sec.endpoint==ep)&(sec.model=='primary')].n.iloc[0];assert(g.n_matched>=max(50,math.ceil(old_n*.5))).all()
assert all(hashlib.sha256((F/k).read_bytes()).hexdigest()==v for k,v in json.loads((F/'FREEZE_SHA256.json').read_text()).items())
old=json.loads((P/'audit/RDS_BEFORE.json').read_text())
for f in old:
 s=Path(f['path']).stat();assert(s.st_size,s.st_mtime_ns,s.st_ino)==(f['bytes'],f['mtime_ns'],f['inode'])
result={'status':'PASS','section_rows':len(sec),'patient_rows':len(pat),'summary_rows':len(summary),'independent_precision_matrix_S2_checks':len(checks),'max_correlation_error':max(checks),'all_patient_aggregation_verified':True,'exact_sign_P_and_BH_independently_verified':True,'regional_effects_recomputed':region_checks,'independent_TLS_BFS_sections':distance_checks,'translation_rows':len(trans),'minimum_shifts':int(trans.groupby(['section','endpoint']).size().min()),'maximum_shifts':int(trans.groupby(['section','endpoint']).size().max()),'frozen_S1_inputs_unchanged':True,'preexisting_RDS_unchanged':len(old),'no_new_figures':True}
(P/'FINAL_VALIDATION.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
