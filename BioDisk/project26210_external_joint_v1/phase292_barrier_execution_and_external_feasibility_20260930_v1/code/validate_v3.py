from pathlib import Path
import hashlib,json,subprocess,datetime,itertools,math
import numpy as np,pandas as pd
from scipy import sparse
P=Path(__file__).resolve().parents[1];R=P.parent
checks=[]
def check(n,v,detail=''):
 checks.append(dict(check=n,passed=bool(v),detail=detail))
 if not v:raise AssertionError(n)
def read(p):return pd.read_csv(p,sep='\t',float_precision='round_trip')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def gate():
 q=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
 with (P/'logs/validation_memory.log').open('a') as f:f.write(q.stdout+q.stderr)
 check('memory_gate_'+str(len(checks)),q.returncode==0)
def independent(z,x,y,model):
 cols=[] if model=='study_only' else ['perivascular_fraction']
 if model.endswith('receiver_depth'):cols+=['log_receiver_depth']
 z=z.dropna(subset=[x,y]+cols);n=len(z)
 design=pd.get_dummies(z.study,dtype=float).to_numpy()
 if cols:design=np.column_stack([design,z[cols].rank(method='average').to_numpy()])
 U,S,V=np.linalg.svd(design,full_matrices=False);rank=int(np.linalg.matrix_rank(design));df=n-rank-1
 if rank<design.shape[1] or n-rank<2:return np.nan,np.nan,rank,df
 residual=np.eye(n)-U[:,:rank]@U[:,:rank].T
 a=residual@z[x].rank().to_numpy();b=residual@z[y].rank().to_numpy()
 den=np.linalg.norm(a)*np.linalg.norm(b)
 if den<1e-10:return np.nan,np.nan,rank,df
 rho=float(a.dot(b)/den)
 if n<10 or df<5:return rho,np.nan,rank,df
 groups=[np.where(z.study.to_numpy()==s)[0] for s in sorted(z.study.unique())]
 possible=math.prod(math.factorial(len(g)) for g in groups);N=min(possible,99999)
 order=np.repeat(np.arange(n)[None,:],N,axis=0)
 if possible<=99999:
  choices=[list(itertools.permutations(g)) for g in groups]
  for i,sets in enumerate(itertools.product(*choices)):
   for ids,v in zip(groups,sets):order[i,ids]=v
 else:
  rng=np.random.default_rng(292)
  for ids in groups:
   if len(ids)>1:order[:,ids]=ids[rng.random((N,len(ids))).argsort(axis=1)]
 null_vectors=b[order]@residual
 den=np.linalg.norm(null_vectors,axis=1)*np.linalg.norm(a)
 good=den>1e-10;null=null_vectors[good]@a/den[good]
 hits=int((np.abs(null)>=abs(rho)-1e-12).sum())
 pv=hits/len(null) if possible<=99999 else (hits+1)/(len(null)+1)
 return rho,pv,rank,df
def main():
 gate()
 f=json.loads((P/'FREEZE.json').read_text())
 check('protocol_hash',sha(P/'PLAN_AND_PROTOCOL.md')==f['protocol_sha256'])
 check('panels_hash',sha(P/'PANELS_FROZEN.json')==f['panels_sha256'])
 for r in json.loads((P/'qc/INPUT_MANIFEST.json').read_text()):check('input_hash_'+Path(r['path']).name,sha(r['path'])==r['sha256'])
 d=read(P/'tables/PATIENT_INPUTS.tsv');a=read(P/'tables/ASSOCIATIONS.tsv');panels=json.loads((P/'PANELS_FROZEN.json').read_text())
 check('unique_patient',len(d)==16 and d.patient.is_unique)
 check('eligibility_13',int(d.eligible.sum())==13 and np.array_equal(d.eligible,(d.stromal_cells>=50)&(d.malignant_cells>=50)&(d.macro_cells>=100)))
 m=read(R/'phase289_state_source_boundary_process_20260929_v1/tables/STROMAL_STATE_MEMBERSHIP.tsv.gz')
 original=read(R/'phase287_stromal_refinement_image_morphology_20260929_v1/inputs/STROMAL_INPUT_METADATA.tsv.gz')
 check('row_alignment',m.cell_id.tolist()==original.cell_id.tolist())
 X=sparse.load_npz(R/'phase287_stromal_refinement_image_morphology_20260929_v1/inputs/STROMAL_COUNTS.npz').tocsr().astype(np.float64)
 genes=np.load(R/'phase287_stromal_refinement_image_morphology_20260929_v1/inputs/GENES.npy',allow_pickle=True).tolist()
 cache=np.load(P/'inputs/PAIRED_STROMAL_PSEUDOBULK.npz');gene=read(P/'tables/PATIENT_GENE_COUNTS.tsv');sub=read(P/'tables/BROAD_CELLTYPE_DESCRIPTIVE.tsv')
 max_score=0.
 for i,r in enumerate(d.itertuples()):
  gate()
  select=m.curated_patient_id.eq(r.patient)&m.sample_id.eq(r.library)&m.source_tissue.eq('Tumor')&m.stromal_pool
  v=np.asarray(X[select.to_numpy()].sum(0)).ravel();total=float(v.sum())
  check('cell_and_UMI_'+r.patient,int(select.sum())==r.stromal_cells and total==r.total_UMI)
  check('cache_counts_'+r.patient,np.array_equal(v[[genes.index(g) for g in cache['genes']]],cache['counts'][i]))
  log=np.log1p(v/max(total,1)*1e4)
  for name,gs in panels.items():
   measured=[genes.index(g) for g in gs if g in genes]
   qualified=len(measured)>=.8*len(gs) and np.count_nonzero(v[measured]>0)>=math.ceil(len(measured)/2) and total>0
   value=float(np.mean(log[measured])) if qualified else np.nan
   delta=abs(value-getattr(r,name))
   if np.isfinite(delta):max_score=max(max_score,delta)
   check('score_'+r.patient+'_'+name,np.isclose(value,getattr(r,name),equal_nan=True,atol=1e-12,rtol=0))
  fraction=m.loc[select,'cell_type_broad_v1'].eq('Pericyte/SMC').mean()
  check('fraction_'+r.patient,np.isclose(fraction,r.perivascular_fraction,equal_nan=True))
  check('gene_counts_'+r.patient,np.array_equal(gene[gene.patient==r.patient].set_index('gene').loc[cache['genes'],'counts'],cache['counts'][i]))
  for celltype in ['Fibroblast','Pericyte/SMC']:
   ss=select&m.cell_type_broad_v1.eq(celltype);rr=sub[(sub.patient==r.patient)&(sub.celltype==celltype)].iloc[0]
   sv=np.asarray(X[ss.to_numpy()].sum(0)).ravel();st=float(sv.sum());sl=np.log1p(sv/max(st,1)*1e4)
   check('subtype_'+r.patient+'_'+celltype,int(ss.sum())==rr.cells and st==rr.total_UMI)
   if ss.sum():
    check('subtype_scores_'+r.patient+'_'+celltype,all(abs(np.mean([sl[genes.index(g)] for g in gs if g in genes])-rr[name])<1e-12 for name,gs in panels.items()))
 # Source values replayed from full original group count vectors.
 G=read(R/'phase157_ordered_bridge_reassessment_20260915_v1/02_inputs/GROUP_INDEX.tsv')
 Z=np.load(R/'phase157_ordered_bridge_reassessment_20260915_v1/02_inputs/GROUP_PSEUDOBULK.npz',allow_pickle=True)
 for r in d.itertuples():
  for label,roles in [('malignant',['Malignant']),('macro',['C4','OtherMacrophage'])]:
   hit=G.patient.eq(r.patient)&G.sample_id.eq(r.library)&G.source_tissue.eq('Tumor')&G.role.isin(roles)
   v=Z['counts'][G.loc[hit,'group_index'].to_numpy(int)].sum(0)
   for g in ['TYMP','SPP1']:
    expected=np.log1p(v[list(Z['genes']).index(g)]/v.sum()*1e4)
    check('source_'+r.patient+'_'+label+'_'+g,np.isclose(expected,getattr(r,label+'_'+g),atol=1e-12,rtol=0))
 gate()
 errors=[];perr=[]
 for i,r in a.iterrows():
  z=d[d.eligible]
  rho,pv,rank,df=independent(z,r.source,r.endpoint,r.model)
  check('model_rank_'+str(i),rank==r.design_rank and df==r.residual_df)
  check('rho_'+str(i),np.isclose(rho,r.rho,equal_nan=True,atol=1e-12,rtol=0))
  check('P_'+str(i),np.isclose(pv,r.P,equal_nan=True,atol=1e-12,rtol=0))
  if np.isfinite(rho):errors.append(abs(rho-r.rho))
  if np.isfinite(pv):perr.append(abs(pv-r.P))
 for key,z in a.groupby('family'):
  pv=z.P.fillna(1).to_numpy();order=np.argsort(pv);q=np.minimum.accumulate((pv[order]*len(pv)/np.arange(1,len(pv)+1))[::-1])[::-1].clip(0,1);back=np.empty(len(pv));back[order]=q
  back[z.P.isna().to_numpy()]=np.nan
  check('BH_'+key,np.allclose(back,z.FDR,equal_nan=True,atol=1e-12))
 check('all_models_27',len(a)==27)
 check('no_cell_pseudorep',a.n.eq(13).all())
 check('constant_mural_fraction_all_qualified',d.loc[d.eligible,'perivascular_fraction'].eq(1).all())
 check('singular_models_kept',a[a.model.ne('study_only')].P.isna().all() and a[a.model.ne('study_only')].method.eq('DESIGN_RANK_OR_RESIDUAL_DIM').all())
 for rr in json.loads((P/'qc/INPUT_MANIFEST.json').read_text()):check('final_hash_'+Path(rr['path']).name,sha(rr['path'])==rr['sha256'])
 pd.DataFrame(checks).to_csv(P/'qc/VALIDATION_CHECKS.tsv',sep='\t',index=False)
 result=dict(status='PASS',checks=len(checks),failed=0,max_score_error=max_score,max_rho_error=max(errors),max_P_error=max(perr),at=datetime.datetime.now().astimezone().isoformat(),interpretation='valid execution; primary composition-adjusted models not estimable')
 (P/'qc/INDEPENDENT_VALIDATION.json').write_text(json.dumps(result,indent=2))
 print(json.dumps(result))
if __name__=='__main__':
 try:main()
 except Exception:
  pd.DataFrame(checks).to_csv(P/'qc/VALIDATION_CHECKS_FAILED_v3.tsv',sep='\t',index=False)
  raise
