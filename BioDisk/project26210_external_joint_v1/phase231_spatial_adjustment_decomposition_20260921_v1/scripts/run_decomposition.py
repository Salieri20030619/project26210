from pathlib import Path
import ast,json,hashlib,traceback
import numpy as np,pandas as pd
from scipy import stats,sparse
P=Path(__file__).resolve().parents[1];R=P.parent;F=R/'phase229_spatial_F1_partition_eligibility_20260921_v1';O=R/'phase230_spatial_F2_program_association_20260921_v1'
plan=json.loads((P/'PLAN_BEFORE_COMPUTE.json').read_text());models=plan['models']
tree=ast.parse((O/'scripts/run_f2.py').read_text());func=[x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='infer'];exec(compile(ast.Module(body=func,type_ignores=[]),'<inherited_infer>','exec'))
def save(n,o):(P/n).write_text(json.dumps(o,indent=2))
def tab(n,d):d.to_csv(P/n,sep='\t',index=False)
def run():
 elig=pd.read_csv(F/'tables/SECTION_ELIGIBILITY.tsv',sep='\t');records=[];checks=[];hashes=[]
 for r in elig.itertuples():
  path=F/'spots'/(r.sample+'_F1_SPATIAL_INPUTS.tsv.gz');d=pd.read_csv(path,sep='\t');hashes.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
  base=R/'phase221_spatial_S1_identity_eligibility_freeze_20260920_v1/frozen_inputs' if r.dataset=='GSE175540' else R/'phase223_spatial_cohort_reinclusion_20260920_v1/inputs'
  G=sparse.load_npz(base/(r.sample+'_adjacency.npz'));deg=np.asarray(G.sum(axis=1)).ravel();W=sparse.diags(1/np.maximum(deg,1))@G
  mask=d.source_eligible.to_numpy(bool);assert mask.sum()>=50
  x=d.SPP1_share.to_numpy()
  for target in plan['targets']:
   source=d['program_'+target].to_numpy();y=W@source
   vals={'x':x,'y':y,'source_program':source}
   for c in set(sum(models.values(),[]))-{'source_program'}:vals[c]=d[c].to_numpy()
   ranks={k:stats.rankdata(v[mask]) for k,v in vals.items()}
   xy=np.column_stack([ranks['x'],ranks['y']])
   for model,cols in models.items():
    if cols:
     cov=np.column_stack([ranks[c] for c in cols]);sd=cov.std(axis=0);assert(sd>0).all()
     D=np.column_stack([np.ones(mask.sum()),(cov-cov.mean(axis=0))/sd])
    else:D=np.ones((int(mask.sum()),1))
    rank=int(np.linalg.matrix_rank(D));condition=float(np.linalg.cond(D))
    assert rank==D.shape[1] and condition<30,(r.sample,target,model,condition)
    res=xy-D@np.linalg.lstsq(D,xy,rcond=None)[0]
    rho=float(np.corrcoef(res.T)[0,1])
    Q=np.linalg.qr(D,mode='reduced')[0];res2=xy-Q@(Q.T@xy);rho2=float(np.corrcoef(res2.T)[0,1])
    assert abs(rho-rho2)<1e-10
    if model=='M0_unadjusted':assert abs(stats.spearmanr(x[mask],y[mask]).statistic-rho)<1e-12
    records.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,target=target,model=model,n_source=int(mask.sum()),rho=rho,rank=rank,condition=condition))
    checks.append(abs(rho-rho2))
 sections=pd.DataFrame(records);assert len(sections)==608
 tab('tables/SECTION_MODELS.tsv',sections)
 prior=pd.read_csv(O/'tables/SECTION_EFFECTS.tsv',sep='\t');prior=prior[prior.model=='primary']
 full=sections[sections.model=='M3_full'].merge(prior,on=['sample','target'],suffixes=('_new','_old'));err=float(abs(full.rho_new-full.rho_old).max());assert len(full)==76 and err<1e-12
 patients=[]
 for (dataset,patient,target,model),g in sections.groupby(['dataset','patient','target','model']):
  rho=float(np.tanh(np.arctanh(g.rho.to_numpy()).mean()))
  patients.append(dict(dataset=dataset,patient=patient,target=target,model=model,rho=rho,n_sections=len(g),sections=';'.join(g['sample'])))
 pat=pd.DataFrame(patients);tab('tables/PATIENT_MODELS.tsv',pat)
 summary=[]
 for i,((dataset,target,model),g) in enumerate(pat.groupby(['dataset','target','model'])):
  summary.append(dict(dataset=dataset,target=target,model=model,**infer(g.rho.to_numpy(),263100+i)))
 su=pd.DataFrame(summary);su['FDR_within_model']=np.nan;su['FDR_diagnostic16']=np.nan
 for (dataset,model),g in su.groupby(['dataset','model']):
  if g.P.notna().all():su.loc[g.index,'FDR_within_model']=stats.false_discovery_control(g.P.to_numpy())
 for dataset,g in su.groupby('dataset'):
  if g.P.notna().all():su.loc[g.index,'FDR_diagnostic16']=stats.false_discovery_control(g.P.to_numpy())
 su['nominal_positive']=(su.P<.05)&(su.median_rho>0)
 su['nominal_negative']=(su.P<.05)&(su.median_rho<0)
 tab('tables/COHORT_MODELS.tsv',su)
 # All12 edges of the three-block cube: add one block while holding others fixed.
 modelsets={m:{k for k,cols in plan['blocks'].items() if set(cols).issubset(set(c))} for m,c in models.items()}
 diffs=[]
 for (dataset,target),g in pat.groupby(['dataset','target']):
  wide=g.pivot(index='patient',columns='model',values='rho');assert wide.notna().all().all()
  for m0,b0 in modelsets.items():
   for m1,b1 in modelsets.items():
    if b0.issubset(b1) and len(b1-b0)==1:
     dv=(wide[m1]-wide[m0]).to_numpy();n=len(dv);low=high=np.nan
     if n>=10:
      rng=np.random.default_rng(263199);boot=np.median(dv[rng.integers(0,n,(10000,n))],axis=1);low,high=np.quantile(boot,[.025,.975])
     diffs.append(dict(dataset=dataset,target=target,from_model=m0,to_model=m1,added_block=next(iter(b1-b0)),n=n,median_before=wide[m0].median(),median_after=wide[m1].median(),paired_median_delta=np.median(dv),positive_delta=int((dv>0).sum()),negative_delta=int((dv<0).sum()),delta_CI_low=low,delta_CI_high=high))
 tab('tables/PAIRED_BLOCK_EFFECT_CHANGES.tsv',pd.DataFrame(diffs))
 for x in hashes:assert hashlib.sha256(Path(x['path']).read_bytes()).hexdigest()==x['sha256']
 for x in json.loads((F/'audit/RDS_BEFORE.json').read_text()):
  st=Path(x['path']).stat();assert(st.st_size,st.st_mtime_ns,st.st_ino)==(x['bytes'],x['mtime_ns'],x['inode'])
 save('INPUT_SHA256.json',hashes)
 save('COMPUTE_VALIDATION.json',dict(status='PASS',models=608,QR_max_error=max(checks),previous_full_model_max_error=err,same_source_spots=True,all_models_estimable=True,old_RDS_unchanged=333,input_hashes_unchanged=True))
 save('STATUS.json',dict(status='COMPUTE_COMPLETE',new_figures=0));(P/'EXIT_STATUS').write_text('0\n')
 print(su[su.model.isin(plan['main_order'])].to_string(index=False),flush=True)
try:run()
except Exception:
 (P/'ERROR.txt').write_text(traceback.format_exc());(P/'EXIT_STATUS').write_text('1\n');raise
