from pathlib import Path
import json,datetime,hashlib,ast,traceback,subprocess,gc,platform
import numpy as np,pandas as pd,h5py,scipy
from scipy import sparse,stats
from anndata.io import read_elem
P=Path(__file__).resolve().parents[1];R=P.parent;OLD=R/'phase234_part4_replan_execution_20260921_v1';GEO=R/'phase276_interface_topology_and_localization_20260928_v1';MARK=R/'phase246_spatial_marker_structure_20260924_v1'
def now():return datetime.datetime.now().astimezone().isoformat()
def dump(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
def tab(d,f):pd.DataFrame(d).to_csv(P/'tables'/f,sep='\t',index=False)
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def status(s,**kw):dump(P/'STATUS.json',dict(status=s,at=now(),**kw))
def gate(stage):
 q=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
 with open(P/'logs/memory_gate.log','a') as f:f.write(stage+'\n'+q.stdout+q.stderr)
 if q.returncode:raise RuntimeError('MEMORY_GATE')
for file,names in [(OLD/'scripts/common.py',['fit','infer']),(GEO/'code/run_tests.py',['neighborhoods'])]:
 tree=ast.parse(file.read_text());exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names],type_ignores=[]),'inherited','exec'),globals())
CORE=['depth_mean','depth_delta','macro_logmean','malignant_mean','malignant_delta','log_side_count_ratio']
KEYS=['dataset','family','endpoint','model','effect_type']
def summarize(rows):
 sec=pd.DataFrame(rows);tab(sec,'SECTION_EFFECTS.tsv');out=[]
 for k,g in sec.groupby(KEYS+['patient'],sort=True):
  v=g.effect.dropna().to_numpy()
  e=(float(np.tanh(np.arctanh(np.clip(v,-.999999999,.999999999)).mean())) if k[4]=='rho' else float(v.mean())) if len(v) else np.nan
  out.append(dict(zip(KEYS+['patient'],k),effect=e,n_sections=len(v),sections=';'.join(g.loc[g.effect.notna(),'sample'])))
 pat=pd.DataFrame(out);tab(pat,'PATIENT_EFFECTS.tsv');out=[]
 for i,(k,g) in enumerate(pat.groupby(KEYS,sort=True)):out.append(dict(zip(KEYS,k),**infer(g.effect,26210277+i)))
 su=pd.DataFrame(out)
 for _,g in su.groupby(['dataset','family']):
  for pc,qc in [('P','FDR'),('P_wilcoxon','FDR_wilcoxon')]:
   su.loc[g.index,qc]=np.where(g[pc].notna(),stats.false_discovery_control(g[pc].fillna(1).to_numpy()),np.nan)
  su.loc[g.index,'family_size_registered']=len(g);su.loc[g.index,'family_size_estimable']=g.P.notna().sum()
 su['evidence']=np.where(su.P.isna(),'DESCRIPTIVE',np.where(su.P<.05,np.where(su.FDR<.05,'FDR_SIGNIFICANT','NOMINAL_SIGNIFICANT'),'NOT_SUPPORTED'))
 tab(su,'COHORT_SUMMARY.tsv');loo=[]
 for k,g in pat[pat.family=='C4_primary_interaction'].groupby(KEYS):
  g=g[g.effect.notna()]
  if len(g)<10:continue
  for row in g.itertuples():
   v=g.loc[g.patient!=row.patient,'effect'].to_numpy();nz=np.count_nonzero(v)
   loo.append(dict(zip(KEYS,k),omitted_patient=row.patient,n=len(v),median=float(np.median(v)),positive=int((v>0).sum()),negative=int((v<0).sum()),P=stats.binomtest(int((v>0).sum()),nz,.5).pvalue if len(v)>=10 and nz else np.nan))
 tab(loo,'C4_LEAVE_ONE_PATIENT_OUT.tsv')
 return sec,pat,su
def interaction(z,M,F,target,covnames):
 cols=[M,F,target]+covnames;mask=z.source_eligible.to_numpy(bool)&np.isfinite(z[cols]).all(1).to_numpy()
 a=z.loc[mask,cols].to_numpy();n=len(a)
 out=dict(effect=np.nan,n=n,reason='',condition=np.nan,beta_M=np.nan,beta_F=np.nan,beta_MF=np.nan,intercept=np.nan,F_slope_M_low=np.nan,F_slope_M_high=np.nan,quadrant_min=0)
 if n<50:out['reason']='FEWER_THAN_50_SEGMENTS';return out
 sd=a.std(0)
 if np.any(sd==0):out['reason']='CONSTANT_PREDICTOR_OR_OUTCOME';return out
 x=(a[:,0]-a[:,0].mean())/sd[0];f=(a[:,1]-a[:,1].mean())/sd[1]
 q=[int((((x>=np.median(x))==mx)&((f>=np.median(f))==ff)).sum()) for mx in [False,True] for ff in [False,True]]
 out['quadrant_min']=min(q)
 if min(q)<5:out['reason']='QUADRANT_SUPPORT_LT5';return out
 cv=(a[:,3:]-a[:,3:].mean(0))/sd[3:]
 X=np.column_stack([np.ones(n),x,f,x*f,cv]);pred=X[:,1:];ds=pred.std(0)
 if (ds==0).any():out['reason']='CONSTANT_DESIGN';return out
 normalized=np.column_stack([np.ones(n),(pred-pred.mean(0))/ds]);cond=float(np.linalg.cond(normalized));out['condition']=cond
 if np.linalg.matrix_rank(X)<X.shape[1] or cond>=30:out['reason']='RANK_OR_CONDITION';return out
 b=np.linalg.lstsq(X,a[:,2],rcond=None)[0];assert np.isfinite(b).all()
 out.update(effect=float(b[3]),beta_M=float(b[1]),beta_F=float(b[2]),beta_MF=float(b[3]),intercept=float(b[0]),F_slope_M_low=float(b[2]-b[3]),F_slope_M_high=float(b[2]+b[3]))
 return out
def identity(programs):
 f=R/'phase200_union_C4_C10_C15_reassessment_20260918_v1/core/contrasts/tumor_C4.npz'
 z=np.load(f,allow_pickle=True);gg=list(z['genes']);lookup={g:i for i,g in enumerate(gg)};rows=[];sm=[];study=[]
 panels=dict(APOE=['APOE'],**{k:programs[k] for k in ['APOEMph_noC4','cholesterol_disjoint','lysosome_disjoint']})
 for idx,(name,genes) in enumerate(panels.items()):
  present=[g for g in genes if g in lookup];assert len(present)>=.8*len(genes);ii=[lookup[g] for g in present]
  a=z['target_log2CPM'][:,ii].mean(1);b=z['reference_log2CPM'][:,ii].mean(1);v=a-b
  for j in range(len(v)):rows.append(dict(patient=str(z['patients'][j]),study=str(z['study'][j]),program=name,C4_value=float(a[j]),other_macro_value=float(b[j]),effect=float(v[j]),genes_used=';'.join(present)))
  sm.append(dict(program=name,**infer(v,26210770+idx)))
  for st in sorted(set(z['study'])):
   vv=v[z['study']==st];study.append(dict(program=name,study=str(st),n=len(vv),positive=int((vv>0).sum()),negative=int((vv<0).sum()),median=float(np.median(vv))))
 su=pd.DataFrame(sm)
 for pc,qc in [('P','FDR'),('P_wilcoxon','FDR_wilcoxon')]:su[qc]=stats.false_discovery_control(su[pc])
 su['family_size']=4;tab(rows,'C4_IDENTITY_PATIENT_EFFECTS.tsv');tab(su,'C4_IDENTITY_SUMMARY.tsv');tab(study,'C4_IDENTITY_STUDY_DIRECTIONS.tsv')
 old=pd.read_csv(R/'phase200_union_C4_C10_C15_reassessment_20260918_v1/core/contrasts/tumor_C4_genes.tsv.gz',sep='\t')
 tab(old[old.gene.isin(['SPP1','APOE','CYP27A1','GPNMB','LPL','FABP5'])],'INHERITED_C4_MARKER_CONTEXT.tsv')
 old=pd.read_csv(R/'phase200_union_C4_C10_C15_reassessment_20260918_v1/functional/FOCUS_PROGRAM_RESULTS.tsv',sep='\t')
 tab(old[(old.contrast=='tumor_C4')&old.term.str.contains('Cholesterol Homeostasis|^Lysosome$',regex=True)],'INHERITED_C4_FUNCTIONAL_CONTEXT.tsv')
 return dict(path=str(f),sha256=sha(f))

def main():
 gate('start');freeze=json.loads((P/'FREEZE.json').read_text())
 assert sha(P/'PROTOCOL_BEFORE_TESTING.md')==freeze['protocol_sha256'] and sha(P/'PROGRAMS_FROZEN.json')==freeze['programs_sha256']
 programs=json.loads((P/'PROGRAMS_FROZEN.json').read_text())['programs'];hashes=[identity(programs)];rawmeta=[]
 regs=pd.read_csv(OLD/'tables/INPUT_SECTION_REGISTER.tsv',sep='\t');regs.to_csv(P/'tables/INPUT_SECTION_REGISTER.tsv',sep='\t',index=False)
 effects=[];eligibility=[];qc=[];design=[]
 for i,r in enumerate(regs.itertuples()):
  gate(r.sample);status('RUNNING',completed=i,total=38,current=r.sample)
  files=[OLD/'inputs'/f'{r.sample}_INPUTS.tsv.gz',OLD/'graphs'/f'{r.sample}_distance800.npz',GEO/'inputs'/f'{r.sample}_nodes.tsv.gz',GEO/'inputs'/f'{r.sample}_membership.tsv.gz',GEO/'inputs'/f'{r.sample}_segments.tsv.gz']
  hashes.extend(dict(path=str(f),sha256=sha(f)) for f in files)
  d=pd.read_csv(files[0],sep='\t');D=sparse.load_npz(files[1]).tocsr();nodes=pd.read_csv(files[2],sep='\t',float_precision='round_trip');mem=pd.read_csv(files[3],sep='\t',float_precision='round_trip');z=pd.read_csv(files[4],sep='\t',float_precision='round_trip')
  f=Path(r.input_path);st=f.stat();rawmeta.append(dict(path=str(f),bytes=st.st_size,mtime_ns=st.st_mtime_ns))
  inherited=json.loads((P/'PROGRAMS_FROZEN.json').read_text())['existing_programs']
  checkpanels={'check_POSTNCAF':inherited['POSTNCAF'],'check_ECM4':inherited['ECM4']}
  allgenes=sorted(set(sum(list(programs.values())+list(checkpanels.values()),[])))
  with h5py.File(f,'r') as h:
   obs=read_elem(h['obs']);var=read_elem(h['var']);assert np.array_equal(obs.index,d.spot_id) and var.index.is_unique
   X=read_elem(h['layers/counts'] if 'counts' in h.get('layers',{}) else h['X']).astype(float)
   assert sparse.issparse(X) and np.isfinite(X.data).all() and (X.data>=0).all()
   total=np.asarray(X.sum(1)).ravel();assert (total>0).all()
   present=[g for g in allgenes if g in var.index];raw=X[:,var.index.get_indexer(present)].toarray()
  log=np.log1p(raw/total[:,None]*10000);lookup={g:j for j,g in enumerate(present)};scores={};maxerr=0.
  for name,panel in {**programs,**checkpanels}.items():
   ids=[lookup[g] for g in panel if g in lookup];vv=log[:,ids]
   if len(ids):
    sd=vv.std(0);score=((vv-vv.mean(0))/np.where(sd>0,sd,1)).mean(1);det=(raw[:,ids]>0).sum(1)>= (3 if name.startswith('C4') else 2)
   else:score=np.full(len(d),np.nan);det=np.zeros(len(d),bool)
   coverage=(len(ids)==len(panel)) if name.startswith('C4') else len(ids)>=max(3,int(np.ceil(.8*len(panel))))
   if name.startswith('check_'):coverage=len(ids)>=max(3,int(np.ceil(.5*len(panel))))
   ok=coverage and det.sum()>=20 and np.nanstd(score)>0
   if name.startswith('check_'):
    original=d[name.replace('check_','')].to_numpy()
    if ok and np.isfinite(original).all():
     err=float(np.max(np.abs(score-original)));assert err<1e-9;maxerr=max(maxerr,err)
    continue
   scores[name]=score if ok else np.full(len(d),np.nan)
   eligibility.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,program=name,genes_required=len(panel),genes_present=len(ids),present_genes=';'.join(g for g in panel if g in lookup),missing_genes=';'.join(g for g in panel if g not in lookup),detected_spots=int(det.sum()),eligible=ok))
  pd.DataFrame(dict(spot_id=d.spot_id,**scores)).to_csv(P/'inputs'/f'{r.sample}_PROGRAMS.tsv.gz',sep='\t',index=False)
  end=nodes[['spot_a','spot_b']].to_numpy(int);ns,ds,_,_=neighborhoods(end,D,len(d))
  source=[ids[dist<=200+1e-6] for ids,dist in zip(ns,ds)]
  z['source_n']=0;z['source_macro_fraction']=np.nan
  sourcecols=list(programs)+['POSTNCAF','ECM4','SPP1_share','TYMP_RNA','source_endothelial']
  for col in sourcecols:z[col]=np.nan
  idsrows=[]
  for sid,gg in mem.groupby('segment',sort=True):
   idx=np.unique(np.concatenate([source[j] for j in gg.node.to_numpy(int)]))
   z.at[sid,'source_n']=len(idx);z.at[sid,'source_macro_fraction']=float(d.total_macro_fraction.iloc[idx].mean())
   for col,values in scores.items():z.at[sid,col]=float(values[idx].mean())
   for col in ['POSTNCAF','ECM4','TYMP_RNA']:z.at[sid,col]=float(d[col].iloc[idx].mean())
   z.at[sid,'SPP1_share']=float(d.V1_SPP1_q05.iloc[idx].sum()/d.total_macro_q05.iloc[idx].sum())
   z.at[sid,'source_endothelial']=float(d.endothelial_fraction.iloc[idx].mean())
   idsrows.extend(dict(segment=int(sid),spot_index=int(j)) for j in idx)
  z['source_eligible']=z.eligible&(z.source_n>=10)&(z.source_macro_fraction>=.05)
  z.to_csv(P/'inputs'/f'{r.sample}_MODEL_INPUT.tsv.gz',sep='\t',index=False)
  pd.DataFrame(idsrows).to_csv(P/'inputs'/f'{r.sample}_SOURCE_MEMBERSHIP.tsv.gz',sep='\t',index=False)
  def add(family,ep,model,e,n,why='',typ='rho'):
   effects.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,family=family,endpoint=ep,model=model,effect_type=typ,effect=e,eligible_segments=n,reason=why))
  for sourcevar in ['APOEMph_noC4','SPP1_share']:
   for model in ['unadjusted','core_adjusted']:
    e,n,why=fit(z[sourcevar].to_numpy(),z.POSTNCAF.to_numpy(),cov=None if model=='unadjusted' else z[CORE].to_numpy(),mask=z.source_eligible.to_numpy(bool),min_n=50)
    add('process_association',sourcevar+'__POSTNCAF',model,e,n,why)
  for sourcevar in ['cholesterol_disjoint','polyamine_MTA_context','TYMP_RNA']:
   for model in ['unadjusted','core_adjusted']:
    e,n,why=fit(z[sourcevar].to_numpy(),z.C4_19.to_numpy(),cov=None if model=='unadjusted' else z[CORE].to_numpy(),mask=z.source_eligible.to_numpy(bool),min_n=50)
    add('metabolic_context',sourcevar+'__C4_19',model,e,n,why)
  configs=[]
  for M in ['APOEMph_noC4','SPP1_share']:
   for model,cov in [('core_adjusted',CORE),('vascular_adjusted',CORE+['source_endothelial'])]:configs.append(('process_interaction',M,'POSTNCAF',model,cov))
  for M,family in [('C4_19','C4_primary_interaction'),('C4_17','C4_disjoint_sensitivity')]:
   for F in ['POSTNCAF','ECM4']:
    for model,cov in [('core_adjusted',CORE),('state_adjusted',CORE+['SPP1_share','APOEMph_noC4','source_endothelial'])]:configs.append((family,M,F,model,cov))
  for family,M,F,model,cov in configs:
   for target in ['CD8_RNA','T_abundance_log2']:
    o=interaction(z,M,F,target+'_delta',cov);ep=M+'x'+F+'__'+target
    add(family,ep,model,o['effect'],o['n'],o['reason'],'interaction_beta')
    design.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,family=family,endpoint=ep,model=model,macrophage=M,stroma=F,target=target,covariates=';'.join(cov),**o))
  qc.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,spots=len(d),all_segments=len(z),geometric_eligible=int(z.eligible.sum()),source_eligible=int(z.source_eligible.sum()),C4_19_qualified=bool(np.isfinite(scores['C4_19']).all()),C4_17_qualified=bool(np.isfinite(scores['C4_17']).all()),marker_reconstruction_error=maxerr))
  print('SECTION_COMPLETE',i+1,r.sample,'source_segments',int(z.source_eligible.sum()),flush=True)
  del X,D,raw,log;gc.collect()
 tab(eligibility,'PROGRAM_ELIGIBILITY.tsv');tab(qc,'SECTION_QC.tsv');tab(design,'INTERACTION_MODEL_DETAILS.tsv')
 sec,pat,su=summarize(effects)
 for row in hashes:assert sha(row['path'])==row['sha256']
 for row in rawmeta:
  st=Path(row['path']).stat();assert st.st_size==row['bytes'] and st.st_mtime_ns==row['mtime_ns']
 dump(P/'qc/READONLY_INPUTS.json',dict(derived_sha256=hashes,raw_h5ad_metadata=rawmeta))
 dump(P/'qc/SOFTWARE.json',dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,h5py=h5py.__version__))
 dump(P/'qc/COMPUTATION_COMPLETED.json',dict(at=now(),sections=38,patients=37,section_effects=len(sec),patient_effects=len(pat),cohort_endpoints=len(su),identity_endpoints=4,figures=0))
 status('COMPUTED_AWAITING_VALIDATION',completed=38,total=38);print('ALL_COMPUTATION_COMPLETED',len(su),flush=True)
if __name__=='__main__':
 try:main()
 except Exception:status('FAILED',traceback=traceback.format_exc());raise
