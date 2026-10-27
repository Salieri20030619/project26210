from pathlib import Path
import os,sys,json,datetime,hashlib,subprocess,traceback,ast,gc
import numpy as np,pandas as pd,h5py
from scipy import sparse,stats
from scipy.sparse.csgraph import dijkstra
from anndata.io import read_elem
P=Path(__file__).resolve().parents[1];ROOT=P.parent;OLD=ROOT/'phase234_part4_replan_execution_20260921_v1'
def now():return datetime.datetime.now().astimezone().isoformat()
def dump(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
def tab(d,f):d.to_csv(P/'tables'/f,sep='\t',index=False)
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def gate(stage):
 r=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
 with open(P/'logs/memory_gate.log','a') as f:f.write(stage+'\n'+r.stdout+r.stderr)
 if r.returncode:raise RuntimeError('memory blocked '+stage)
def state(s,**kw):dump(P/'STATUS.json',dict(status=s,at=now(),**kw))
tree=ast.parse((OLD/'scripts/common.py').read_text())
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['fit','infer']],type_ignores=[]),'inherited_numeric_functions','exec'),globals())
def summarize(rows,prefix):
 d=pd.DataFrame(rows);tab(d,prefix+'_SECTION_EFFECTS.tsv')
 keys=['dataset','family','endpoint','model','effect_type'];out=[]
 for key,g in d.groupby(keys+['patient'],dropna=False,sort=True):
  v=g.effect.dropna().to_numpy()
  e=(float(np.tanh(np.arctanh(np.clip(v,-.999999999,.999999999)).mean())) if key[4]=='rho' else float(np.mean(v))) if len(v) else np.nan
  out.append(dict(zip(keys+['patient'],key),effect=e,n_sections=len(v)))
 patients=pd.DataFrame(out);tab(patients,prefix+'_PATIENT_EFFECTS.tsv');s=[]
 for i,(key,g) in enumerate(patients.groupby(keys,dropna=False,sort=True)):s.append(dict(zip(keys,key),**infer(g.effect,26210246+i)))
 su=pd.DataFrame(s);su['FDR']=np.nan
 for k,g in su.groupby(['dataset','family']):
  q=stats.false_discovery_control(g.P.fillna(1).to_numpy())
  su.loc[g.index,'FDR']=np.where(g.P.notna(),q,np.nan)
  su.loc[g.index,'family_size']=len(g)
 tab(su,prefix+'_COHORT_SUMMARY.tsv')
 return su
def main():
 gate('freeze')
 markers={'Macrophage':['CD68','CSF1R','FCER1G','AIF1'],'Tcell':['CD3D','CD3E','TRAC'],'CD8lineage':['CD3D','CD3E','TRAC','CD8A','CD8B'],'Endothelial':['PECAM1','VWF','EMCN','KDR'],'Fibroblast':['COL1A1','COL1A2','COL3A1','DCN']}
 plan=dict(at=now(),known_previous_results=True,posthoc_exploratory=True,markers=markers,identity='RNA-marker concordance with existing mapping, not independent protein/cell identity; markers not held out during original mapping',marker_gate='all listed genes present; >=20 spots with >=2 marker genes detected; CD8 target zones additionally require CD8B and >=1 T-lineage transcript',neighborhood='fixed400um graph excludes self; source mature-macrophage fraction>=0.05, >=3 neighbors, >=50 sources; old V3 SPP1-state only',models=['A0','A3_rank_depth_composition_source_target'],A3_covariates='source+neighbor log depth, source+neighbor mature-macrophage fraction, source+neighbor malignant fraction, source target score',distance='SPP1 macro share top20% minus bottom20% within macro-rich spots; graph distance to target marker top20% among marker-qualified target spots; >=20 targets, >=20 reachable sources per group; no cross-component bridges',patient='same-patient sections Fisher-z mean for rho or arithmetic mean for distance; GLMF2 descriptive; GSE17554023 and TJ12 separate; n>=10 exact sign primary, inherited Wilcoxon companion, bootstrap10000',FDR='per dataset and module family, all models and fixed endpoints included; no endpoint-based sample selection',composite='original 4x4/5x5 block medians; original T_balance, balance excluding SPP1 from reference, direct log1p T-cell abundance; all three and both grids in one6 family. SPP1-only versus vascular compartment separate2 family',old_results='read only',new_figures=0)
 dump(P/'PLAN_BEFORE_ANALYSIS.json',plan)
 regs=pd.read_csv(OLD/'tables/INPUT_SECTION_REGISTER.tsv',sep='\t')
 assert len(regs)==38 and regs['sample'].is_unique
 genes=sorted(set(sum(markers.values(),[])))
 effects=[];coverage=[];inputs=[];geometry=[]
 for i,r in enumerate(regs.itertuples()):
  gate('sample_'+r.sample);state('MARKER_AND_STRUCTURE',completed=i,total=38,current=r.sample)
  f=Path(r.input_path);info=f.stat();inputs.append(dict(path=str(f),bytes=info.st_size,mtime_ns=info.st_mtime_ns))
  d=pd.read_csv(OLD/'inputs'/f'{r.sample}_INPUTS.tsv.gz',sep='\t')
  with h5py.File(f,'r') as h:
   obs=read_elem(h['obs']);var=read_elem(h['var']);assert np.array_equal(obs.index,d.spot_id) and var.index.is_unique
   X=read_elem(h['layers/counts'] if 'counts' in h.get('layers',{}) else h['X']).astype(float)
   assert sparse.issparse(X) and np.isfinite(X.data).all() and (X.data>=0).all()
   tot=np.asarray(X.sum(1)).ravel();assert (tot>0).all()
   present=[g for g in genes if g in var.index];counts=X[:,var.index.get_indexer(present)].toarray()
   log=np.log1p(counts/tot[:,None]*10000);lookup={g:i for i,g in enumerate(present)}
  scores={};detected={};ok={}
  for name,m in markers.items():
   used=[g for g in m if g in lookup]
   jj=[lookup[g] for g in used]
   if jj:
    raw=counts[:,jj];v=log[:,jj];sd=v.std(0);s=((v-v.mean(0))/np.where(sd>0,sd,1)).mean(1);detect=(raw>0).sum(1)>=2
   else:s=np.full(len(d),np.nan);detect=np.zeros(len(d),bool)
   if name=='CD8lineage':
    detect &= (counts[:,lookup['CD8B']]>0) if 'CD8B' in lookup else False
    detect &= (counts[:,[lookup[g] for g in markers['Tcell'] if g in lookup]]>0).any(1)
   qualified=len(used)==len(m) and detect.sum()>=20 and np.nanstd(s)>0
   scores[name]=s if qualified else np.full(len(d),np.nan);detected[name]=detect;ok[name]=qualified
   coverage.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,marker=name,present=';'.join(used),missing=';'.join(g for g in m if g not in lookup),detected_spots=int(detect.sum()),eligible=qualified,independent_cell_identity=False))
  saved=pd.DataFrame({'spot_id':d.spot_id,**scores});saved.to_csv(P/'inputs'/f'{r.sample}_MARKERS.tsv.gz',sep='\t',index=False)
  def add(fam,ep,model,e,n,reason='',typ='rho'):
   effects.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,family=fam,endpoint=ep,model=model,effect_type=typ,effect=e,eligible_spots=n,reason=reason))
  for target,ab in [('Macrophage','total_macro_fraction'),('Tcell','Tcell_fraction'),('Endothelial','endothelial_fraction')]:
   e,n,reason=fit(d[ab],scores[target]);add('mapping_RNA_concordance',ab+'__'+target,'A0',e,n,reason)
  e,n,reason=fit(d.V3_SPP1_macro_share,d.SPP1_RNA,mask=d.total_macro_fraction.to_numpy()>=.05)
  add('mapping_RNA_concordance','SPP1_macro_share__SPP1_RNA','macro_ge005',e,n,reason)
  D=sparse.load_npz(OLD/'graphs'/f'{r.sample}_distance800.npz').tocsr();assert D.shape==(len(d),len(d))
  H=D.copy();H.data=(H.data<=400).astype(float);H.eliminate_zeros();nn=np.asarray(H.sum(1)).ravel();W=sparse.diags(1/np.maximum(nn,1))@H
  G=D.copy();G.data=np.where(np.isclose(G.data,100),100.,0.);G.eliminate_zeros()
  assert not G.diagonal().any() and (G-G.T).nnz==0
  src=(d.total_macro_fraction.to_numpy()>=.05)&(nn>=3)
  for target in ['CD8lineage','Endothelial','Fibroblast']:
   y=scores[target];ny=W@y
   for model in ['A0','A3']:
    cov=None if model=='A0' else np.column_stack([d.log1p_UMI,W@d.log1p_UMI.to_numpy(),d.total_macro_fraction,W@d.total_macro_fraction.to_numpy(),d.malignant_fraction,W@d.malignant_fraction.to_numpy(),y])
    e,n,reason=fit(d.V3_SPP1_macro_share,ny,cov=cov,mask=src)
    add('gene_disjoint_neighborhood','SPP1_share__'+target+'_400um',model,e,n,reason)
   mask=(d.total_macro_fraction.to_numpy()>=.05)
   v=d.V3_SPP1_macro_share.to_numpy();targetmask=detected[target]&np.isfinite(y)
   reason='';delta=np.nan;n=0;nt=0;nh=0;nl=0
   if mask.sum()<100 or targetmask.sum()<20:reason='insufficient_source_or_marker_qualified_targets'
   else:
    hi,lo=np.quantile(v[mask],[.8,.2])
    ids=np.flatnonzero(targetmask&(y>=np.quantile(y[targetmask],.8)));nt=len(ids)
    if nt<20 or hi<=lo:reason='insufficient_top_marker_targets_or_source_variation'
    else:
     distances=dijkstra(G,directed=False,indices=ids,min_only=True)
     highs=mask&(v>=hi)&np.isfinite(distances);lows=mask&(v<=lo)&np.isfinite(distances);nh=int(highs.sum());nl=int(lows.sum());n=nh+nl
     if min(nh,nl)<20:reason='insufficient_reachable_sources'
     else:delta=float(np.median(distances[highs])-np.median(distances[lows]))
   geometry.append(dict(sample=r.sample,target=target,target_spots=nt,high_sources=nh,low_sources=nl,delta_distance_um=delta,reason=reason,definition='expression rich region, not segmented cells'))
   add('expression_region_distance','SPP1_high_minus_low_distance_to_'+target,'fixed_20percent',delta,n,reason,'difference_um')
  print('SAMPLE',i+1,r.sample,flush=True)
  del X,counts,log,D,H,G,W;gc.collect()
 tab(pd.DataFrame(coverage),'MARKER_ELIGIBILITY.tsv');tab(pd.DataFrame(geometry),'DISTANCE_ELIGIBILITY.tsv')
 gate('patient_statistics');su=summarize(effects,'MARKER_STRUCTURE')
 gate('shared_component_decomposition')
 f=Path('/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/tymp_vascular_niche_v1/tables/spot_mechanism_scores_v1.tsv.gz')
 before=sha(f);sp=pd.read_csv(f,sep='\t');rows=[];checks=[]
 for sid,g in sp.groupby('sample'):
  gate('composite_'+sid)
  def z(x):return (x-x.mean())/x.std(ddof=0)
  g=g.copy();raw=g.total_immune_q05-g.T_cell_q05-g.SPP1_TAM_q05
  assert raw.min()>-1e-5
  g['T_balance_no_SPP1']=z(np.log1p(g.T_cell_q05))-z(np.log1p(np.maximum(raw,0)))
  g['T_absolute_log']=np.log1p(g.T_cell_q05)
  reg=regs[regs['sample']==sid];assert len(reg)==1
  rr=reg.iloc[0]
  for grid in ['4x4','5x5']:
   b=g.groupby('block_'+grid)[['spp1_vascular_niche','T_cell_balance','T_balance_no_SPP1','T_absolute_log','spp1_percentile','vascular_compartment']].median()
   for target in ['T_cell_balance','T_balance_no_SPP1','T_absolute_log']:
    rho=float(stats.spearmanr(b.spp1_vascular_niche,b[target]).statistic)
    rows.append(dict(sample=sid,dataset=rr.dataset,patient=rr.patient,family='shared_composite_recheck',endpoint='niche__'+target,model=grid,effect_type='rho',effect=rho,eligible_spots=len(b),reason=''))
   rho=float(stats.spearmanr(b.spp1_percentile,b.vascular_compartment).statistic)
   rows.append(dict(sample=sid,dataset=rr.dataset,patient=rr.patient,family='direct_vascular_recheck',endpoint='SPP1__vascular_compartment',model=grid,effect_type='rho',effect=rho,eligible_spots=len(b),reason=''))
  checks.append(dict(sample=sid,negative_residual_min=float(raw.min()),original_rows=len(g)))
 assert sha(f)==before
 cs=summarize(rows,'COMPOSITE_RECHECK');tab(pd.DataFrame(checks),'COMPOSITE_DENOMINATOR_QC.tsv')
 for r in inputs:
  f=Path(r['path']);assert f.stat().st_size==r['bytes'] and f.stat().st_mtime_ns==r['mtime_ns']
 dump(P/'INPUT_READONLY_MANIFEST.json',inputs)
 dump(P/'FINAL_VALIDATION.json',dict(status='PASS',at=now(),sections=38,marker_eligibility_rows=len(coverage),structure_section_endpoints=len(effects),composite_section_endpoints=len(rows),old_source_metadata_unchanged=len(inputs),old_composite_sha256_unchanged=True,new_figures=0,independent_cell_identity_claim=False))
 state('COMPLETE_VERIFIED');(P/'EXIT_STATUS').write_text('0\n');(P/'RUN_FINISHED').write_text(now())
 print(cs.to_string(index=False),flush=True)
try:main()
except BaseException:
 (P/'logs/ERROR.txt').write_text(traceback.format_exc());state('FAILED',error=traceback.format_exc());(P/'EXIT_STATUS').write_text('1\n');(P/'RUN_FINISHED').write_text(now());raise

