from pathlib import Path
import json,datetime,hashlib,ast,traceback,gc,subprocess,sys
import numpy as np,pandas as pd
from scipy import sparse,stats
P=Path(__file__).resolve().parents[1];R=P.parent;SRC=R/'phase200_union_C4_C10_C15_reassessment_20260918_v1';OLD=R/'phase234_part4_replan_execution_20260921_v1'
sys.path.insert(0,str(R/'phase213_immune_restart_execution_20260919_v1/scripts'))
from stats_shared import main_test,bootstrap_effect
def now():return datetime.datetime.now().astimezone().isoformat()
def dump(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
def tab(d,f):pd.DataFrame(d).to_csv(P/'tables'/f,sep='\t',index=False)
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def status(s,**kw):dump(P/'STATUS.json',dict(status=s,at=now(),**kw))
def gate(stage):
 q=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
 with open(P/'logs/memory_gate.log','a') as f:f.write(stage+'\n'+q.stdout+q.stderr)
 if q.returncode:raise RuntimeError('MEMORY_GATE')
tree=ast.parse((OLD/'scripts/common.py').read_text());exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['fit','infer']],type_ignores=[]),'inherited','exec'),globals())
def main():
 gate('start');seal=json.loads((P/'FREEZE.json').read_text());assert sha(P/'PROTOCOL_BEFORE_TESTING.md')==seal['protocol_sha256']
 cfg=json.loads((P/'PROGRAMS_FROZEN.json').read_text());programs=cfg['programs'];outputs=cfg['outputs'];assert all(not set(gs)&set(outputs) for gs in programs.values())
 paths=[SRC/'inputs/MACRO_METADATA.tsv.gz',SRC/'inputs/MACRO_COUNTS.npz',SRC/'inputs/GENES.tsv',SRC/'core/contrasts/tumor_C4_genes.tsv.gz',R/'phase213_immune_restart_execution_20260919_v1/I1_patient_inputs/LIBRARY_METRICS_ALL_VARIANTS.tsv']
 manifest=[dict(path=str(f),sha256=sha(f)) for f in paths]
 m=pd.read_csv(paths[0],sep='\t',index_col=0,dtype={'candidate':str});X=sparse.load_npz(paths[1]).tocsr();genes=pd.read_csv(paths[2],sep='\t').gene.tolist();lookup={g:i for i,g in enumerate(genes)}
 assert len(m)==22444 and int((m.candidate=='4').sum())==2674 and X.shape==(len(m),len(genes))
 total=np.asarray(X.sum(1)).ravel();assert np.allclose(total,m.library_counts,atol=1e-6,rtol=0)
 used={name:[g for g in gs if g in lookup] for name,gs in programs.items()}
 assert all(len(used[n])>=.8*len(gs) for n,gs in programs.items()) and all(g in lookup for g in outputs)
 union=sorted(set(outputs+sum(list(used.values()),[])));cols=[lookup[g] for g in union];raw=X[:,cols].toarray().astype(float);log=np.log1p(raw/total[:,None]*10000);idx={g:i for i,g in enumerate(union)}
 old=pd.read_csv(paths[3],sep='\t');source=old[old.gene.isin(outputs)].copy();assert len(source)==15
 source['source_nominal_enriched']=(source.P<.05)&(source.median_delta>0);source['P_FDR_scope']='original Phase200 tumor_C4 all-gene family; not recomputed'
 tab(source,'INHERITED_OUTPUT_ENRICHMENT.tsv');tab([dict(program=n,requested=len(programs[n]),present=len(gs),genes=';'.join(gs)) for n,gs in used.items()],'PROGRAM_GENE_COVERAGE.tsv')
 tumor=m.source_tissue.str.lower().eq('tumor').to_numpy();effects=[];cellcache=[];patientmeta=[]
 for k,(pat,gg) in enumerate(m[tumor].groupby('curated_patient_id',sort=True)):
  ids=m.index.get_indexer(gg.index);lv=log[ids];rv=raw[ids];dep=np.log1p(total[ids]);ng=gg.n_genes_raw.to_numpy()
  scores={};det={}
  for name,gs in used.items():
   jj=[idx[g] for g in gs];v=lv[:,jj];sd=v.std(0);scores[name]=((v-v.mean(0))/np.where(sd>0,sd,1)).mean(1);det[name]=(rv[:,jj]>0).sum(1)>=2
  cache=pd.DataFrame(dict(cell_id=gg.index,patient=pat,study=gg.study.iloc[0],candidate=gg.candidate.to_numpy(),log_depth=dep,n_genes=ng,**scores,**{g:lv[:,idx[g]] for g in outputs},**{'detected_'+n:v for n,v in det.items()}))
  cellcache.append(cache)
  for group,mask in [('C4',gg.candidate.eq('4').to_numpy()),('WholeMacro',np.ones(len(gg),bool))]:
   nc=int(mask.sum());cov=np.column_stack([dep,ng])
   if group=='WholeMacro':
    dummy=pd.get_dummies(gg.candidate,dtype=float,drop_first=True).to_numpy();cov=np.column_stack([cov,dummy])
   for program,sc in scores.items():
    for gene in outputs:
     yy=lv[:,idx[gene]];nout=int(((rv[:,idx[gene]]>0)&mask).sum());nprog=int((det[program]&mask).sum())
     eligible=nc>=(50 if group=='C4' else 100) and nout>=10 and nprog>=20
     for model in ['unadjusted','depth_complexity_adjusted' if group=='C4' else 'depth_complexity_cluster_adjusted']:
      if eligible:e,n,why=fit(sc,yy,cov=None if model=='unadjusted' else cov,mask=mask,min_n=50 if group=='C4' else 100)
      else:e=np.nan;n=nc;why='CELL_OR_DETECTION_GATE'
      effects.append(dict(patient=pat,study=gg.study.iloc[0],group=group,program=program,gene=gene,model=model,effect=e,cells=nc,output_detected_cells=nout,program_detected_cells=nprog,reason=why))
   patientmeta.append(dict(patient=pat,study=gg.study.iloc[0],group=group,cells=nc,libraries=';'.join(sorted(set(gg.sample_id[mask])))))
  print('PATIENT_COUPLING',k+1,pat,flush=True)
 pd.concat(cellcache,ignore_index=True).to_csv(P/'inputs/CELL_SCORE_CACHE.tsv.gz',sep='\t',index=False);tab(patientmeta,'SOURCE_PATIENT_QUALIFICATION.tsv')
 ef=pd.DataFrame(effects);tab(ef,'WITHIN_PATIENT_EFFECTS.tsv');keys=['group','program','gene','model'];rows=[]
 for i,(key,g) in enumerate(ef.groupby(keys,sort=True)):rows.append(dict(zip(keys,key),**infer(g.effect,262102780+i)))
 su=pd.DataFrame(rows)
 for _,g in su.groupby(['group','program']):
  assert len(g)==30
  for pc,qc in [('P','FDR'),('P_wilcoxon','FDR_wilcoxon')]:su.loc[g.index,qc]=np.where(g[pc].notna(),stats.false_discovery_control(g[pc].fillna(1).to_numpy()),np.nan)
  su.loc[g.index,'family_size']=30
 su['evidence']=np.where(su.P.isna(),'DESCRIPTIVE_OR_UNESTIMABLE',np.where(su.P<.05,np.where(su.FDR<.05,'FDR_SIGNIFICANT','NOMINAL_SIGNIFICANT'),'NOT_SUPPORTED'))
 tab(su,'WITHIN_PATIENT_SUMMARY.tsv')
 studyrows=[]
 for key,g in ef.groupby(keys+['study']):
  v=g.effect.dropna().to_numpy();studyrows.append(dict(zip(keys+['study'],key),n=len(v),positive=int((v>0).sum()),negative=int((v<0).sum()),median=float(np.median(v)) if len(v) else np.nan))
 tab(studyrows,'STUDY_DIRECTION_SUMMARY.tsv');leave=[]
 nominated=su[(su.group=='C4')&(su.model=='depth_complexity_adjusted')&(su.P<.05)&(su['median']>0)]
 for r in nominated.itertuples():
  g=ef[(ef.group==r.group)&(ef.program==r.program)&(ef.gene==r.gene)&(ef.model==r.model)]
  for st in sorted(g.study.unique()):
   v=g.loc[g.study!=st,'effect'].dropna().to_numpy();leave.append(dict(program=r.program,gene=r.gene,omitted_study=st,n=len(v),positive=int((v>0).sum()),negative=int((v<0).sum()),median=float(np.median(v)) if len(v) else np.nan))
 tab(leave,'NOMINAL_C4_LEAVE_ONE_STUDY.tsv')
 # Matched libraries: sender summaries in the exact same frozen library.
 libs=pd.read_csv(paths[4],sep='\t');selected=libs[(libs.variant=='primary')&libs.endpoint.isin(['CD8_cytotoxicity','CD8_dysfunction'])&libs.C4_intrinsic_eligible]
 rr=[]
 for r in selected.itertuples():
  mask=tumor&m.curated_patient_id.eq(r.patient).to_numpy()&m.sample_id.eq(r.library).to_numpy()&m.candidate.eq('4').to_numpy()
  assert int(mask.sum())==r.n_C4 and mask.sum()>=20
  den=total[mask].sum()
  for program,gs in used.items():
   counts=raw[np.ix_(mask,[idx[g] for g in gs])].sum(0);value=float(np.log2(counts/den*1e6+1).mean())
   rr.append(dict(patient=r.patient,study=r.study,library=r.library,program=program,endpoint=r.endpoint,sender_C4_cells=int(mask.sum()),receiver_cells=r.receiver_cells,x=value,y=r.outcome))
 lib=pd.DataFrame(rr);tab(lib,'MATCHED_CD8_LIBRARY_INPUTS.tsv')
 pats=lib.groupby(['patient','study','program','endpoint']).agg(x=('x','mean'),y=('y','mean'),libraries=('library',lambda v:';'.join(sorted(v))),n_libraries=('library','size')).reset_index();tab(pats,'MATCHED_CD8_PATIENT_INPUTS.tsv')
 summaries=[]
 for i,((program,endpoint),g) in enumerate(pats.groupby(['program','endpoint'],sort=True)):
  for model in ['unadjusted','study_adjusted']:
   st=g.study.to_numpy() if model=='study_adjusted' else np.repeat('one_block',len(g))
   result=dict(program=program,endpoint=endpoint,model=model,n=len(g),rho=np.nan,P=np.nan,CI_low=np.nan,CI_high=np.nan,method='DESCRIPTIVE_N_LT10',permutations=0,possible_permutations='',residual_df=np.nan)
   if len(g)>=10:
    out,dist=main_test(g.x.to_numpy(),g.y.to_numpy(),st,seed=2621278+i,B=99999);result.update(rho=out['rho'],P=out['p'],method=out['method'],permutations=out['permutations'],possible_permutations=out['possible_permutations'],residual_df=out['residual_df'])
    low,high,valid,boot=bootstrap_effect(g.x,g.y,st,B=5000,seed=2622278+i);result.update(CI_low=low,CI_high=high,valid_bootstraps=valid)
    np.savez_compressed(P/'inputs'/f'NULL_{program}_{endpoint}_{model}.npz',permutation=dist,bootstrap=boot)
   elif len(g)>=3:
    from stats_shared import effect
    result['rho']=effect(g.x,g.y,st)
   summaries.append(result)
 cs=pd.DataFrame(summaries);assert len(cs)==8;cs['FDR']=np.where(cs.P.notna(),stats.false_discovery_control(cs.P.fillna(1)),np.nan);cs['family_size']=8;tab(cs,'MATCHED_CD8_FUNCTION_SUMMARY.tsv')
 for rec in manifest:assert sha(rec['path'])==rec['sha256']
 dump(P/'qc/SOURCE_READONLY_MANIFEST.json',manifest)
 dump(P/'qc/SOURCE_COMPUTATION_COMPLETED.json',dict(at=now(),within_rows=len(ef),within_endpoints=len(su),matched_endpoints=len(cs),source_cells=len(m),tumor_cells=int(tumor.sum()),figures=0))
 status('SOURCE_COMPLETE_RECEIVER_PENDING');print('SOURCE_COMPLETE',len(su),len(cs),flush=True)
if __name__=='__main__':
 try:main()
 except Exception:status('FAILED',traceback=traceback.format_exc());raise
