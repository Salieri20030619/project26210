from pathlib import Path
import json,datetime,hashlib,ast,sys,subprocess,itertools,math
import numpy as np,pandas as pd,scipy
from scipy import sparse,stats
P=Path(__file__).resolve().parents[1];R=P.parent
M=R/'phase200_union_C4_C10_C15_reassessment_20260918_v1'
S=R/'phase157_ordered_bridge_reassessment_20260915_v1'
OLD=R/'phase234_part4_replan_execution_20260921_v1'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def dump(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
def tab(d,f):pd.DataFrame(d).to_csv(P/'tables'/f,sep='\t',index=False)
def read(f,**kw):return pd.read_csv(f,sep='\t',**kw)
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def gate(stage):
 q=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
 with (P/'logs/memory_gate.log').open('a') as f:f.write(stage+'\n'+q.stdout+q.stderr)
 if q.returncode:raise RuntimeError('MEMORY_GATE')
tree=ast.parse((OLD/'scripts/common.py').read_text());exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['fit','infer']],type_ignores=[]),'inherited','exec'),globals())
def bh(d,keys):
 for _,g in d.groupby(keys,dropna=False):
  for pc,qc in [('P','FDR'),('P_wilcoxon','FDR_wilcoxon')]:
   if pc in d:d.loc[g.index,qc]=np.where(g[pc].notna(),stats.false_discovery_control(g[pc].fillna(1).to_numpy()),np.nan)
  d.loc[g.index,'family_size']=len(g)
 return d
def cross(x,y,study,extra=None,seed=280,B=99999):
 x=np.asarray(x,float);y=np.asarray(y,float);study=np.asarray(study)
 if extra is None:extra=np.empty((len(x),0))
 extra=np.asarray(extra,float);extra=extra[:,None] if extra.ndim==1 else extra
 keep=np.isfinite(x)&np.isfinite(y)&np.isfinite(extra).all(1);x=x[keep];y=y[keep];study=study[keep];extra=extra[keep]
 n=len(x);o=dict(n=n,rho=np.nan,P=np.nan,method='DESCRIPTIVE_N_LT10',design_rank=0,residual_df=0,permutations=0,movable_patients=0)
 if n<3 or min(len(np.unique(x)),len(np.unique(y)))<3:return o
 D=np.column_stack([np.ones(n)]+[(study==s).astype(float) for s in sorted(set(study))[1:]]+[stats.rankdata(extra[:,j]) for j in range(extra.shape[1])])
 rk=np.linalg.matrix_rank(D);o.update(design_rank=rk,residual_df=n-rk-1)
 if rk<D.shape[1] or n-rk<3:o['method']='RANK_OR_RESIDUAL_DF';return o
 Q=np.linalg.qr(D,mode='reduced')[0];a=stats.rankdata(x);b=stats.rankdata(y);a=a-Q@(Q.T@a);b=b-Q@(Q.T@b)
 den=np.linalg.norm(a)*np.linalg.norm(b)
 if den<1e-10:o['method']='RESIDUAL_CONSTANT';return o
 obs=float(a@b/den);o['rho']=obs
 if n<10:return o
 gs=[np.flatnonzero(study==s) for s in sorted(set(study))];possible=math.prod(math.factorial(len(g)) for g in gs)
 o.update(possible_permutations=str(possible),movable_patients=sum(len(g) for g in gs if len(g)>1))
 exact=possible<=B
 if exact:
  perms=np.tile(np.arange(n),(possible,1))
  for i,choices in enumerate(itertools.product(*[list(itertools.permutations(g.tolist())) for g in gs])):
   for g,ch in zip(gs,choices):perms[i,g]=ch
 else:
  rng=np.random.default_rng(seed);perms=np.tile(np.arange(n),(B,1))
  for g in gs:
   if len(g)>1:perms[:,g]=g[np.argsort(rng.random((B,len(g))),axis=1)]
 V=b[perms];V-=V@Q@Q.T
 den=np.linalg.norm(V,axis=1)*np.linalg.norm(a);null=np.divide(V@a,den,out=np.full(len(V),np.nan),where=den>1e-10)
 valid=np.isfinite(null);hits=np.count_nonzero(np.abs(null[valid])>=abs(obs)-1e-12)
 o.update(P=hits/valid.sum() if exact else (hits+1)/(valid.sum()+1),permutations=int(valid.sum()),method=('exact' if exact else 'MC_plus1')+'_study_blocked'+('_residual_model_approximate' if extra.shape[1] else ''))
 return o
def main():
 gate('start')
 assert sha(P/'PROTOCOL_BEFORE_TESTING.md')==json.loads((P/'FREEZE.json').read_text())['protocol_sha256']
 paths=[M/'inputs/MACRO_METADATA.tsv.gz',M/'inputs/MACRO_COUNTS.npz',M/'inputs/GENES.tsv',M/'core/contrasts/tumor_C4.npz',S/'02_inputs/GROUP_INDEX.tsv',S/'02_inputs/GROUP_PSEUDOBULK.npz',R/'phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv',R/'phase213_immune_restart_execution_20260919_v1/I1_patient_inputs/LIBRARY_METRICS_ALL_VARIANTS.tsv',R/'phase278_metabolic_state_output_receiver_20260928_v1/PROGRAMS_FROZEN.json']
 manifest=[dict(path=str(f),sha256=sha(f)) for f in paths]
 cfg=json.loads(paths[-1].read_text());programs={k:[g for g in v if g not in {'TYMP','SPP1'}] for k,v in cfg['programs'].items()};dump(P/'resources/PROGRAMS.json',programs)
 G=read(paths[4]);Z=np.load(paths[5],allow_pickle=True);gg=list(Z['genes']);gi={g:i for i,g in enumerate(gg)};C=Z['counts'];assert len(G)==len(C)
 coverage=[dict(program=k,requested=len(v),present=sum(g in gi for g in v),missing=';'.join(g for g in v if g not in gi)) for k,v in programs.items()];tab(coverage,'SOURCE_PROGRAM_COVERAGE.tsv')
 assert all(r['present']>=.8*r['requested'] for r in coverage)
 programs={k:[g for g in v if g in gi] for k,v in programs.items()};dump(P/'resources/PROGRAMS_MEASURED.json',programs)
 b=read(paths[6]);assert len(b)==16 and b.patient.is_unique
 ag=[];summary=[]
 for rr in b.itertuples():
  rec=dict(patient=rr.patient,study=rr.study,library=rr.sample_id,C4_fraction=rr.C4_fraction)
  for name,roles in [('tumor',['Malignant']),('macro',['C4','OtherMacrophage']),('C4',['C4'])]:
   rows=G[(G.patient==rr.patient)&(G.sample_id==rr.sample_id)&G.role.isin(roles)]
   assert len(rows) and rows.source_tissue.eq('Tumor').all()
   v=C[rows.group_index.to_numpy(int)].sum(0);nc=int(rows.cells.sum());tot=float(v.sum());log=np.log2(v/tot*1e6+1)
   rec[name+'_cells']=nc
   for gene in ['TYMP','SPP1','LGALS3','CTSB','CTSD','APOE']:rec[name+'_'+gene]=float(log[gi[gene]]) if name!='C4' or nc>=20 else np.nan
   for namep,gs in programs.items():rec[name+'_'+namep]=float(log[[gi[g] for g in gs]].mean()) if name!='C4' or nc>=20 else np.nan
   ag.append(dict(patient=rr.patient,study=rr.study,library=rr.sample_id,compartment=name,cells=nc,total_counts=tot,TYMP_counts=float(v[gi['TYMP']]),SPP1_counts=float(v[gi['SPP1']]),TYMP_log2CPM=float(log[gi['TYMP']])))
   if name=='tumor':assert tot==rr.malignant_total_counts and nc==rr.malignant_cells and abs(log[gi['TYMP']]-rr.TYMP_pseudobulk)<1e-9
   if name=='macro':assert tot==rr.macrophage_total_counts and nc==rr.macrophage_cells and abs(log[gi['SPP1']]-rr.SPP1_pseudobulk)<1e-9
  summary.append(rec)
 b=pd.DataFrame(summary);tab(b,'PAIRED_SOURCE_PATIENT_INPUTS.tsv');ag=pd.DataFrame(ag);tab(ag,'SOURCE_COUNTS_AND_PROVENANCE.tsv')
 # Paired source abundance difference; not tissue source fraction.
 piv=ag.pivot(index='patient',columns='compartment',values='TYMP_log2CPM')
 desc=pd.DataFrame(dict(patient=piv.index,macro_minus_tumor=piv['macro']-piv['tumor']));tab(desc,'SOURCE_EXPRESSION_PAIRED_EFFECTS.tsv')
 tab([dict(endpoint='macro_minus_tumor_TYMP_log2CPM',**infer(desc.macro_minus_tumor,280))],'SOURCE_EXPRESSION_PAIRED_SUMMARY.tsv')
 zz=np.load(paths[3],allow_pickle=True);j=list(zz['genes']).index('TYMP');e=zz['target_log2CPM'][:,j]-zz['reference_log2CPM'][:,j]
 tab(pd.DataFrame(dict(patient=zz['patients'],study=zz['study'],C4=zz['target_log2CPM'][:,j],OtherMacro=zz['reference_log2CPM'][:,j],effect=e)),'C4_TYMP_PAIRED_EFFECTS.tsv')
 tab([dict(endpoint='C4_minus_other_TYMP',**infer(e,281))],'C4_TYMP_PAIRED_SUMMARY.tsv')
 endpoints=['C4_fraction','macro_SPP1','macro_LGALS3','macro_lysosome_disjoint','macro_polyamine_MTA_context','C4_SPP1']
 out=[];leave=[]
 for source in ['tumor_TYMP','macro_TYMP']:
  other='macro_TYMP' if source=='tumor_TYMP' else 'tumor_TYMP'
  for y in endpoints:
   for model in ['study_adjusted','study_other_source_adjusted']:
    ext=b[other].to_numpy() if model.endswith('other_source_adjusted') else None
    o=cross(b[source],b[y],b.study,ext,seed=2621280+len(out));out.append(dict(source=source,endpoint=y,model=model,**o))
    for st in sorted(b.study.unique()):
     z=b[b.study!=st]
     le=cross(z[source],z[y],z.study,z[other].to_numpy() if ext is not None else None,B=0)
     leave.append(dict(source=source,endpoint=y,model=model,omitted_study=st,n=le['n'],rho=le['rho']))
 tab(bh(pd.DataFrame(out),['source']),'PAIRED_SOURCE_ASSOCIATIONS.tsv');tab(leave,'PAIRED_SOURCE_LEAVE_STUDY.tsv')
 print('PAIRED_COMPLETE',flush=True)
 # Patient-within-macrophage TYMP links.
 m=read(paths[0],index_col=0,dtype={'candidate':str});X=sparse.load_npz(paths[1]).tocsr();genes=read(paths[2]).gene.tolist();lookup={g:i for i,g in enumerate(genes)}
 total=np.asarray(X.sum(1)).ravel();assert len(m)==22444 and m.candidate.eq('4').sum()==2674 and np.allclose(total,m.library_counts)
 selected=['TYMP','SPP1','LGALS3','CTSB','CTSD'];union=sorted(set(selected+sum(list(programs.values()),[])));jj=[lookup[g] for g in union]
 raw=X[:,jj].toarray().astype(float);log=np.log1p(raw/total[:,None]*1e4);idx={g:i for i,g in enumerate(union)};effects=[];caches=[]
 for patient,z in m[m.source_tissue.str.lower().eq('tumor')].groupby('curated_patient_id'):
  ids=m.index.get_indexer(z.index);rv=raw[ids];lv=log[ids];vals={g:lv[:,idx[g]] for g in selected};dets={g:rv[:,idx[g]]>0 for g in selected}
  for pn,gs in programs.items():
   vv=lv[:,[idx[g] for g in gs]];sd=vv.std(0);vals[pn]=((vv-vv.mean(0))/np.where(sd>0,sd,1)).mean(1);dets[pn]=(rv[:,[idx[g] for g in gs]]>0).sum(1)>=2
  dep=np.log1p(total[ids]);ng=z.n_genes_raw.to_numpy()
  caches.append(pd.DataFrame(dict(cell_id=z.index,patient=patient,study=z.study.iloc[0],candidate=z.candidate.to_numpy(),log_depth=dep,n_genes=ng,**vals,**{'det_'+k:v for k,v in dets.items()})))
  for group,mask,minn in [('C4',z.candidate.eq('4').to_numpy(),50),('WholeMacro',np.ones(len(z),bool),100)]:
   cov=np.column_stack([dep,ng])
   if group=='WholeMacro':cov=np.column_stack([cov,pd.get_dummies(z.candidate,dtype=float,drop_first=True).to_numpy()])
   for y in ['SPP1','LGALS3','CTSB','CTSD']+list(programs):
    nx=int((dets['TYMP']&mask).sum());ny=int((dets[y]&mask).sum());eligible=mask.sum()>=minn and nx>=10 and ny>=(20 if y in programs else 10)
    for model in ['unadjusted','adjusted']:
     if eligible:eff,n,reason=fit(vals['TYMP'],vals[y],None if model=='unadjusted' else cov,mask,min_n=minn)
     else:eff,n,reason=np.nan,int(mask.sum()),'CELL_DETECTION_GATE'
     effects.append(dict(patient=patient,study=z.study.iloc[0],group=group,endpoint=y,model=model,effect=eff,cells=int(mask.sum()),TYMP_detected=nx,outcome_detected=ny,reason=reason))
 ef=pd.DataFrame(effects);tab(ef,'WITHIN_MACRO_PATIENT_EFFECTS.tsv');pd.concat(caches,ignore_index=True).to_csv(P/'inputs/CELL_CACHE.tsv.gz',sep='\t',index=False)
 su=[]
 for key,z in ef.groupby(['group','endpoint','model']):su.append(dict(zip(['group','endpoint','model'],key),**infer(z.effect,262280+len(su))))
 tab(bh(pd.DataFrame(su),['group']),'WITHIN_MACRO_SUMMARY.tsv')
 study=[];leave=[]
 for key,z in ef.groupby(['group','endpoint','model']):
  for st in sorted(z.study.unique()):
   v=z.loc[z.study==st,'effect'].dropna();study.append(dict(zip(['group','endpoint','model'],key),study=st,n=len(v),positive=int((v>0).sum()),negative=int((v<0).sum()),median=v.median()))
   v=z.loc[z.study!=st,'effect'].dropna();leave.append(dict(zip(['group','endpoint','model'],key),omitted_study=st,n=len(v),positive=int((v>0).sum()),negative=int((v<0).sum()),median=v.median()))
 tab(study,'WITHIN_MACRO_STUDY_DIRECTIONS.tsv');tab(leave,'WITHIN_MACRO_LEAVE_STUDY.tsv')
 print('WITHIN_MACRO_COMPLETE',flush=True)
 # Existing frozen receiver eligibility, same-library source matching.
 libs=read(paths[7]);libs=libs[(libs.variant=='primary')&libs.eligible&libs.endpoint.isin(['CD8_cytotoxicity','CD8_dysfunction'])]
 joined=libs.merge(b,left_on=['patient','study','library'],right_on=['patient','study','library'],validate='many_to_one')
 cols=['patient','study','library','endpoint','outcome','receiver_cells','tumor_TYMP','macro_TYMP'];tab(joined[cols],'TYMP_CD8_MATCHED_LIBRARIES.tsv')
 pa=joined.groupby(['patient','study','endpoint']).agg(y=('outcome','mean'),tumor_TYMP=('tumor_TYMP','mean'),macro_TYMP=('macro_TYMP','mean'),libraries=('library',lambda v:';'.join(sorted(set(v))))).reset_index();tab(pa,'TYMP_CD8_MATCHED_PATIENTS.tsv');ou=[]
 for endpoint,z in pa.groupby('endpoint'):
  for source in ['tumor_TYMP','macro_TYMP']:
   other='macro_TYMP' if source=='tumor_TYMP' else 'tumor_TYMP'
   for model in ['study_adjusted','study_other_source_adjusted']:
    res=cross(z[source],z.y,z.study,z[other].to_numpy() if model.endswith('other_source_adjusted') else None,seed=2623280+len(ou))
    ou.append(dict(source=source,endpoint=endpoint,model=model,**res))
 tab(bh(pd.DataFrame(ou),['source']),'TYMP_CD8_FUNCTION_SUMMARY.tsv')
 for rec in manifest:assert sha(rec['path'])==rec['sha256']
 dump(P/'qc/SOURCE_INPUTS.json',manifest);dump(P/'qc/SOURCE_COMPLETED.json',dict(at=now(),within_patient_rows=len(ef),matched_patients=pa.patient.nunique(),python=sys.version,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__))
 print('SOURCE_FUNCTION_COMPLETE',flush=True)
if __name__=='__main__':main()
