from pathlib import Path
import json,hashlib,datetime,subprocess,math,itertools,sys
import numpy as np,pandas as pd,scipy
from scipy import sparse,stats
P=Path(__file__).resolve().parents[1];R=P.parent
def read(p):return pd.read_csv(p,sep='\t',float_precision='round_trip')
def tab(x,n):pd.DataFrame(x).to_csv(P/'tables'/n,sep='\t',index=False)
def dump(n,x):(P/n).write_text(json.dumps(x,indent=2,ensure_ascii=False,default=lambda z:z.item() if hasattr(z,'item') else str(z)))
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def gate(stage):
 q=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
 with (P/'logs/memory_gate.log').open('a') as f:f.write(stage+'\n'+q.stdout+q.stderr)
 if q.returncode:raise RuntimeError('MEMORY_GATE')
def cross(d,x,y,model,infer=True):
 cov=[] if model=='study_only' else ['perivascular_fraction']
 if model.endswith('receiver_depth'):cov+=['log_receiver_depth']
 cols=[x,y]+cov
 z=d.replace([np.inf,-np.inf],np.nan).dropna(subset=cols)
 n=len(z);o=dict(n=n,studies=z.study.nunique(),rho=np.nan,P=np.nan,method='NOT_ESTIMABLE',residual_df=0,design_rank=0,permutations=0,seed=292)
 if n<3:return o
 D=np.column_stack([np.ones(n)]+[(z.study.to_numpy()==s).astype(float) for s in sorted(z.study.unique())[1:]]+[stats.rankdata(z[c]) for c in cov])
 rank=np.linalg.matrix_rank(D);df=n-rank-1;o.update(design_rank=int(rank),residual_df=int(df))
 if rank!=D.shape[1] or n-rank<2:o['method']='DESIGN_RANK_OR_RESIDUAL_DIM';return o
 Q=np.linalg.qr(D,mode='reduced')[0]
 a=stats.rankdata(z[x]);b=stats.rankdata(z[y]);a=a-Q@(Q.T@a);b=b-Q@(Q.T@b)
 den=np.linalg.norm(a)*np.linalg.norm(b)
 if den<1e-10:o['method']='CONSTANT_RESIDUAL';return o
 rho=float(a@b/den);o['rho']=rho
 if n<10 or df<5 or not infer:o['method']='DESCRIPTIVE';return o
 gs=[np.flatnonzero(z.study.to_numpy()==s) for s in sorted(z.study.unique())]
 possible=math.prod(math.factorial(len(g)) for g in gs);exact=possible<=99999
 N=possible if exact else 99999
 perms=np.tile(np.arange(n),(N,1))
 if exact:
  for i,c in enumerate(itertools.product(*[list(itertools.permutations(g.tolist())) for g in gs])):
   for g,v in zip(gs,c):perms[i,g]=v
 else:
  rng=np.random.default_rng(292)
  for g in gs:
   if len(g)>1:perms[:,g]=g[np.argsort(rng.random((N,len(g))),axis=1)]
 v=b[perms];v=v-v@Q@Q.T;den=np.linalg.norm(v,axis=1)*np.linalg.norm(a)
 null=np.divide(v@a,den,out=np.full(N,np.nan),where=den>1e-10)
 ok=np.isfinite(null);hits=int(np.count_nonzero(np.abs(null[ok])>=abs(rho)-1e-12));valid=int(ok.sum())
 o.update(P=hits/valid if exact else (hits+1)/(valid+1),method=('exact' if exact else 'MC_plus1')+'_study_blocked'+('_residual_approximate' if cov else ''),permutations=valid,possible_permutations=str(possible))
 return o
def main():
 gate('load_inputs')
 f=json.loads((P/'FREEZE.json').read_text())
 assert sha(P/'PLAN_AND_PROTOCOL.md')==f['protocol_sha256'] and sha(P/'PANELS_FROZEN.json')==f['panels_sha256']
 manifest=json.loads((P/'qc/INPUT_MANIFEST.json').read_text())
 for q in manifest:assert sha(q['path'])==q['sha256']
 panels=json.loads((P/'PANELS_FROZEN.json').read_text())
 old=R/'phase287_stromal_refinement_image_morphology_20260929_v1'
 m=read(R/'phase289_state_source_boundary_process_20260929_v1/tables/STROMAL_STATE_MEMBERSHIP.tsv.gz')
 original=read(old/'inputs/STROMAL_INPUT_METADATA.tsv.gz')
 assert m.cell_id.is_unique and original.cell_id.is_unique
 assert m.cell_id.tolist()==original.cell_id.tolist()
 X=sparse.load_npz(old/'inputs/STROMAL_COUNTS.npz').tocsr()
 genes=np.load(old/'inputs/GENES.npy',allow_pickle=True).tolist();gi={g:i for i,g in enumerate(genes)}
 assert len(gi)==len(genes) and X.shape==(len(m),len(genes))
 assert np.isfinite(X.data).all() and (X.data>=0).all() and (X.data==np.floor(X.data)).all()
 totals=np.asarray(X.sum(axis=1)).ravel()
 assert np.array_equal(totals,m.total_counts_raw.to_numpy())
 pool=m.stromal_pool.astype(str).str.lower().eq('true')
 assert int(pool.sum())==15438
 coverage=[]
 for name,gs in panels.items():
  kept=[g for g in gs if g in gi]
  coverage.append(dict(program=name,requested=len(gs),present=len(kept),eligible=len(kept)>=.8*len(gs),genes=';'.join(kept),missing=';'.join(set(gs)-set(kept))))
 tab(coverage,'PANEL_COVERAGE.tsv')
 union=sorted(set(sum(panels.values(),[]))&set(genes))
 counts=X[:,[gi[g] for g in union]].toarray().astype(float)
 lookup={g:i for i,g in enumerate(union)}
 paired=read(R/'phase289_state_source_boundary_process_20260929_v1/tables/MATCHED_SOURCE_RECEIVER_INPUTS.tsv')
 assert len(paired)==16 and paired.patient.is_unique
 G=read(R/'phase157_ordered_bridge_reassessment_20260915_v1/02_inputs/GROUP_INDEX.tsv')
 Z=np.load(R/'phase157_ordered_bridge_reassessment_20260915_v1/02_inputs/GROUP_PSEUDOBULK.npz',allow_pickle=True)
 gg=Z['genes'].tolist();CC=Z['counts'];assert np.array_equal(G.group_index,np.arange(len(G)))
 rows=[];gene_rows=[];subrows=[];pb=[];subpb=[]
 for r in paired.itertuples():
  gate('aggregate_'+r.patient)
  use=(m.curated_patient_id.eq(r.patient)&m.sample_id.eq(r.library)&m.source_tissue.eq('Tumor')&pool).to_numpy()
  n=int(use.sum());assert n==r.stromal_pool_cells
  rec=dict(patient=r.patient,study=r.study,library=r.library,stromal_cells=n)
  for label,roles in [('malignant',['Malignant']),('macro',['C4','OtherMacrophage'])]:
   hit=G.patient.eq(r.patient)&G.sample_id.eq(r.library)&G.source_tissue.eq('Tumor')&G.role.isin(roles)
   vec=CC[G.loc[hit,'group_index'].to_numpy(int)].sum(axis=0)
   nc=int(G.loc[hit,'cells'].sum());total=float(vec.sum())
   assert nc==getattr(r,label+'_cells')
   rec[label+'_cells']=nc
   for gene in ['TYMP','SPP1']:
    val=float(np.log1p(vec[gg.index(gene)]/total*1e4)) if total>0 else np.nan
    assert np.isclose(val,getattr(r,label+'_'+gene),equal_nan=True)
    rec[label+'_'+gene]=val
  rec['eligible']=(n>=50 and rec['malignant_cells']>=50 and rec['macro_cells']>=100)
  assert rec['eligible']==bool(r.eligible)
  v=counts[use].sum(0);total=float(totals[use].sum())
  log=np.log1p(v/max(total,1)*1e4)
  rec.update(total_UMI=total,perivascular_fraction=float(m.loc[use,'cell_type_broad_v1'].eq('Pericyte/SMC').mean()) if n else np.nan,log_receiver_depth=float(np.log1p(total/n)) if n else np.nan)
  for name,gs in panels.items():
   present=[lookup[g] for g in gs if g in lookup];det=int(np.count_nonzero(v[present]>0))
   ok=len(present)>=.8*len(gs) and det>=math.ceil(len(present)/2) and total>0
   rec[name]=float(log[present].mean()) if ok else np.nan
   rec[name+'_genes_detected']=det
  rows.append(rec);pb.append(v)
  for j,g in enumerate(union):gene_rows.append(dict(patient=r.patient,study=r.study,library=r.library,gene=g,counts=v[j],positive_cells=int(np.count_nonzero(counts[use,j])),logCP10K=log[j]))
  for celltype in ['Fibroblast','Pericyte/SMC']:
   sel=use&m.cell_type_broad_v1.eq(celltype).to_numpy();sn=int(sel.sum());stotal=float(totals[sel].sum());sv=counts[sel].sum(0)
   slog=np.log1p(sv/max(stotal,1)*1e4)
   sr=dict(patient=r.patient,study=r.study,library=r.library,celltype=celltype,cells=sn,total_UMI=stotal,interpretation='descriptive_not_pure_CAF')
   for name,gs in panels.items():sr[name]=float(slog[[lookup[g] for g in gs if g in lookup]].mean()) if sn else np.nan
   subrows.append(sr);subpb.append(sv)
 d=pd.DataFrame(rows);assert d.eligible.sum()==13
 tab(d,'PATIENT_INPUTS.tsv');tab(gene_rows,'PATIENT_GENE_COUNTS.tsv');tab(subrows,'BROAD_CELLTYPE_DESCRIPTIVE.tsv')
 np.savez_compressed(P/'inputs/PAIRED_STROMAL_PSEUDOBULK.npz',genes=np.array(union),counts=np.array(pb),patients=d.patient.to_numpy(str),total_UMI=d.total_UMI.to_numpy(),subtype_counts=np.array(subpb))
 gate('patient_associations')
 models=['study_only','study_plus_perivascular_fraction','study_plus_perivascular_fraction_plus_receiver_depth']
 sources=['malignant_TYMP','macro_TYMP','macro_SPP1'];out=[];leave=[];stud=[];loo_gene=[]
 for x in sources:
  for model in models:
   for y in panels:
    z=d[d.eligible&d[y].notna()].copy()
    result=cross(z,x,y,model);out.append(dict(source=x,endpoint=y,model=model,family=x+'__'+model,family_size=3,**result))
    for st in sorted(z.study.unique()):
     w=z[z.study.ne(st)];s=cross(w,x,y,model,False);leave.append(dict(source=x,endpoint=y,model=model,omission_type='study',omitted=st,**s))
     w=z[z.study.eq(st)];s=cross(w,x,y,'study_only',False);stud.append(dict(source=x,endpoint=y,study=st,**s))
    for pat in z.patient:
     s=cross(z[z.patient.ne(pat)],x,y,model,False);leave.append(dict(source=x,endpoint=y,model=model,omission_type='patient',omitted=pat,**s))
    # Fixed leave-one-gene diagnostics for every endpoint/source; no P or selection.
    if model=='study_plus_perivascular_fraction':
     gd=pd.DataFrame(gene_rows).pivot(index='patient',columns='gene',values='logCP10K')
     for gene in panels[y]:
      gs=[g for g in panels[y] if g in gd.columns and g!=gene]
      w=z.copy();w['gene_removed_score']=gd.loc[w.patient,gs].mean(axis=1).to_numpy()
      s=cross(w,x,'gene_removed_score',model,False);loo_gene.append(dict(source=x,endpoint=y,removed_gene=gene,**s))
 a=pd.DataFrame(out);a['FDR']=np.nan
 for fam,ix in a.groupby('family').groups.items():
  pv=a.loc[ix,'P'].fillna(1).to_numpy();q=stats.false_discovery_control(pv,method='bh')
  a.loc[ix,'FDR']=np.where(a.loc[ix,'P'].notna(),q,np.nan)
 a['evidence']=np.where(a.P.isna(),'not_estimable',np.where(a.FDR.lt(.05),'FDR_significant',np.where(a.P.lt(.05),'nominal','not_supported')))
 tab(a,'ASSOCIATIONS.tsv');tab(leave,'LEAVE_OUT_DESCRIPTIVE.tsv');tab(pd.DataFrame(stud).drop_duplicates(),'WITHIN_STUDY_DESCRIPTIVE.tsv');tab(loo_gene,'LEAVE_GENE_DESCRIPTIVE.tsv')
 gate('postrun_integrity')
 for q in manifest:assert sha(q['path'])==q['sha256']
 assert sha(P/'PLAN_AND_PROTOCOL.md')==f['protocol_sha256'] and sha(P/'PANELS_FROZEN.json')==f['panels_sha256']
 dump('resources/SOFTWARE.json',dict(executable=sys.executable,python=sys.version,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__))
 dump('qc/COMPUTED.json',dict(status='COMPUTED_PENDING_INDEPENDENT_VALIDATION',at=datetime.datetime.now().astimezone().isoformat(),input_cells=len(m),pool_cells=int(pool.sum()),patients=int(d.eligible.sum()),study_counts=d[d.eligible].study.value_counts().to_dict(),models=len(a),nominal_models=int(a.P.lt(.05).sum()),fdr_models=int(a.FDR.lt(.05).sum())))
 print(a[['source','endpoint','model','n','rho','P','FDR','evidence']].to_string(index=False))
 print('PHASE292_COMPUTED',flush=True)
if __name__=='__main__':main()
