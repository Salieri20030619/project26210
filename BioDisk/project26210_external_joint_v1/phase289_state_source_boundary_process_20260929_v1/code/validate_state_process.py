from common import *
from scipy import stats,sparse,linalg
import math,itertools,h5py
def blocked(x,y,study,extra,seed=289,B=9999):
 x=np.asarray(x,float);y=np.asarray(y,float);study=np.asarray(study)
 e=np.empty((len(x),0)) if extra is None else np.asarray(extra,float)
 good=np.isfinite(x)&np.isfinite(y)&np.isfinite(e).all(1);x=x[good];y=y[good];e=e[good];study=study[good];n=len(x)
 D=np.column_stack([np.ones(n)]+[(study==s).astype(float) for s in sorted(set(study))[1:]]+[stats.rankdata(e[:,j]) for j in range(e.shape[1])])
 xr=stats.rankdata(x);yr=stats.rankdata(y)
 proj=D@np.linalg.pinv(D);a=xr-proj@xr;b=yr-proj@yr
 rho=float(a@b/np.linalg.norm(a)/np.linalg.norm(b))
 if n<10:return rho,np.nan
 gs=[np.flatnonzero(study==s) for s in sorted(set(study))];possible=math.prod(math.factorial(len(g)) for g in gs)
 if possible<=B:
  per=[]
  for choices in itertools.product(*[list(itertools.permutations(g.tolist())) for g in gs]):
   inds=np.arange(n)
   for g,ch in zip(gs,choices):inds[g]=ch
   per.append(inds)
  perms=np.array(per);exact=True
 else:
  rng=np.random.default_rng(seed);perms=np.tile(np.arange(n),(B,1))
  for g in gs:
   if len(g)>1:perms[:,g]=g[np.argsort(rng.random((B,len(g))),axis=1)]
  exact=False
 v=b[perms];v=v-v@proj
 null=v@a/(np.linalg.norm(v,axis=1)*np.linalg.norm(a));keep=np.isfinite(null);hits=int((np.abs(null[keep])>=abs(rho)-1e-12).sum())
 p=hits/keep.sum() if exact else (hits+1)/(keep.sum()+1)
 return rho,p
def main():
 gate('validate_state_process')
 old=R/'phase287_stromal_refinement_image_morphology_20260929_v1'
 m=read(P/'tables/STROMAL_STATE_MEMBERSHIP.tsv.gz');o=read(old/'tables/STROMAL_ANNOTATION_V1.tsv.gz')
 assert m.cell_id.tolist()==o.cell_id.tolist()
 z=np.load(P/'inputs/STROMAL_PANEL_COUNTS.npz');raw=sparse.load_npz(old/'inputs/STROMAL_COUNTS.npz');gg=np.load(old/'inputs/GENES.npy').tolist();genes=z['genes'].tolist();idx={g:i for i,g in enumerate(genes)}
 assert np.array_equal(raw[:,[gg.index(g) for g in genes]].toarray(),z['counts'])
 x=z['counts'];total=z['total_counts'];log=np.log1p(x/total[:,None]*1e4);panels=json.loads((P/'PANELS_FROZEN.json').read_text())
 for prog,gs in panels.items():
  jj=[idx[g] for g in gs if g in idx]
  if len(jj)>=.8*len(gs):
   assert np.allclose(log[:,jj].mean(1),m['program_'+prog],atol=1e-12)
   assert np.array_equal((x[:,jj]>0).sum(1),m['det_'+prog])
 pool=m.cell_type_broad_v1.isin(['Fibroblast','Pericyte/SMC'])&~(m.competing_immune|m.competing_epi|m.competing_endo)&~m.malignant_conservative_v1.astype(str).str.lower().isin(['true','1'])
 act=pool&m.det_CAF_activation.ge(2)&(m.fibro_candidate|m.cluster_primary.isin([9,15]))
 trans=pool&m.det_transition_RCC.ge(2)&m.det_perivascular.ge(1)
 for name,value in [('stromal_pool',pool),('activated_candidate',act),('transition_candidate',trans)]:assert np.array_equal(m[name],value)
 mm=read(R/'phase200_union_C4_C10_C15_reassessment_20260918_v1/inputs/MACRO_METADATA.tsv.gz',dtype={'candidate':str})
 assert len(mm)==22444 and mm.candidate.eq('4').sum()==2674
 assert not (set(mm.cell_id)&set(m.loc[pool,'cell_id']))
 pb=read(P/'tables/STROMAL_STATE_PSEUDOBULK.tsv');zp=np.load(P/'inputs/STROMAL_STATE_PSEUDOBULK_COUNTS.npz')
 masks={'stromal_pool':m.stromal_pool,'activated':m.activated_candidate,'transition':m.transition_candidate,'conservative':m.conservative_fibro_reference}
 for i,r in enumerate(pb.itertuples()):
  keep=m.study.eq(r.study)&m.curated_patient_id.eq(r.patient)&m.sample_id.eq(r.library)&m.source_tissue.eq(r.tissue)&masks[r.population]
  assert keep.sum()==r.cells and np.array_equal(x[keep].sum(0),zp['counts'][i]) and total[keep].sum()==zp['total_counts'][i]
  val=np.log1p(zp['counts'][i]/zp['total_counts'][i]*1e4)
  for prog,gs in panels.items():
   jj=[idx[g] for g in gs if g in idx]
   if len(jj)>=.8*len(gs):assert abs(float(val[jj].mean())-getattr(r,prog))<1e-12
 # Source gene table: frozen source counts and H5 fingerprint checks, plus independent exact group spot checks.
 src=read(P/'tables/SOURCE_GENE_PATIENT_LIBRARY.tsv')
 assert src.loc[(src.source_group=='Frozen_C4')&(src.gene=='TYMP'),'cells'].sum()==2674
 assert src.loc[src.gene=='TYMP','cells'].sum()==258120
 assert src.loc[src.gene=='SPP1','cells'].sum()==258120
 d=read(P/'tables/MATCHED_SOURCE_RECEIVER_INPUTS.tsv');s=read(P/'tables/SOURCE_RECEIVER_ASSOCIATIONS.tsv');use=d.eligible
 assert np.array_equal(use,(d.malignant_cells>=50)&(d.macro_cells>=100)&(d.stromal_pool_cells>=50))
 maxrho=0.;maxp=0.
 for r in s.itertuples():
  extra=None
  if r.model=='cell_count_sensitivity':
   cell='malignant_cells' if r.predictor.startswith('malignant') else 'macro_cells'
   extra=np.column_stack([np.log1p(d.loc[use,cell]),np.log1p(d.loc[use,'stromal_pool_cells'])])
  rho,p=blocked(d.loc[use,r.predictor],d.loc[use,r.endpoint],d.loc[use,'study'],extra)
  maxrho=max(maxrho,abs(rho-r.rho));maxp=max(maxp,abs(p-r.P));assert abs(rho-r.rho)<1e-10 and abs(p-r.P)<1e-12
 for _,g in s.groupby(['family','model']):assert np.allclose(stats.false_discovery_control(g.P.to_numpy()),g.FDR)
 ref=json.loads((old/'qc/STROMAL_CLUSTER_COMPLETE.json').read_text());st=H.stat();assert st.st_size==ref['source_bytes'] and st.st_mtime_ns==ref['source_mtime_ns']
 freeze=json.loads((P/'FREEZE.json').read_text());assert sha(P/'PROTOCOL_BEFORE_ANALYSIS.md')==freeze['protocol_sha256'] and sha(P/'PANELS_FROZEN.json')==freeze['panels_sha256']
 for r in json.loads((P/'qc/PROCESS_INPUT_MANIFEST.json').read_text()):assert sha(r['path'])==r['sha256']
 for n in ['states','process']:assert (P/'logs'/f'{n}_exit.txt').read_text().strip()=='0'
 anchors=[]
 for name,mask,gs in [('activated',act,['FAP','POSTN','CTHRC1']),('transition',trans,['CYSLTR2','LRRC17'])]:
  detected=(x[:,[idx[g] for g in gs]]>0).any(1)
  anchors.append(dict(candidate=name,cells=int(mask.sum()),with_subtype_anchor=int((mask&detected).sum()),anchor_genes=';'.join(gs)))
 tab(anchors,'STATE_CANDIDATE_ANCHOR_AUDIT.tsv')
 dump(P/'qc/STATE_PROCESS_INDEPENDENT_VALIDATION.json',dict(status='PASS',full_panel_counts_match=True,all_scores_flags_pseudobulks_match=True,C4_members=2674,source_cells=258120,process_models=len(s),max_rho_difference=maxrho,max_P_difference=maxp,all_FDR_match=True,source_and_protocol_unchanged=True))
 print('STATE_PROCESS_VALIDATION_PASS',flush=True)
if __name__=='__main__':main()
