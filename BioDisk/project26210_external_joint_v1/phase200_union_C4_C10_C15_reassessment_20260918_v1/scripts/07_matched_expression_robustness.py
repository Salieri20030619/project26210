from common import *
import sys,json
from scipy.stats import spearmanr
sys.path.insert(0,str(R/'phase157_ordered_bridge_reassessment_20260915_v1'))
from bridge_stats import bridge,rankz,corrp
B=pd.read_csv(R/'phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv',sep='\t').set_index('patient')
M=pd.read_csv(D/'bridge/SAME_LIBRARY_GROUP_METRICS.tsv',sep='\t')
X=pd.read_csv(R/'phase157_ordered_bridge_reassessment_20260915_v1/04_other_X_screen/PATIENT_CANDIDATE_VALUES.tsv',sep='\t',index_col=0)
reg=pd.read_csv(R/'phase157_ordered_bridge_reassessment_20260915_v1/04_other_X_screen/CANDIDATE_REGISTRY.tsv',sep='\t').set_index('candidate')
eligible={g:set(M.loc[(M.group==g)&(M.cells>=20)&(M.other_cells>=20),'patient']) for g in ['U','C4']}
common=eligible['U']&eligible['C4']
rows=[];direct=[]
for g in ['U','C4']:
 table=M[M.group==g].set_index('patient')
 for pop,ps in [('own_eligible',eligible[g]),('common_eligible',common)]:
  ix=B.index[B.index.isin(ps)];base=B.loc[ix];xx=X.loc[ix];tab=table.loc[ix]
  for ep in ['SPP1_log2CPM','SPP1_detection','SPP1_RNA_share']:
   yy=tab[ep]
   for adj in ['raw','study_adjusted']:
    cov=np.ones((len(ix),1)) if adj=='raw' else pd.get_dummies(base.study,dtype=float).to_numpy()
    rho,p=corrp(rankz(base.TYMP_pseudobulk),rankz(yy),cov)
    direct.append(dict(group=g,population=pop,endpoint=ep,adjustment=adj,n=len(ix),rho=rho,P=p))
    for candidate in xx:
     fit=bridge(base.TYMP_pseudobulk,xx[candidate],yy,base.study,adj)
     rows.append(dict(group=g,population=pop,endpoint=ep,adjustment=adj,candidate=candidate,family=reg.loc[candidate,'family'],**fit))
out=pd.DataFrame(rows)
for edge in ['a','XY','b']:add_fdr(out,['group','population','endpoint','adjustment','family'],edge+'_P',edge+'_FDR')
out.to_csv(D/'bridge/MATCHED_GROUP_SPP1_FIXED_CANDIDATES.tsv.gz',sep='\t',index=False)
add_fdr(pd.DataFrame(direct),['group','population','adjustment']).to_csv(D/'bridge/MATCHED_GROUP_SPP1_DIRECT_TYMP.tsv',sep='\t',index=False)
# Compare selected gene/program effects directly on identical patients and references.
core=[]
selected=['SPP1','TYMP','GPNMB','CCL18','CYP27A1','LPL','MKI67','TOP2A']
for scope in ['all','tumor']:
 for kind in ['common','shared_background']:
  u=np.load(D/'core/contrasts'/f'{scope}_U_{kind}.npz',allow_pickle=True)
  c=np.load(D/'core/contrasts'/f'{scope}_C4_{kind}.npz',allow_pickle=True)
  assert np.array_equal(u['patients'],c['patients'])
  for gene in selected:
   j=list(u['genes']).index(gene);v=u['effects'][:,j]-c['effects'][:,j]
   rng=np.random.default_rng(262204);bs=np.median(v[rng.integers(0,len(v),(2000,len(v)))],axis=1)
   core.append(dict(scope=scope,reference=kind,gene=gene,U_median=float(np.median(u['effects'][:,j])),C4_median=float(np.median(c['effects'][:,j])),CI_low=np.quantile(bs,.025),CI_high=np.quantile(bs,.975),**summary(v)))
add_fdr(pd.DataFrame(core),['scope','reference']).to_csv(D/'core/DIRECT_U_MINUS_C4_MARKER_EFFECTS.tsv',sep='\t',index=False)
# Full functional leave-one-study effects for U/C4; n remains patient count.
meta=pd.read_csv(D/'inputs/MACRO_METADATA.tsv.gz',sep='\t')
study=meta.groupby('curated_patient_id').study.first().to_dict()
los=[];compare=[]
for f in (D/'functional/patient_effects').glob('*.npz'):
 contrast,lib,version=f.stem.split('__')
 if contrast not in ['all_U','all_C4','tumor_U','tumor_C4']:continue
 z=np.load(f);studies=np.array([study[p] for p in z['patients']])
 for st in sorted(set(studies)):
  hit=studies!=st;e=z['effects'][hit]
  for j,term in enumerate(z['terms']):
   if not any(s in term.upper() for s in ['CHOLESTEROL','LYSOSOM','PHAGOS','ANTIGEN','HYPOXIA','GLYCOL','E2F','G2-M']):continue
   los.append(dict(contrast=contrast,library=lib,version=version,term=term,omitted_study=st,**summary(e[:,j])))
pd.DataFrame(los).to_csv(D/'functional/LEAVE_STUDY_FOCUS.tsv',sep='\t',index=False)
for scope in ['all','tumor']:
 for kind in ['common','shared_background']:
  for version in ['complete','without_SPP1_TYMP']:
   for lib in ['KEGG_2021_Human','MSigDB_Hallmark_2020','Reactome_2022']:
    up=np.load(D/'functional/patient_effects'/f'{scope}_U_{kind}__{lib}__{version}.npz')
    cp=np.load(D/'functional/patient_effects'/f'{scope}_C4_{kind}__{lib}__{version}.npz')
    assert np.array_equal(up['patients'],cp['patients'])
    shared=set(up['terms'])&set(cp['terms'])
    for term in sorted(shared):
     j=list(up['terms']).index(term);k=list(cp['terms']).index(term)
     a=up['effects'][:,j];b=cp['effects'][:,k]
     compare.append(dict(scope=scope,reference=kind,library=lib,version=version,term=term,U_median=float(np.median(a)),C4_median=float(np.median(b)),**summary(a-b)))
cc=add_fdr(pd.DataFrame(compare),['scope','reference','library','version'])
cc.to_csv(D/'functional/U_MINUS_C4_SHARED_PATIENT_PROGRAMS.tsv.gz',sep='\t',index=False)
# Frozen trajectory coverage and common-patient dynamic contrasts; no graph changes.
T=pd.read_csv(R/'part2_rebuild_v1/07_trajectory/monocle3_v3/ALL_CELL_PSEUDOTIME.tsv.gz',sep='\t',dtype={'candidate':str,'partition':str},low_memory=False)
coverage=[]
for (mode,st),z in T.groupby(['mode','study']):
 for group,cs in {'U':['4','10','15'],'C4':['4'],'C10':['10'],'C15':['15']}.items():
  a=z[z.candidate.isin(cs)]
  for pv in ['pt_primary','pt_root_q25','pt_root_q50','dpt_pseudotime'] if mode=='primary' else ['pt_primary']:
   coverage.append(dict(mode=mode,study=st,group=group,pseudotime=pv,cells=len(a),valid=int(np.isfinite(a[pv]).sum()),patients=a.curated_patient_id.nunique(),C15_removed_by_model=mode!='primary'))
pd.DataFrame(coverage).to_csv(D/'trajectory/ALL_GROUP_COVERAGE.tsv',sep='\t',index=False)
u=pd.read_csv(D/'trajectory/U_v2/PATIENT_DYNAMIC_TRENDS.tsv',sep='\t')
c=pd.read_csv(R/'part2_rebuild_v1/07_trajectory/dynamics_v2/PATIENT_DYNAMIC_TRENDS.tsv',sep='\t')
keys=['mode','scope','pt_version','study','patient','family','feature']
j=u.merge(c,on=keys,suffixes=('_U','_C4'))
diffs=[]
for key,z in j.groupby(['mode','scope','pt_version','family','feature']):
 for adj in ['rho','depth_adjusted_rho']:
  diffs.append(dict(zip(['mode','scope','pt_version','family','feature'],key),adjustment=adj,U_median=z[adj+'_U'].median(),C4_median=z[adj+'_C4'].median(),**summary(z[adj+'_U']-z[adj+'_C4'])))
add_fdr(pd.DataFrame(diffs),['mode','scope','pt_version','family','adjustment']).to_csv(D/'trajectory/SHARED_PATIENT_DYNAMIC_CHANGE.tsv',sep='\t',index=False)
(D/'qc/MATCHED_ROBUSTNESS_STATUS.json').write_text(json.dumps({'status':'PASS','own_eligible_expression_patients':{k:len(v) for k,v in eligible.items()},'shared_expression_patients':len(common),'fixed_candidate_expression_models':len(out),'shared_gene_effect_comparisons':len(core),'shared_program_comparisons':len(cc),'no_graph_refit':True},indent=2))
print('MATCHED_ROBUSTNESS_COMPLETE',flush=True)
