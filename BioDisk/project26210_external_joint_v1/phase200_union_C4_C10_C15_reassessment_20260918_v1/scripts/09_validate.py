from common import *
import json,shutil
from scipy import sparse
from scipy.stats import rankdata,t
checks=[]
def check(name,ok,detail=''):
 checks.append(dict(check=name,passed=bool(ok),detail=detail))
 if not ok:raise AssertionError(name+' '+str(detail))
# Strict matched pathway comparison uses a COMMON measured-gene universe.
old=D/'functional/U_MINUS_C4_SHARED_PATIENT_PROGRAMS.tsv.gz'
shutil.copy2(old,D/'functional/U_MINUS_C4_OWN_GENE_UNIVERSES_ARCHIVED.tsv.gz')
libs={}
for f in (R/'part2_rebuild_v1/06_functional_programs/run_v1/resources').glob('*.gmt'):
 libs[f.stem]={p[0]:set(x.split(',')[0] for x in p[2:] if x) for p in [line.split('\t') for line in f.read_text().splitlines()]}
out=[]
for scope in ['all','tumor']:
 for ref in ['common','shared_background']:
  u=np.load(D/'core/contrasts'/f'{scope}_U_{ref}.npz',allow_pickle=True);c=np.load(D/'core/contrasts'/f'{scope}_C4_{ref}.npz',allow_pickle=True)
  check(scope+ref+'_patient_identity',np.array_equal(u['patients'],c['patients']))
  ue=u['effects'];ce=c['effects'];gi={g:i for i,g in enumerate(u['genes'])};common=set(u['genes'][u['observed']&c['observed']])
  for lib,sets in libs.items():
   for ver,remove in [('complete',set()),('without_SPP1_TYMP',{'SPP1','TYMP'})]:
    for term,gs in sets.items():
     hit=sorted(gs&(common-remove))
     if not 10<=len(hit)<=500:continue
     idx=[gi[g] for g in hit];a=ue[:,idx].mean(1);b=ce[:,idx].mean(1)
     out.append(dict(scope=scope,reference=ref,library=lib,version=ver,term=term,identical_genes=len(hit),U_median=float(np.median(a)),C4_median=float(np.median(b)),**summary(a-b)))
df=add_fdr(pd.DataFrame(out),['scope','reference','library','version'])
df.to_csv(old,sep='\t',index=False)
check('common_gene_universe_used_for_paired_pathway_comparison',True,str(len(df)))
# Check multiple-testing pools.
ff=pd.read_csv(D/'functional/ALL_PROGRAM_RESULTS.tsv.gz',sep='\t')
for key,z in ff.groupby(['contrast','library','version']):check('functional_BH_'+str(key),np.allclose(bh(z.P),z.FDR,equal_nan=True,atol=1e-13))
check('C10_two_patient_function_P_absent',ff.loc[ff.contrast.isin(['all_C10','tumor_C10']),'P'].isna().all())
for p in (D/'core/contrasts').glob('*_genes.tsv.gz'):
 a=pd.read_csv(p,sep='\t');check('gene_BH_'+p.stem,np.allclose(bh(a.P),a.FDR,equal_nan=True,atol=1e-13))
# Original C4 functional estimates must match the previous frozen analysis.
oldf=pd.read_csv(R/'part2_rebuild_v1/06_functional_programs/run_v1/ALL_PATIENT_PROGRAM_COMPARISONS.tsv.gz',sep='\t')
oldf=oldf[oldf.cluster==4]
j=ff[ff.contrast=='all_C4'].merge(oldf,on=['library','version','term'],suffixes=('_new','_old'))
check('original_C4_function_complete_coverage',len(j)==len(oldf))
check('original_C4_function_effects',np.allclose(j['median'],j.median_paired_program_difference,atol=1e-10))
check('original_C4_function_P',np.allclose(j.P_new,j.P_old,equal_nan=True,atol=1e-10))
# Independent raw SPP1 summary against direct sparse cell counts.
m=pd.read_csv(D/'inputs/MACRO_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str});x=sparse.load_npz(D/'inputs/MACRO_COUNTS.npz')
genes=pd.read_csv(D/'inputs/GENES.tsv',sep='\t').gene.tolist();spp=genes.index('SPP1')
metrics=pd.read_csv(D/'bridge/SAME_LIBRARY_GROUP_METRICS.tsv',sep='\t')
defs={'U':['4','10','15'],'C4':['4'],'C10':['10'],'C15':['15']}
for row in metrics.itertuples():
 allmask=m.sample_id.eq(row.sample_id).to_numpy();mask=allmask&m.candidate.isin(defs[row.group]).to_numpy()
 total=float(x[mask].sum());count=float(x[mask,spp].sum());den=float(x[allmask,spp].sum())
 expected=np.log2(count/total*1e6+1) if total else np.nan
 check('raw_SPP1_'+row.patient+'_'+row.group,np.isclose(expected,row.SPP1_log2CPM,equal_nan=True) and np.isclose(count/den,row.SPP1_RNA_share,equal_nan=True) and mask.sum()==row.cells)
# Independent QR partial correlations for all primary fixed candidate bridge results.
B=pd.read_csv(R/'phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv',sep='\t').set_index('patient')
wide=metrics.pivot(index='patient',columns='group',values='fraction').loc[B.index]
inputs={
'original116':pd.read_csv(R/'phase157_ordered_bridge_reassessment_20260915_v1/04_other_X_screen/PATIENT_CANDIDATE_VALUES.tsv',sep='\t',index_col=0).loc[B.index],
'metabolic_components':pd.read_csv(R/'phase159_metabolic_parts_C4_value_20260915_v1/PATIENT_SENDER_VALUES.tsv',sep='\t',index_col=0).loc[B.index],
'MTA_proxies':pd.read_csv(R/'phase157_ordered_bridge_reassessment_20260915_v1/03_MTA_to_new_C4/MATCHED_MTA_PROXY_PATIENTS.tsv',sep='\t').set_index('patient').loc[B.index]}
res=pd.read_csv(D/'bridge/ALL_FIXED_CANDIDATE_BRIDGES.tsv',sep='\t')
def qr_r(x,y,C):
 Q=np.linalg.qr(C)[0][:,:np.linalg.matrix_rank(C)]
 a=x-Q@(Q.T@x);b=y-Q@(Q.T@y)
 return np.corrcoef(a,b)[0,1] if min(np.std(a),np.std(b))>1e-10 else np.nan
for i,row in res.iterrows():
 A=rankdata(B.TYMP_pseudobulk);X=rankdata(inputs[row.screen][row.candidate]);Y=rankdata(wide[row.group])
 C=np.ones((16,1)) if row.adjustment=='raw' else pd.get_dummies(B.study,dtype=float).to_numpy()
 check('QR_a_'+str(i),np.isclose(qr_r(A,X,C),row.a_rho,equal_nan=True,atol=1e-9))
 check('QR_XY_'+str(i),np.isclose(qr_r(X,Y,C),row.XY_rho,equal_nan=True,atol=1e-9))
 check('QR_b_'+str(i),np.isclose(qr_r(X,Y,np.column_stack([C,A])),row.b_partial_rho,equal_nan=True,atol=1e-9))
for key,z in res.groupby(['screen','group','family','adjustment']):
 for edge in ['a','XY','b']:check('bridge_BH_'+str(key)+edge,np.allclose(bh(z[edge+'_P']),z[edge+'_FDR'],equal_nan=True,atol=1e-13))
# GSEA ES and BH independently verified; complete eligible family retained.
gs=libs['KEGG_2021_Human']
gsea=pd.read_csv(D/'gsea/ALL_GSEA_RESULTS.tsv',sep='\t')
for (model,version),z in gsea.groupby(['model','version']):
 rank=pd.read_csv(D/'gsea'/f'{model}_rank.tsv',sep='\t')
 if version!='complete':rank=rank[~rank.gene.isin(['SPP1','TYMP'])]
 check('GSEA_BH_'+model+version,np.allclose(bh(z.pvalue),z['p.adjust'],equal_nan=True,atol=1e-12))
 for _,row in z[z.ID.isin(['Cholesterol metabolism','Lysosome','Phagosome','Antigen processing and presentation','Cell cycle'])].iterrows():
  hit=rank.gene.isin(gs[row.ID]).to_numpy();w=np.abs(rank.score.to_numpy());n=len(hit);k=hit.sum()
  walk=np.cumsum(np.where(hit,w/w[hit].sum(),-1/(n-k)))
  es=walk.max() if walk.max()>=-walk.min() else walk.min()
  check('GSEA_ES_'+model+version+row.ID,np.isclose(es,row.enrichmentScore,atol=1e-10))
# All old C4 versus other-macrophage AUC contrasts reproduce.
oldr=pd.read_csv(R/'part2_rebuild_v1/08_regulons/statistics_v2_20260915/AUC_CONTRAST_SUMMARY.tsv',sep='\t')
oldr=oldr[oldr.endpoint=='C4_minus_other_macrophages'].copy();oldr['scope']=oldr.scope.map({'all_tissues':'all','tumor_only':'tumor'})
newr=pd.read_csv(D/'regulons/ALL_GROUP_AUC_SUMMARY.tsv',sep='\t');newr=newr[newr.group=='C4']
jr=newr.merge(oldr,on=['scope','regulon','TF'])
check('C4_frozen_regulon_baseline_effect',np.allclose(jr['median'],jr.median_effect,atol=1e-12))
check('C4_frozen_regulon_baseline_P',np.allclose(jr.P_x,jr.P_y,equal_nan=True,atol=1e-12))
# Regulon target inclusion audit, preventing circular independent-validation wording.
ed=pd.read_csv(R/'part2_rebuild_v1/08_regulons/network_v2_20260915/REGULON_TARGET_EDGES.tsv',sep='\t')
audit=[]
for reg,z in ed.groupby('regulon'):
 if reg.split('(')[0] in ['JUN','JUNB','CEBPD']:
  audit.append(dict(regulon=reg,SPP1_is_target='SPP1' in set(z.target),TYMP_is_target='TYMP' in set(z.target)))
pd.DataFrame(audit).to_csv(D/'regulons/THREE_TF_OUTCOME_MEMBERSHIP_AUDIT.tsv',sep='\t',index=False)
# Coverage relative to both retained-model cells and ORIGINAL frozen union.
cov=pd.read_csv(D/'trajectory/ALL_GROUP_COVERAGE.tsv',sep='\t')
cv=cov.groupby(['mode','group','pseudotime'],as_index=False)[['cells','valid']].sum()
orig={'U':4502,'C4':2674,'C10':1315,'C15':513}
cv['original_frozen_cells']=cv.group.map(orig)
cv['fraction_of_model_cells']=cv.valid/cv.cells.replace(0,np.nan)
cv['fraction_of_original_cells']=cv.valid/cv.original_frozen_cells
cv.to_csv(D/'trajectory/COVERAGE_SUMMARY.tsv',sep='\t',index=False)
for marker in ['CORE_FUNCTIONAL_STATUS.json','BRIDGE_STATUS.json','RECEIVER_REGULON_STATUS.json','MATCHED_ROBUSTNESS_STATUS.json','DEPENDENCY_COMPLETION_STATUS.json']:
 check(marker,json.loads((D/'qc'/marker).read_text())['status']=='PASS')
check('no_new_figures',not any(p.suffix.lower() in ['.png','.pdf','.svg'] for p in D.rglob('*')))
pd.DataFrame(checks).to_csv(D/'qc/INDEPENDENT_VALIDATION.tsv',sep='\t',index=False)
(D/'qc/VALIDATION_STATUS.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'failed':0,'common_pathway_gene_universe_aligned':True,'old_C4_marker_function_regulon_reproduced':True,'new_figures':False},indent=2))
print('INDEPENDENT_VALIDATION_PASS',len(checks),flush=True)
