from common import *
from spatial import CORE,ORG,IMM
from process import NEW,LINKS,CTX
def ref(x,y,cv):
 A=pd.DataFrame(np.column_stack([x,y])).rank(method='average').to_numpy()
 c=pd.DataFrame(cv).rank(method='average').to_numpy();sd=c.std(0);D=np.column_stack([np.ones(len(A)),(c-c.mean(0))/sd])
 U=np.linalg.svd(D,full_matrices=False)[0];e=A-U@(U.T@A)
 return float(stats.pearsonr(e[:,0],e[:,1]).statistic)
def main():
 for name in ['spatial','process']:assert (P/'logs'/f'{name}_exit.txt').read_text().strip()=='0'
 checks=[];cache=read(P/'inputs/PROCESS_CELL_CACHE.tsv.gz',dtype={'candidate':str});e=read(P/'tables/PROCESS_PATIENT_EFFECTS.tsv');mx=0.;count=0
 for (pat,scope,model),g in e[e.effect.notna()].groupby(['patient','scope','model']):
  d=cache[cache.patient==pat]
  if scope=='C4':d=d[d.candidate=='4']
  cc=[d.log_depth.to_numpy(),d.n_genes.to_numpy()]
  if scope=='WholeMacro':cc += list(pd.get_dummies(d.candidate,dtype=float,drop_first=True).to_numpy().T)
  if model=='context':cc += [d[x].to_numpy() for x in CTX]
  cv=np.column_stack(cc)
  for r in g.itertuples():
   value=ref(d[r.source],d[r.target],cv);assert abs(value-r.effect)<1e-10,(r,value)
   mx=max(mx,abs(value-r.effect));count+=1
 checks.append(dict(check='process_rhos_independent_SVD',n=count,max_error=mx))
 configs={'unadjusted':[],'core_adjusted':CORE,'inside_TYMP_adjusted':CORE+['inside_TYMP'],'vascular_stromal_sensitivity':CORE+['inside_TYMP','inside_endothelial','outside_endothelial','fibroblast_mean'],'inside_TYMP_and_coverage':CORE+['inside_TYMP','SPP1_POSTN_coverage']}
 sec=read(P/'tables/OUTER_TYMP_SECTION_EFFECTS.tsv');mx=0.;count=0
 for sample,g in sec[sec.effect.notna()].groupby('sample'):
  z=read(PRE/'inputs'/f'{sample}_SPLIT_INTERFACE.tsv.gz')
  for r in g.itertuples():
   cols=['outside_TYMP',r.endpoint]+configs[r.model];d=z.loc[z.split_eligible,cols].dropna()
   if configs[r.model]:value=ref(d.outside_TYMP,d[r.endpoint],d[configs[r.model]].to_numpy())
   else:value=float(stats.spearmanr(d.outside_TYMP,d[r.endpoint]).statistic)
   assert abs(value-r.effect)<1e-10;mx=max(mx,abs(value-r.effect));count+=1
 checks.append(dict(check='spatial_rhos_independent_SVD',n=count,max_error=mx))
 pa=read(P/'tables/OUTER_TYMP_PATIENT_EFFECTS.tsv');keys=['dataset','family','endpoint','model','patient']
 for row in pa.to_dict('records'):
  mask=np.ones(len(sec),bool)
  for k in keys:mask &= sec[k].eq(row[k]).to_numpy()
  vals=sec.loc[mask,'effect'].dropna().to_numpy()
  if len(vals):
   x=np.arctanh(np.clip(vals,-.99999999,.99999999));value=(np.exp(2*x.mean())-1)/(np.exp(2*x.mean())+1);assert abs(value-row['effect'])<1e-10
  else:assert pd.isna(row['effect'])
 checks.append(dict(check='patient_section_aggregation',n=len(pa)))
 for sf,ef,keys,family in [('PROCESS_SUMMARY.tsv','PROCESS_PATIENT_EFFECTS.tsv',['scope','source','target','model'],['scope']),('OUTER_TYMP_COHORT_SUMMARY.tsv','OUTER_TYMP_PATIENT_EFFECTS.tsv',['dataset','family','endpoint','model'],['dataset','family'])]:
  su=read(P/'tables'/sf);ee=read(P/'tables'/ef)
  for r in su.to_dict('records'):
   mask=np.ones(len(ee),bool)
   for k in keys:mask &= ee[k].eq(r[k]).to_numpy()
   v=ee.loc[mask,'effect'].dropna().to_numpy();assert len(v)==r['n'] and int((v>0).sum())==r['positive']
   if len(v):assert abs(np.median(v)-r['median'])<1e-12
   if len(v)>=10:
    pp=stats.binomtest(int((v>0).sum()),int((v!=0).sum()),.5).pvalue if np.any(v!=0) else 1
    assert abs(pp-r['P'])<1e-12
   else:assert pd.isna(r['P'])
  for _,d in su.groupby(family):
   q=stats.false_discovery_control(d.P.fillna(1).to_numpy());ok=d.P.notna().to_numpy();assert np.allclose(q[ok],d.FDR.to_numpy()[ok],atol=1e-12)
  checks.append(dict(check='P_BH_'+sf,n=len(su)))
 # Reconstruct all same-library selected counts and all target correlations.
 base=read(P/'tables/MALIGNANT_TO_MACRO_PATIENT_INPUTS.tsv');S=R/'phase157_ordered_bridge_reassessment_20260915_v1/02_inputs';G=read(S/'GROUP_INDEX.tsv');A=np.load(S/'GROUP_PSEUDOBULK.npz',allow_pickle=True);genes=list(A['genes'])
 for r in base.itertuples():
  for comp,roles in [('malignant',['Malignant']),('macro',['C4','OtherMacrophage'])]:
   rows=G[(G.patient==r.patient)&(G.sample_id==r.library)&G.role.isin(roles)];v=A['counts'][rows.group_index.to_numpy(int)].sum(0);lp=np.log2(v/v.sum()*1e6+1)
   assert abs(lp[genes.index('TYMP')]-getattr(r,comp+'_TYMP'))<1e-10
   if comp=='macro':
    for gene in NEW:assert abs(lp[genes.index(gene)]-getattr(r,gene))<1e-10
 su=read(P/'tables/MALIGNANT_TO_MACRO_SUMMARY.tsv');mx=0.
 for r in su.itertuples():
  d=base;cols=[np.ones(len(d))]+[d.study.eq(s).to_numpy(float) for s in sorted(d.study.unique())[1:]]
  if r.model=='study_macroTYMP':cols.append(d.macro_TYMP.rank().to_numpy())
  D=np.column_stack(cols);U=np.linalg.svd(D,full_matrices=False)[0];V=d[['malignant_TYMP',r.target]].rank().to_numpy();V=V-U@(U.T@V);rho=stats.pearsonr(V[:,0],V[:,1]).statistic
  assert abs(rho-r.rho)<1e-10;mx=max(mx,abs(rho-r.rho))
 assert np.allclose(stats.false_discovery_control(su.P.fillna(1)),su.FDR)
 checks.append(dict(check='same_library_128_gene_values_and_12_rhos',n=len(su),max_error=mx))
 # Validate newly rebuilt environment programs against original counts.
 cfg=json.loads((P/'resources/PROCESS_CONTEXT_PROGRAMS.json').read_text());assert all(not(set(gs)&set(cfg['excluded'])) for gs in cfg['programs'].values())
 m=read(M/'inputs/MACRO_METADATA.tsv.gz',index_col=0);X=sparse.load_npz(M/'inputs/MACRO_COUNTS.npz').tocsr();gene=read(M/'inputs/GENES.tsv').gene.tolist();gi={g:i for i,g in enumerate(gene)};total=np.asarray(X.sum(1)).ravel();mx=0.;nv=0
 for pat,d in cache.groupby('patient'):
  ids=m.index.get_indexer(d.cell_id)
  for name,gs in cfg['programs'].items():
   a=X[ids][:,[gi[g] for g in gs]].toarray().astype(float);v=pd.DataFrame(np.log1p(a/total[ids,None]*1e4));sd=v.std(ddof=0).replace(0,1);score=((v-v.mean())/sd).mean(1).to_numpy()
   er=np.max(np.abs(score-d[name].to_numpy()));assert er<1e-10;mx=max(mx,float(er));nv+=len(d)
 checks.append(dict(check='disjoint_environment_scores',n=nv,max_error=mx))
 for file in ['PROCESS_INPUTS.json','SPATIAL_INPUTS.json']:
  for rec in json.loads((P/'qc'/file).read_text()):assert sha(rec['path'])==rec['sha256']
 assert sha(P/'PROTOCOL_BEFORE_TESTING.md')==json.loads((P/'FREEZE.json').read_text())['protocol_sha256']
 # Archive evidence from raw cell annotation confirming no independent CAF cohort.
 q={'source':'part1_scvi_integrated_cell_annotated_cnv_v1.h5ad obs/cell_type_broad_v1','tumor_fibroblast_cells':39,'patients':1,'patient':'GSE156632__RCC1','annotation_confidence':'Medium','use':'descriptive_asset_only_no_inferential_CAF_response'}
 dump(P/'qc/FIBROBLAST_ASSET_LIMIT.json',q)
 dump(P/'qc/FINAL_VALIDATION.json',dict(status='PASS',checks=checks,source_files_unchanged=True,figures=0,software=dict(python=sys.version,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__)))
 print('VALIDATION_PASS',checks,flush=True)
if __name__=='__main__':main()
