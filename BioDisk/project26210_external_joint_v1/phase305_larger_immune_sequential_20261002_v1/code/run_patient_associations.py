from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='4'
import json,subprocess,datetime,gc
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import rankdata,binomtest,t as tdist
P=Path(__file__).resolve().parents[1];R=P.parent
context=json.loads((R/'phase280_TYMP_source_niche_function_20260928_v1/resources/CONTEXT_PROGRAMS.json').read_text())['genes'];context={k:[g for g in gs if g not in ['TYMP','SPP1']] for k,gs in context.items()}
stroma=json.loads((R/'phase292_barrier_execution_and_external_feasibility_20260930_v1/PANELS_FROZEN.json').read_text());im=pd.read_csv(R/'phase213_immune_restart_execution_20260919_v1/I1_patient_inputs/PROGRAM_GENE_COVERAGE.tsv',sep='\t');immune={k:im.loc[im.endpoint.eq(k),'gene'].unique().tolist() for k in ['CD8_cytotoxicity','CD8_dysfunction']}
identity=json.loads((P/'resources/LINEAGE_IDENTITY_FROZEN.json').read_text());gates=identity['n_gates']
protocol={'created':datetime.datetime.now().astimezone().isoformat(),'primary':'OBRADOVIC untreated CD45+ mature-macrophage within-patient TYMP-SPP1 rank association; no new C4 classification','within_models':['depth_complexity_library_subcluster','depth_complexity_library_subcluster_IFN_hypoxia'],'within_gate':'>=100 qualified macrophages; TYMP and SPP1 each >=10 expressing cells; residual dimensions >=8. No cell-level P.','within_inference':'Two-sided exact binomial sign test on independent patients if n>=10; models form BH family of2. n<10 descriptive effects/directions only. HCRN baseline and progression separate; donor-first-available summary is descriptive sensitivity, never duplicates a donor.','aggregation':'Sum raw RNA counts per patient and treatment stratum within qualified lineage; log1p CP10K; fixed-program equal gene mean. OBRADOVIC CD45+ and CD45- fractions paired by patient, not the same sequencing library.','single_gene_gate':'gene measured in every contributing library and >=5 expressing cells','program_gate':'at least80% fixed genes measured in every contributing library, at least50% full panel detected in >=5% lineage cells; no missing gene treated as measured zero','programs':stroma|immune,'context':context,'across_models':['unadjusted','prespecified_adjusted'],'across_adjustment':'OBRADOVIC ranked response/predictor controlled for author processing batch1 versus batch2 and author pT1b versus pT3a. Stored batch labels batch1_v2/batch2_v3 reflect matrix-format groups, not verified 10x chemistry. HCRN evaluated within each timepoint, no cross-stratum pooling. Constant design columns have no effect.','across_inference':'Partial Spearman rho, two-sided t approximation df=N-rank(D)-1 when N>=10 and df>=5; otherwise effects only. Families per cohort,stratum,source,receiver,stage include planned endpoints and both models: stroma3x2,immune2x2,upstream1x2. BH on evaluable hypotheses; report rawP and FDR; rawP<.05 nominal screening.','HCRN_donor_first':'Baseline preferred regardless of cell or target eligibility; otherwise Progression; not a treatment-naive cohort. No cross-cohort pooling.','major_limit':'RNA programs measure transcription, not migration or physical exclusion; tumour-like identity is transcriptomic, not newly CNV-validated.'}
fr=P/'resources/PATIENT_ASSOCIATION_PROTOCOL.json'
if fr.exists():
 previous=json.loads(fr.read_text());assert {k:v for k,v in previous.items() if k!='created'}=={k:v for k,v in protocol.items() if k!='created'},'Protocol changed: version required'
else:fr.write_text(json.dumps(protocol,indent=2))
def design(n,num=(),cat=()):
 cols=[np.ones(n)]
 for v in num:
  if len(v) and np.ptp(v)>0:cols.append(rankdata(v))
 for v in cat:
  d=pd.get_dummies(pd.Series(v).astype(str),drop_first=True,dtype=float)
  cols.extend(d.to_numpy().T)
 return np.column_stack(cols)
def corr(x,y,D):
 xy=np.column_stack([rankdata(x),rankdata(y)]);res=xy-D@np.linalg.lstsq(D,xy,rcond=None)[0];rank=int(np.linalg.matrix_rank(D));rho=float(np.corrcoef(res.T)[0,1]) if np.std(res[:,0])>1e-9 and np.std(res[:,1])>1e-9 else np.nan
 return rho,rank
def bh(df,groups):
 df['FDR']=np.nan;df['family_tested_n']=0
 for key,idx in df.groupby(groups,dropna=False).groups.items():
  ids=[i for i in idx if pd.notna(df.at[i,'P'])];m=len(ids);df.loc[list(idx),'family_tested_n']=m
  if m:
   order=sorted(ids,key=lambda i:df.at[i,'P']);p=df.loc[order,'P'].to_numpy();q=np.minimum.accumulate((p*m/np.arange(1,m+1))[::-1])[::-1];df.loc[order,'FDR']=np.clip(q,0,1)
 return df
measures=[];coverage=[];within=[];caches=[]
for ds in ['OBRADOVIC','HCRN']:
 subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True,stdout=subprocess.DEVNULL)
 a=ad.read_h5ad(P/f'inputs/{ds}_counts_qc.h5ad');meta=pd.read_csv(P/f'tables/{ds}_QUALIFIED_CELL_METADATA.tsv.gz',sep='\t',index_col=0);assert a.obs_names.equals(meta.index);x=a.X.tocsr().astype(np.float64);norm=x.multiply((1e4/np.asarray(x.sum(axis=1)).ravel())[:,None]).tocsr();norm.data=np.log1p(norm.data)
 present={}
 for lib in meta.library.unique():
  b=ad.read_h5ad(P/f'inputs/{ds}_{lib}_qc.h5ad',backed='r');present[lib]=set(b.var_names);b.file.close()
 for unit,g in meta.groupby('analysis_unit',observed=True):
  ui=meta.analysis_unit.eq(unit).to_numpy();patient=str(g.patient.iloc[0]);pt3=patient.split(':')[-1] in ['PatientA','Patient2','Patient3','Patient5','Patient7']
  base={'dataset':ds,'analysis_unit':unit,'patient':patient,'stratum':str(g.stratum.iloc[0]),'technical_batch':str(g.technical_batch.iloc[0]),'pT_category':('pT3a' if pt3 else 'pT1b') if ds=='OBRADOVIC' else 'not_used','donor_first_available':bool(g.donor_first_available.iloc[0])}
  for lin,minn in gates.items():
   sel=ui&meta['lineage_'+lin].to_numpy(dtype=bool);n=int(sel.sum());libraries=meta.loc[sel,'library'].unique().tolist();available=set.intersection(*[present[l] for l in libraries]) if libraries else set();summ=np.asarray(x[sel].sum(axis=0)).ravel();pb=np.log1p(1e4*summ/max(summ.sum(),1));det=np.asarray((x[sel]>0).sum(axis=0)).ravel()
   panels={'TYMP':['TYMP'],'SPP1':['SPP1']} if lin=='Macrophage' else ({'TYMP':['TYMP']} if lin=='Tumour_like' else (immune if lin=='CD8' else stroma))
   for name,gs in panels.items():
    keep=[z for z in gs if z in available];jj=a.var_names.get_indexer(keep);score=float(pb[jj].mean()) if len(jj) else np.nan
    detectable=int(det[jj].sum()) if len(gs)==1 else int((det[jj]>=max(1,.05*n)).sum());ok=n>=minn and len(keep)>=.8*len(gs) and (detectable>=5 if len(gs)==1 else detectable>=.5*len(gs))
    reason=[]
    if n<minn:reason.append('CELL_COUNT_GATE')
    if len(keep)<.8*len(gs):reason.append('GENE_MEASUREMENT_GATE')
    if not (detectable>=5 if len(gs)==1 else detectable>=.5*len(gs)):reason.append('DETECTION_GATE')
    measures.append(base|{'lineage':lin,'measure':name,'cells':n,'value':score,'eligible':ok,'reason':';'.join(reason),'genes_expected':len(gs),'genes_measured':len(keep),'detected_cells_or_genes':detectable,'libraries':';'.join(libraries),'median_UMI':float(meta.loc[sel,'total_counts'].median()) if n else np.nan})
    for gene,jj0 in zip(keep,jj):coverage.append(base|{'lineage':lin,'measure':name,'gene':gene,'cells':n,'detected_cells':int(det[jj0]),'pseudobulk_logCP10K':float(pb[jj0])})
   if lin!='Macrophage':continue
   target_ok=all(z in available for z in ['TYMP','SPP1']);vals={z:norm[sel,a.var_names.get_loc(z)].toarray().ravel() for z in ['TYMP','SPP1']} if target_ok else {};good=n>=100 and target_ok and all((v>0).sum()>=10 for v in vals.values())
   numeric=[meta.loc[sel,'total_counts'].to_numpy(),meta.loc[sel,'n_genes'].to_numpy()];cats=[meta.loc[sel,'library'].to_numpy(),meta.loc[sel,'myeloid_cluster'].to_numpy()];env=[];env_ok=all(len(set(gs)&available)>=.8*len(gs) for gs in context.values())
   if env_ok:
    for gs in context.values():env.append(np.asarray(norm[sel][:,a.var_names.isin(set(gs)&available)].mean(axis=1)).ravel())
   cache={'TYMP':vals.get('TYMP',np.array([])),'SPP1':vals.get('SPP1',np.array([]))}
   for model in protocol['within_models']:
    D=design(n,numeric+(env if model.endswith('hypoxia') else []),cats);rho=np.nan;rank=int(np.linalg.matrix_rank(D)) if n else 0;reason=''
    if not good:reason='CELL_OR_TARGET_DETECTION_GATE'
    elif model.endswith('hypoxia') and not env_ok:reason='ENVIRONMENT_GENE_COVERAGE_GATE'
    elif n-rank<8:reason='RESIDUAL_DIMENSION_GATE'
    else:rho,rank=corr(vals['TYMP'],vals['SPP1'],D);reason='' if np.isfinite(rho) else 'NO_RESIDUAL_VARIATION'
    within.append(base|{'model':model,'cells':n,'TYMP_detected':int((vals.get('TYMP',np.array([]))>0).sum()),'SPP1_detected':int((vals.get('SPP1',np.array([]))>0).sum()),'rho':rho,'design_rank':rank,'reason':reason,'P':np.nan})
    cache[model+'_D']=D
   cachefile=P/'inputs'/('WITHIN_'+ds+'_'+unit.replace(':','_')+'.npz');np.savez_compressed(cachefile,**cache);caches.append({'analysis_unit':unit,'path':str(cachefile)})
 print('PATIENT_INPUTS_DONE',ds,flush=True);del a,x,norm;gc.collect()
m=pd.DataFrame(measures);m.to_csv(P/'tables/PATIENT_MEASURES.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(P/'tables/PROGRAM_GENE_COVERAGE.tsv',sep='\t',index=False);w=pd.DataFrame(within);w.to_csv(P/'tables/WITHIN_MACRO_PATIENT_EFFECTS.tsv',sep='\t',index=False);(P/'inputs/WITHIN_CACHE_MANIFEST.json').write_text(json.dumps(caches,indent=2))
ws=[]
for ds in ['OBRADOVIC','HCRN']:
 scopes=['Untreated'] if ds=='OBRADOVIC' else ['Baseline','Progression','donor_first_available']
 for scope in scopes:
  selected=w[w.dataset.eq(ds)&(w.donor_first_available if scope=='donor_first_available' else w.stratum.eq(scope))]
  for model,g in selected.groupby('model'):
   g=g[g.rho.notna()];assert g.patient.nunique()==len(g);n=len(g);positive=int((g.rho>0).sum());negative=int((g.rho<0).sum());pval=binomtest(positive,positive+negative,.5).pvalue if n>=10 and positive+negative>=10 else np.nan
   ws.append({'dataset':ds,'scope':scope,'model':model,'patients':n,'positive':positive,'negative':negative,'median_rho':float(g.rho.median()) if n else np.nan,'min_rho':float(g.rho.min()) if n else np.nan,'max_rho':float(g.rho.max()) if n else np.nan,'P':pval,'family':'within_TYMP_SPP1_two_models','family_planned_n':2})
ws=bh(pd.DataFrame(ws),['dataset','scope','family']);ws.to_csv(P/'tables/WITHIN_MACRO_SUMMARY.tsv',sep='\t',index=False)
(P/'qc/STEP1_WITHIN_COMPLETE.json').write_text(json.dumps({'status':'COMPLETE'}))
edges=[];pairs=[]
def edge(stage,source,target,family,endpoints):
 for ds,scope in [('OBRADOVIC','Untreated'),('HCRN','Baseline'),('HCRN','Progression')]:
  d=m[m.dataset.eq(ds)&m.stratum.eq(scope)&m.eligible];cols=['analysis_unit','patient','technical_batch','pT_category','value']
  aa=d[d.lineage.eq(source[0])&d.measure.eq(source[1])][cols];bb=d[d.lineage.eq(target[0])&d.measure.eq(target[1])][['analysis_unit','patient','value']];g=aa.merge(bb,on=['analysis_unit','patient'],suffixes=('_source','_target'),validate='one_to_one');n=len(g);assert g.patient.nunique()==n
  for model in protocol['across_models']:
   cats=[g.technical_batch.to_numpy(),g.pT_category.to_numpy()] if ds=='OBRADOVIC' and model=='prespecified_adjusted' else []
   D=design(n,cat=cats);rho=np.nan;rank=int(np.linalg.matrix_rank(D)) if n else 0;df=n-rank-1
   if n>=3 and g.value_source.nunique()>1 and g.value_target.nunique()>1 and df>=1:rho,rank=corr(g.value_source,g.value_target,D)
   pval=float(2*tdist.sf(abs(rho)*np.sqrt(df/max(1-rho*rho,1e-15)),df)) if n>=10 and df>=5 and np.isfinite(rho) else np.nan
   rec={'stage':stage,'dataset':ds,'scope':scope,'source':':'.join(source),'target':':'.join(target),'model':model,'patients':n,'rho':rho,'design_rank':rank,'df':df,'P':pval,'family':family,'family_planned_n':endpoints*2,'status':'inferential' if pd.notna(pval) else ('descriptive' if pd.notna(rho) else 'not_evaluable')};edges.append(rec)
   for row in g.to_dict('records'):pairs.append(rec|row)
edge(1,('Tumour_like','TYMP'),('Macrophage','SPP1'),'tumour_TYMP_macro_SPP1',1)
for lin in ['Fibro_strict_candidate','Stromal_comarker','Mural']:
 for src in [('Tumour_like','TYMP'),('Macrophage','TYMP'),('Macrophage','SPP1')]:
  for ep in stroma:edge(2,src,(lin,ep),'stroma_'+lin+'_'+'_'.join(src),3)
for src in [('Tumour_like','TYMP'),('Macrophage','TYMP'),('Macrophage','SPP1')]+[(lin,'collagen_maturation') for lin in ['Fibro_strict_candidate','Stromal_comarker','Mural']]:
 for ep in immune:edge(3,src,('CD8',ep),'immune_'+'_'.join(src),2)
e=bh(pd.DataFrame(edges),['dataset','scope','family']);e.to_csv(P/'tables/ALL_PATIENT_ASSOCIATIONS.tsv',sep='\t',index=False);pd.DataFrame(pairs).to_csv(P/'tables/ASSOCIATION_PATIENT_PAIRS.tsv',sep='\t',index=False)
for step in [1,2,3]:e[e.stage.eq(step)].to_csv(P/f'tables/STEP{step}_PATIENT_ASSOCIATIONS.tsv',sep='\t',index=False);(P/f'qc/STEP{step}_ASSOCIATIONS_COMPLETE.json').write_text(json.dumps({'status':'COMPLETE'}))
print('ASSOCIATIONS_COMPLETE',len(e),flush=True)
