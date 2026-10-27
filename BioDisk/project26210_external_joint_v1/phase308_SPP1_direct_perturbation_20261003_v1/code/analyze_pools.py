from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='4'
import json,datetime,hashlib,itertools
import numpy as np,pandas as pd,anndata as ad
P=Path(__file__).resolve().parents[1];R=P.parent
frozen=json.loads((P/'resources/LINEAGE_IDENTITY_FROZEN.json').read_text());a=ad.read_h5ad(P/frozen['input_object']);x=a.X.tocsr();meta=a.obs.copy();names=a.var_names;groups=frozen['cluster_labels'];meta['broad_lineage']=meta.cluster.astype(str).map(groups);assert meta.broad_lineage.notna().all()
def n_support(gs):
 jj=names.get_indexer([g for g in gs if g in names]);return np.asarray((x[:,jj]>0).sum(1)).ravel()
t_support=n_support(['Cd3d','Cd3e','Trac']);my_support=n_support(['Lyz2','Tyrobp','Lst1','Fcer1g']);ep_support=n_support(['Epcam','Krt8','Krt18','Krt19']);en_support=n_support(['Pecam1','Cdh5','Vwf','Esam']);fib_support=n_support(['Col1a1','Col1a2','Col3a1','Dcn','Lum','Pdgfra']);mural_support=n_support(['Rgs5','Cspg4','Mcam','Myh11','Des']);macro_support=n_support(['C1qa','C1qb','C1qc','Adgre1','Csf1r','Cd68','Mrc1']);cd8_support=n_support(['Cd8a','Cd8b1'])
labels=meta.broad_lineage.to_numpy();flags={
'T':np.isin(labels,frozen['T_cluster_categories'])&(t_support>=2)&(my_support<2)&(ep_support<2)&(en_support<2),
'Macrophage':np.isin(labels,frozen['Macro_cluster_categories'])&(macro_support>=2)&(t_support<2)&(ep_support<2)&(en_support<2),
'Fibroblast':np.isin(labels,frozen['Fibro_cluster_categories'])&(fib_support>=2)&(my_support<2)&(t_support<2)&(ep_support<2)&(en_support<2)&(mural_support<2),
'Mural':np.isin(labels,frozen['Mural_cluster_categories'])&(mural_support>=2)&(my_support<2)&(t_support<2)&(ep_support<2)&(en_support<2),
'Mesothelial_like':np.isin(labels,frozen['Meso_cluster_categories'])&(n_support(['Msln','Upk3b','Wt1','Krt19'])>=2)&(my_support<2)&(t_support<2)&(en_support<2),
'Epithelial_like':np.isin(labels,frozen['Epi_cluster_categories'])&(ep_support>=2)&(my_support<2)&(t_support<2)&(en_support<2)}
flags['CD8']=flags['T']&(cd8_support>=1);flags['Immune']=np.isin(labels,frozen['Immune_cluster_categories'])&(n_support(['Ptprc','Cd3d','Cd3e','Cd79a','Cd79b','Tyrobp','Lst1','Fcer1g','Nkg7','Jchain','Mzb1','Cpa3'])>=1)
flags['Immune']|=flags['T']|flags['Macrophage'];assert np.all(~flags['CD8']|flags['T'])
for k,v in flags.items():meta['lineage_'+k]=v
meta.to_csv(P/'tables/QUALIFIED_CELL_METADATA.tsv.gz',sep='\t')
human=json.loads((R/'phase292_barrier_execution_and_external_feasibility_20260930_v1/PANELS_FROZEN.json').read_text());orth=pd.read_csv(R/'phase307_spatial_perturbation_and_RCC_protein_20261002_v1/tables/MGI_ONE_TO_ONE_ORTHOLOGS.tsv',sep='\t');mapping=orth.set_index('human_gene').mouse_gene.to_dict();stromal={k:[mapping[g] for g in gs if g in mapping] for k,gs in human.items()};assert all(len(stromal[k])==len(gs) for k,gs in human.items())
panels=stromal|{'myCAF_markers':['Postn','Acta2','Lrrc15'],'steady_state_markers':['Pi16','Dpt'],'CD8_cytotoxicity':['Gzmb','Prf1','Nkg7','Ctsw'],'CD8_dysfunction':['Pdcd1','Ctla4','Lag3','Havcr2','Tigit','Tox','Cxcl13','Entpd1'],'CXCL9_10':['Cxcl9','Cxcl10']}
protocol={'created':datetime.datetime.now().astimezone().isoformat(),'panels':panels,'stromal_translation':'Phase292humanprograms,allone-to-oneMGImapped; nohumanGNLYorthologinmouseCD8panel','perturbation_zero_rule':'Measured zeros retained for all predefined readouts. Perpool eligibility uses gene measurement and lineage cellcount only. The inherited observational >=5detectedcells gate would censor real treatment-induced loss/gain and is replaced BEFORE effects by descriptive expression_supported/coverage fields. Sparse/allzero endpoints are not evidence of strong effects; noP/FDR. Initial protocol and pre-effect amendment retained.','statistics':'3unpairedpools/arm,nopooleffectP/noFDR/nocellP. MedianKO-minus-medianWT andall9cross-pooldirections descriptive only; missing qualification doesnotmeanzerobiology.'};(P/'resources/ENDPOINT_PROGRAMS_FROZEN.json').write_text(json.dumps(protocol,indent=2))
count_rows=[];measures=[];coverage=[]
for lib,m in meta.groupby('library',observed=True):
 ii=meta.library.eq(lib).to_numpy();group=str(m.group.iloc[0]);base={'library':str(lib),'group':group,'pool_tumors':2};n=int(ii.sum());rec=base|{'QC_cells':n}
 for lin,mask in flags.items():rec[lin]=int((ii&mask).sum())
 count_rows.append(rec)
 for endpoint,num,den in [('CD8_fraction_all','CD8','QC_cells'),('CD8_fraction_immune','CD8','Immune'),('T_fraction_all','T','QC_cells'),('Fibroblast_fraction_all','Fibroblast','QC_cells'),('Macrophage_fraction_all','Macrophage','QC_cells')]:
  v=rec[num]/rec[den] if rec[den]>0 else np.nan;measures.append(base|{'lineage':'Composition','endpoint':endpoint,'cells':rec[num],'denominator_cells':rec[den],'value':v,'eligible':rec[den]>=100,'reason':'' if rec[den]>=100 else 'DENOMINATOR_GATE','scale':'fraction'})
 for lin in ['Macrophage','Fibroblast','Mural','Mesothelial_like','Epithelial_like','CD8']:
  sel=ii&flags[lin];nc=int(sel.sum());sums=np.asarray(x[sel].sum(0)).ravel();tot=float(sums.sum());det=np.asarray((x[sel]>0).sum(0)).ravel();pb=np.log1p(sums*1e4/max(tot,1));cellgate=100 if lin in ['Macrophage','CD8','Epithelial_like'] else 20
  use={'Spp1':['Spp1'],'Tymp':['Tymp']}
  if lin in ['Fibroblast','Mural','Mesothelial_like']:use|={k:v for k,v in panels.items() if k in [*stromal,'myCAF_markers','steady_state_markers']}
  if lin=='Macrophage':use['CXCL9_10']=panels['CXCL9_10']
  if lin=='CD8':use|={k:panels[k] for k in ['CD8_cytotoxicity','CD8_dysfunction']}
  for ep,gs in use.items():
   keep=[g for g in gs if g in names];jj=names.get_indexer(keep);value=float(pb[jj].mean()) if len(jj) else np.nan;reason=[]
   if nc<cellgate:reason.append('CELL_GATE')
   if len(keep)<.8*len(gs):reason.append('MEASUREMENT_GATE')
   expression_supported=bool((det[jj]>=5).sum()>=.5*len(gs))
   measures.append(base|{'lineage':lin,'endpoint':ep,'cells':nc,'value':value,'eligible':not reason,'reason':';'.join(reason),'scale':'mean_logCP10K','measured_genes':len(keep),'expected_genes':len(gs),'supported_genes':int((det[jj]>=5).sum()),'expression_supported':expression_supported})
   for gene,j in zip(keep,jj):coverage.append(base|{'lineage':lin,'endpoint':ep,'gene':gene,'cells':nc,'detected_cells':int(det[j]),'gene_counts':float(sums[j]),'total_lineage_counts':tot,'logCP10K':float(pb[j])})
m=pd.DataFrame(measures);m.to_csv(P/'tables/POOL_ENDPOINTS_ALL.tsv',sep='\t',index=False);pd.DataFrame(count_rows).to_csv(P/'tables/LINEAGE_COUNTS_BY_POOL.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(P/'tables/PROGRAM_GENE_COVERAGE.tsv',sep='\t',index=False);effects=[]
for (lin,ep),g in m.groupby(['lineage','endpoint']):
 ko=g[g.group.eq('KO')&g.eligible].value.to_numpy();wt=g[g.group.eq('WT')&g.eligible].value.to_numpy();complete=len(ko)==len(wt)==3;mk=float(np.median(ko)) if len(ko) else np.nan;mw=float(np.median(wt)) if len(wt) else np.nan;diff=mk-mw if complete else np.nan;comparisons=np.array([a-b for a,b in itertools.product(ko,wt)]) if complete else np.array([])
 effects.append({'lineage':lin,'endpoint':ep,'KO_qualified_pools':len(ko),'WT_qualified_pools':len(wt),'KO_median':mk,'WT_median':mw,'median_difference':diff,'ratio_of_medians':mk/mw if complete and mw!=0 and g.scale.iloc[0]=='fraction' else np.nan,'KO_range_min':min(ko) if len(ko) else np.nan,'KO_range_max':max(ko) if len(ko) else np.nan,'WT_range_min':min(wt) if len(wt) else np.nan,'WT_range_max':max(wt) if len(wt) else np.nan,'cross_pool_comparisons':len(comparisons),'positive_comparisons':int((comparisons>0).sum()),'negative_comparisons':int((comparisons<0).sum()),'scale':g.scale.iloc[0],'eligible':complete,'reason':'' if complete else 'NOT_ALL3POOLS_QUALIFIED_PER_ARM','P':np.nan,'FDR':np.nan})
e=pd.DataFrame(effects);e.to_csv(P/'tables/PERTURBATION_EFFECTS_ALL.tsv',sep='\t',index=False);(P/'qc/POOLED_ANALYSIS_COMPLETE.json').write_text(json.dumps({'status':'COMPLETE','time':datetime.datetime.now().astimezone().isoformat(),'pool_rows':len(m),'effect_rows':len(e),'identity_sha256':hashlib.sha256((P/'resources/LINEAGE_IDENTITY_FROZEN.json').read_bytes()).hexdigest()},indent=2));print('POOLED_EFFECTS',e[e.lineage.isin(['Composition','Fibroblast','CD8'])].to_string(index=False),flush=True)
