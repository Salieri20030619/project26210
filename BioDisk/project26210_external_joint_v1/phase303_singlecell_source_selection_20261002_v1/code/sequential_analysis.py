from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='4'
import json,datetime,hashlib,subprocess
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import rankdata,spearmanr
P=Path(__file__).resolve().parents[1];R=P.parent
pan=json.loads((P/'resources/QUALIFICATION_PROTOCOL.json').read_text())
func=pan['function_panels'];ctx=json.loads((R/'phase280_TYMP_source_niche_function_20260928_v1/resources/CONTEXT_PROGRAMS.json').read_text())['genes']
im=pd.read_csv(R/'phase213_immune_restart_execution_20260919_v1/I1_patient_inputs/PROGRAM_GENE_COVERAGE.tsv',sep='\t')
immune={k:im.loc[im.endpoint.eq(k),'gene'].unique().tolist() for k in ['CD8_cytotoxicity','CD8_dysfunction']}
ctx={k:[g for g in gs if g not in ['TYMP','SPP1']] for k,gs in ctx.items()}
rules={'created':datetime.datetime.now().astimezone().isoformat(),'selection':'Frozen lineage qualification; primary untreated GSE210038; GSE224630 untreated; GSE222703 therapy unreported. No cross-study patient pooling; treated GSE210038 individual records only. No outcome-based selection.','lineage_gates':{'Macrophage':100,'RCC_like':50,'Fibro_strict':20,'Fibro_transition':20,'Mural':50,'CD8':20},'single_gene_gate':'Measured in original library; >=5 expressing cells for patient aggregate; within-macrophage TYMP/SPP1 each >=10 expressing cells.','program_gate':'At least 80% panel genes measured in original library; at least 50% full-panel genes detected in >=5% cells of compartment.','programs':func|immune,'context':ctx,'summary':'Pseudobulk sum raw counts per patient-lineage; log1p(CP10K) per gene; equal-weight mean fixed panel; same-library pair only. Sparse missing library genes never treated as measured zero.','within_macrophage':'Patient-specific partial Spearman: residual ranks of logCP10K TYMP/SPP1 on ranks(logUMI,nGenes); environment adds ranks of mean logCP10K IFNalpha,IFNgamma,Hypoxia. No cell-level P.','across_patient':'Per study Spearman effect only if >=3 unique eligible patients and both measures vary; N<10 no P/FDR; no pooled inference in this source qualification round.','stages':['1_TYMP_SPP1','2_sources_stroma','3_immune_function'],'interpretation':'RCC_like denotes transcriptome-supported tumour-like identity without new CNV proof. Fibro_transition is marker-overlap perivascular ECM-producing pool, not pure CAF; expression is not protein action, physical matrix alignment, migration or exclusion.'}
fr=P/'resources/SEQUENTIAL_ANALYSIS_PROTOCOL.json'
if not fr.exists():fr.write_text(json.dumps(rules,ensure_ascii=False,indent=2))
else:raise RuntimeError('Protocol already exists; preserve previous run and explicitly version any rerun')
def gate():subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True,stdout=subprocess.DEVNULL)
def partial(x,y,cols):
 d=np.column_stack([np.ones(len(x))]+[rankdata(c) for c in cols if np.ptp(c)>0]);rx=rankdata(x);ry=rankdata(y)
 rx=rx-d@np.linalg.lstsq(d,rx,rcond=None)[0];ry=ry-d@np.linalg.lstsq(d,ry,rcond=None)[0]
 return float(np.corrcoef(rx,ry)[0,1]) if np.std(rx)>1e-9 and np.std(ry)>1e-9 else np.nan
rows=[];within=[];cov=[]
for ds in ['GSE210038','GSE224630','GSE222703']:
 gate();a=ad.read_h5ad(P/f'inputs/{ds}_counts_qc.h5ad');meta=pd.read_csv(P/f'tables/{ds}_QUALIFIED_CELL_METADATA.tsv.gz',sep='\t',index_col=0)
 assert a.obs_names.equals(meta.index);x=a.X.tocsr().astype(np.float64)
 norm=x.multiply((1e4/np.asarray(x.sum(axis=1)).ravel())[:,None]).tocsr();norm.data=np.log1p(norm.data)
 for sample in a.obs['sample'].unique():
  lib=ad.read_h5ad(P/f'inputs/{sample}_qc.h5ad',backed='r');present=set(lib.var_names);lib.file.close()
  ix=meta['sample'].eq(sample).to_numpy();base={'dataset':ds,'sample':sample,'patient':str(meta.loc[ix,'patient'].iloc[0]),'main':bool(meta.loc[ix,'main'].iloc[0]),'treatment':str(meta.loc[ix,'treatment'].iloc[0])}
  for lin,nmin in rules['lineage_gates'].items():
   sel=ix&meta['lineage_'+lin].to_numpy(dtype=bool);n=int(sel.sum());summ=np.asarray(x[sel].sum(axis=0)).ravel();pb=np.log1p(summ*1e4/max(summ.sum(),1));det=np.asarray((x[sel]>0).sum(axis=0)).ravel()
   measures={'TYMP':['TYMP'],'SPP1':['SPP1']} if lin in ['Macrophage','RCC_like'] else (immune if lin=='CD8' else func)
   for name,genes in measures.items():
    available=[g for g in genes if g in present];jj=a.var_names.get_indexer(available);assert (jj>=0).all()
    if len(genes)==1:expr=int(det[jj].sum());ok=(n>=nmin and len(available)==1 and expr>=5)
    else:expr=int((det[jj]>=max(1,.05*n)).sum());ok=n>=nmin and len(available)>=.8*len(genes) and expr>=.5*len(genes)
    val=float(pb[jj].mean()) if len(jj) else np.nan
    rows.append(base|{'lineage':lin,'measure':name,'n_cells':n,'minimum_cells':nmin,'genes_expected':len(genes),'genes_measured':len(available),'expressing_cells_or_genes':expr,'eligible':bool(ok),'value':val,'missing_genes':';'.join(sorted(set(genes)-present))})
    for g,j in zip(available,jj):cov.append(base|{'lineage':lin,'measure':name,'gene':g,'n_cells':n,'detected_cells':int(det[j]),'pseudobulk_logCP10K':float(pb[j])})
   if lin!='Macrophage' or n<nmin:continue
   vals={g:norm[sel,a.var_names.get_loc(g)].toarray().ravel() for g in ['TYMP','SPP1'] if g in present}
   if len(vals)<2 or any((v>0).sum()<10 for v in vals.values()):continue
   cols=[meta.loc[sel,'total_counts'].to_numpy(),meta.loc[sel,'n_genes'].to_numpy()]
   rr=base|{'n_cells':n,'TYMP_detected':int((vals['TYMP']>0).sum()),'SPP1_detected':int((vals['SPP1']>0).sum())}
   within.append(rr|{'model':'depth_complexity','rho':partial(vals['TYMP'],vals['SPP1'],cols),'P':np.nan,'FDR':np.nan})
   env=[]
   for k,gs in ctx.items():
    keep=[g for g in gs if g in present]
    if len(keep)<.8*len(gs):break
    env.append(np.asarray(norm[sel][:,a.var_names.isin(keep)].mean(axis=1)).ravel())
   if len(env)==len(ctx):within.append(rr|{'model':'depth_complexity_environment','rho':partial(vals['TYMP'],vals['SPP1'],cols+env),'P':np.nan,'FDR':np.nan})
 print('AGGREGATED',ds,flush=True)
df=pd.DataFrame(rows);df.to_csv(P/'tables/PATIENT_LINEAGE_MEASURES.tsv',sep='\t',index=False)
pd.DataFrame(cov).to_csv(P/'tables/FUNCTION_GENE_COVERAGE.tsv',sep='\t',index=False)
w=pd.DataFrame(within);w.to_csv(P/'tables/STEP1_WITHIN_MACRO_PATIENT.tsv',sep='\t',index=False)
summary=[]
for (ds,model),g in w[w.main].groupby(['dataset','model']):
 summary.append({'dataset':ds,'model':model,'n_patients':g.patient.nunique(),'positive':int((g.rho>0).sum()),'negative':int((g.rho<0).sum()),'median_rho':float(g.rho.median()),'P':np.nan,'FDR':np.nan,'status':'exploratory_descriptive_N_lt_10'})
pd.DataFrame(summary).to_csv(P/'tables/STEP1_WITHIN_MACRO_SUMMARY.tsv',sep='\t',index=False)
edges=[];pairs=[]
def edge(stage,src,tar,family):
 for ds in df.dataset.unique():
  d=df[df.dataset.eq(ds)&df.main&df.eligible];aa=d[d.lineage.eq(src[0])&d.measure.eq(src[1])][['sample','patient','value']];bb=d[d.lineage.eq(tar[0])&d.measure.eq(tar[1])][['sample','patient','value']]
  g=aa.merge(bb,on=['sample','patient'],suffixes=('_source','_target'),validate='one_to_one');n=g.patient.nunique();assert n==len(g)
  good=n>=3 and g.value_source.nunique()>1 and g.value_target.nunique()>1
  rho=float(spearmanr(g.value_source,g.value_target).statistic) if good else np.nan
  rec={'stage':stage,'dataset':ds,'source':':'.join(src),'target':':'.join(tar),'family':family,'n_patients':n,'rho':rho,'P':np.nan,'FDR':np.nan,'status':'exploratory_descriptive_N_lt_10' if good else 'not_evaluable_insufficient_patients_or_variation'}
  assert n<10,'Unexpected eligible N: freeze inferential model before calculation'
  edges.append(rec)
  for rr in g.to_dict('records'):pairs.append(rec|rr)
edge(1,('RCC_like','TYMP'),('Macrophage','SPP1'),'tumor_TYMP_macro_SPP1_primary')
pd.DataFrame([e for e in edges if e['stage']==1]).to_csv(P/'tables/STEP1_TUMOR_MACRO.tsv',sep='\t',index=False)
(P/'qc/STEP1_COMPLETE.json').write_text(json.dumps({'status':'COMPLETE','next':'stroma'}))
for lin in ['Fibro_strict','Fibro_transition','Mural']:
 for src in [('RCC_like','TYMP'),('Macrophage','TYMP'),('Macrophage','SPP1')]:
  for end in func:edge(2,src,(lin,end),'stroma_'+lin+'_'+'_'.join(src)+'_3endpoints')
pd.DataFrame([e for e in edges if e['stage']==2]).to_csv(P/'tables/STEP2_STROMA.tsv',sep='\t',index=False)
(P/'qc/STEP2_COMPLETE.json').write_text(json.dumps({'status':'COMPLETE','next':'immune function'}))
for src in [('RCC_like','TYMP'),('Macrophage','TYMP'),('Macrophage','SPP1')]+[(lin,'collagen_maturation') for lin in ['Fibro_strict','Fibro_transition','Mural']]:
 for end in immune:edge(3,src,('CD8',end),'immune_'+'_'.join(src)+'_2endpoints')
pd.DataFrame([e for e in edges if e['stage']==3]).to_csv(P/'tables/STEP3_IMMUNE.tsv',sep='\t',index=False)
pd.DataFrame(edges).to_csv(P/'tables/ALL_SEQUENTIAL_EDGES.tsv',sep='\t',index=False)
pd.DataFrame(pairs).to_csv(P/'tables/ALL_MATCHED_PATIENT_PAIRS.tsv',sep='\t',index=False)
(P/'qc/STEP3_COMPLETE.json').write_text(json.dumps({'status':'COMPLETE','next':'validation and report'}))
print('SEQUENTIAL_ANALYSIS_COMPLETE',len(edges),'edges',flush=True)
