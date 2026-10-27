from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='4'
import json,datetime,math,hashlib,itertools,subprocess
import numpy as np,pandas as pd,anndata as ad
P=Path(__file__).resolve().parents[1];checks=[]
def ck(name,ok):
 checks.append({'check':name,'passed':bool(ok)})
 if not ok:raise AssertionError(name)
def eq(a,b):return np.allclose(a,b,atol=1e-10,equal_nan=True)
for stage in ['prepare_inputs','qc_cluster','refine_identity','analyze_pools']:ck(stage+'_exit',(P/f'logs/{stage}.exit').read_text().strip()=='0')
plan=pd.read_csv(P/'tables/SAMPLE_PLAN.tsv',sep='\t');ck('three_pools_each',plan.groupby('group').size().to_dict()=={'KO':3,'WT':3});ck('two_tumors_per_pool_not_paired',plan.pool_tumors.eq(2).all() and (~plan.paired).all());q=pd.read_csv(P/'tables/QC_BY_LIBRARY.tsv',sep='\t');m=pd.read_csv(P/'tables/POOL_ENDPOINTS_ALL.tsv',sep='\t');e=pd.read_csv(P/'tables/PERTURBATION_EFFECTS_ALL.tsv',sep='\t');cov=pd.read_csv(P/'tables/PROGRAM_GENE_COVERAGE.tsv',sep='\t');count=pd.read_csv(P/'tables/LINEAGE_COUNTS_BY_POOL.tsv',sep='\t');meta=pd.read_csv(P/'tables/QUALIFIED_CELL_METADATA.tsv.gz',sep='\t',index_col=0);frozen=json.loads((P/'resources/LINEAGE_IDENTITY_FROZEN.json').read_text());panels=json.loads((P/'resources/ENDPOINT_PROGRAMS_FROZEN.json').read_text())['panels'];allgene=sorted(set(cov.gene));gi={g:i for i,g in enumerate(allgene)};acc={};scount={}
for row in plan.itertuples():
 subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True,stdout=subprocess.DEVNULL)
 a=ad.read_h5ad(P/f'inputs/{row.library}_counts_qc.h5ad');mm=meta.loc[a.obs_names];ck(row.library+'_metadata_alignment',mm.library.eq(row.library).all() and mm.group.eq(row.group).all());x=a.X.tocsr();ck(row.library+'_raw_integer',np.isfinite(x.data).all() and (x.data>=0).all() and (x.data==np.floor(x.data)).all());tot=np.asarray(x.sum(1)).ravel();ng=np.asarray((x>0).sum(1)).ravel();ck(row.library+'_QC',np.all(tot>=500) and np.all(ng>=200) and np.all(ng<=6000) and mm.pct_mito.le(20).all() and (~mm.predicted_doublet).all());ck(row.library+'_QC_count',len(mm)==q.loc[q.library.eq(row.library),'final_cells'].iloc[0]);scount[row.library]={'QC_cells':len(mm)}
 for lin in ['T','CD8','Immune','Macrophage','Fibroblast','Mural','Mesothelial_like','Epithelial_like']:
  sel=mm['lineage_'+lin].to_numpy(bool);scount[row.library][lin]=int(sel.sum());ck(row.library+lin+'_cluster_categories',set(mm.loc[sel,'broad_lineage'])<=set(frozen['T_cluster_categories'] if lin in ['T','CD8'] else frozen['Macro_cluster_categories'] if lin=='Macrophage' else frozen['Fibro_cluster_categories'] if lin=='Fibroblast' else frozen['Mural_cluster_categories'] if lin=='Mural' else frozen['Meso_cluster_categories'] if lin=='Mesothelial_like' else frozen['Epi_cluster_categories'] if lin=='Epithelial_like' else frozen['Immune_cluster_categories']))
  if lin not in ['Macrophage','Fibroblast','Mural','Mesothelial_like','Epithelial_like','CD8']:continue
  idx=a.var_names.get_indexer(allgene);ck(row.library+lin+'_genes_measured',(idx>=0).all());xx=x[sel];sums=np.asarray(xx[:,idx].sum(0)).ravel();det=np.asarray((xx[:,idx]>0).sum(0)).ravel();acc[(row.library,lin)]={'cells':int(sel.sum()),'total':float(xx.sum()),'sum':sums,'det':det}
 print('ORIGINAL_LIBRARY_VERIFIED',row.library,flush=True)
ck('allQC_retained_in_metadata',len(meta)==q.final_cells.sum()==count.QC_cells.sum())
for r in count.itertuples():
 for lin,v in scount[r.library].items():ck(r.library+lin+'_count_table',v==getattr(r,lin))
for i,r in cov.iterrows():
 a=acc[(r.library,r.lineage)];j=gi[r.gene];ck(f'gene{i}_count',a['sum'][j]==r.gene_counts and a['det'][j]==r.detected_cells and a['total']==r.total_lineage_counts);calc=math.log1p(a['sum'][j]*1e4/max(a['total'],1));ck(f'gene{i}_normalization',eq(calc,r.logCP10K))
for i,r in m.iterrows():
 if r.lineage=='Composition':
  num,den={'CD8_fraction_all':('CD8','QC_cells'),'CD8_fraction_immune':('CD8','Immune'),'T_fraction_all':('T','QC_cells'),'Fibroblast_fraction_all':('Fibroblast','QC_cells'),'Macrophage_fraction_all':('Macrophage','QC_cells')}[r.endpoint];q0=scount[r.library];ck(f'measure{i}_composition',eq(q0[num]/q0[den] if q0[den] else np.nan,r.value) and r.cells==q0[num] and r.denominator_cells==q0[den]);continue
 a=acc[(r.library,r.lineage)];gs=[r.endpoint] if r.endpoint in ['Spp1','Tymp'] else panels[r.endpoint];keep=[g for g in gs if g in gi];vv=[math.log1p(a['sum'][gi[g]]*1e4/max(a['total'],1)) for g in keep];ck(f'measure{i}_raw_pseudobulk',eq(math.fsum(vv)/len(keep) if keep else np.nan,r.value));minimum=100 if r.lineage in ['Macrophage','CD8','Epithelial_like'] else 20;ck(f'measure{i}_qualification',r.eligible==(a['cells']>=minimum and r.measured_genes>=.8*r.expected_genes));ck(f'measure{i}_no_detection_exclusion',pd.isna(r.reason) or 'DETECTION' not in r.reason)
for i,r in e.iterrows():
 g=m[m.lineage.eq(r.lineage)&m.endpoint.eq(r.endpoint)&m.eligible];ko=g[g.group.eq('KO')].value.tolist();wt=g[g.group.eq('WT')].value.tolist();ok=len(ko)==len(wt)==3;ck(f'effect{i}_N',r.KO_qualified_pools==len(ko) and r.WT_qualified_pools==len(wt) and r.eligible==ok)
 if ok:
  mk=sorted(ko)[1];mw=sorted(wt)[1];ck(f'effect{i}_medians',eq([mk,mw,mk-mw],[r.KO_median,r.WT_median,r.median_difference]));pairs=[a-b for a,b in itertools.product(ko,wt)];ck(f'effect{i}_directions',r.cross_pool_comparisons==9 and sum(v>0 for v in pairs)==r.positive_comparisons and sum(v<0 for v in pairs)==r.negative_comparisons)
 else:ck(f'effect{i}_not_imputed',pd.isna(r.median_difference) and r.cross_pool_comparisons==0)
ck('no_inferential_statistics',e.P.isna().all() and e.FDR.isna().all());ck('complete_pool_rows',len(m)==len(e)*6);ck('measured_zeroes_retained',m.loc[m.endpoint.eq('Spp1')&m.cells.ge(100),'eligible'].all());ck('no_CAF_relabeling',set(meta.loc[meta.lineage_Fibroblast,'broad_lineage'])=={'Fibroblast'})
(P/'qc/VALIDATION_CHECKS.json').write_text(json.dumps(checks,indent=2));(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'time':datetime.datetime.now().astimezone().isoformat(),'method':'Original per-library raw RNA reaggregation, exact pooled sampling units, independent math normalization and sorted median/cross-pool comparisons; no cell-level inference'},indent=2));print('VALIDATION_PASS',len(checks),flush=True)
