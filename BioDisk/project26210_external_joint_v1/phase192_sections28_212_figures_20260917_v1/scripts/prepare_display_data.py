from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy.sparse import load_npz
from scipy.stats import rankdata
from scipy.cluster.hierarchy import linkage,leaves_list
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');q=r/'part2_rebuild_v1';d=r/'phase192_sections28_212_figures_20260917_v1'
for sub in ['inputs','tables','figures','qc','scripts']:(d/sub).mkdir(exist_ok=True)
plan={'authorization':'用户2026-09-17明确要求2.8到2.12各1-2代表图，拒绝散点/箱线等纯统计图；本批独立交付，不自动插报告。','design_reference':'sciadv.aed0733 user PDF pp3/5/7/10 inspected','figures':{'2.8':['01_MTA_metabolic_map'],'2.9':['02_macrophage_density_atlas','03_macrophage_marker_heatmap'],'2.10':['04_MTA_SPP1_tissue_landscape'],'2.11':['05_macrophage_program_heatmap','06_C4_expression_landscapes'],'2.12':['07_partitioned_trajectory_maps','08_pseudotime_expression_heatmap']},'formats':['PDF vector','SVG vector','PNG 600dpi'],'freeze':'No new clustering, trajectory fit, or inferential tests. Density smoothing is display only. New display aggregations do not replace published statistics.'}
(d/'PLAN_BEFORE_PLOTTING.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2))
m=pd.read_csv(q/'08_regulons/input_v1/CELL_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str},low_memory=False)
genes=pd.read_csv(q/'08_regulons/input_v1/GENES.tsv',sep='\t').gene.tolist();X=load_npz(q/'08_regulons/input_v1/GENES_BY_CELLS_RAW.npz').tocsr()
assert X.shape==(len(genes),len(m))
clusters=['0','3','4','5','6','7','8','9','10','11','12','15']
ma=m[m.candidate.isin(clusters)].copy();assert len(ma)==22444 and (ma.candidate=='4').sum()==2674
uv=pd.read_csv(q/'03_macrophage_clustering/run_v4_exploratory/UMAP_COORDINATES.tsv.gz',sep='\t',index_col=0).loc[ma.index]
ma[['UMAP1','UMAP2']]=uv[['UMAP1','UMAP2']]
sc=pd.read_csv(q/'07_trajectory/dynamics_v2/CELL_EXPRESSION_PROGRAM_SCORES.tsv.gz',sep='\t',index_col=0)
assert ma.index.isin(sc.index).all()
programs=['Cholesterol metabolism','Cholesterol Homeostasis','Fatty Acid Metabolism','Lysosome','Phagosome','Antigen processing and presentation','Inflammatory Response','Interferon Gamma Response','Hypoxia','Angiogenesis']
for name in programs:ma['score:'+name]=sc.loc[ma.index,'program:'+name]
marker_groups={'C0':['CCL3','CCL4','IL1B'],'C3':['CX3CR1','FCGR3A','LILRB1'],'C4':['GPNMB','CCL18','SPP1','CYP27A1','LPL','APOE'],'C5':['TNFSF10','GBP4','IFIT3'],'C6':['ICAM1','CXCL2','CXCL3'],'C7':['HSPA6','HSPA1A','HSPA1B'],'C8':['MRC1','FOLR2','LYVE1'],'C9':['HBEGF','VEGFA','ADM'],'C10':['USP53','ALOX5AP'],'C11':['BAG3','KLF2','DUSP1'],'C12':['MARCO','APOC1','C1QC'],'C15':['MKI67','TOP2A','STMN1']}
mg=[g for group in marker_groups.values() for g in group]
dynamic_groups={'Monocyte':['FCN1','VCAN','S100A8','S100A9','LYZ'],'C4 identity':['SPP1','GPNMB','CCL18'],'Lipid handling':['APOE','LPL','CYP27A1','ABCA1','ABCG1','NPC1','SOAT1'],'Lysosome':['CTSB','CTSD','CTSL','LAMP1','LAMP2','LIPA','ATP6V0D1'],'Antigen presentation':['HLA-DRA','HLA-DPA1','HLA-DPB1','CD74'],'Inflammation':['IL1B','TNF','CCL3','CCL4'],'Macrophage backbone':['C1QA','C1QB','C1QC','TYROBP','CSF1R']}
dg=[g for v in dynamic_groups.values() for g in v]
needed=sorted(set(mg+dg));missing=[g for g in needed if g not in genes]
assert not missing,missing
E=np.log1p(X[[genes.index(g) for g in needed]].toarray().T/m.total_counts_raw.to_numpy()[:,None]*1e4)
for g in ['SPP1','GPNMB','CCL18']:ma['expr:'+g]=E[m.index.get_indexer(ma.index),needed.index(g)]
ma.to_csv(d/'inputs/MACRO_DISPLAY_METADATA.tsv.gz',sep='\t',compression='gzip')
# Patient-equal expression matrices; each patient x candidate has >=20 cells.
pbroot=q/'04_markers_annotation/run_v4_macrophage_only'
pb=load_npz(pbroot/'PATIENT_CLUSTER_PSEUDOBULK_COUNTS.npz').tocsr();ri=pd.read_csv(pbroot/'PSEUDOBULK_ROW_INDEX.tsv',sep='\t',dtype={'cluster':str});pg=pd.read_csv(pbroot/'PSEUDOBULK_GENE_INDEX.tsv',sep='\t').gene.tolist()
assert len(ri)==pb.shape[0]
tot=np.asarray(pb.sum(axis=1)).ravel()
V=np.log2(pb[:,[pg.index(g) for g in mg]].toarray()/tot[:,None]*1e6+1)
mat=[];sizes=[]
for c in clusters:
 ix=(ri.cluster==c)&(ri.n_cells>=20);assert ix.sum()>0
 mat.append(V[ix].mean(0));sizes.append(int(ix.sum()))
hm=pd.DataFrame(np.array(mat).T,index=mg,columns=['C'+c for c in clusters]);hm.to_csv(d/'tables/MARKER_PATIENT_EQUAL_LOG2CPM.tsv',sep='\t')
pmeans=ma.groupby(['curated_patient_id','candidate'])[[('score:'+s) for s in programs]].mean()
ns=ma.groupby(['curated_patient_id','candidate']).size();pmeans=pmeans[ns>=20]
pmean=pmeans.groupby('candidate').mean().reindex(clusters).T;pmean.index=programs;pmean.columns=['C'+c for c in clusters];pmean.to_csv(d/'tables/PROGRAM_PATIENT_EQUAL_SCORES.tsv',sep='\t')
labels={'0':'CCL3 / CCL4','3':'CX3CR1 / FCGR3A','4':'GPNMB / CCL18 / SPP1','5':'TNFSF10 / GBP4','6':'ICAM1 / CCL3','7':'HSPA6 / BAG3','8':'MRC1 / FOLR2','9':'HBEGF / VEGFA','10':'USP53 / ALOX5AP','11':'BAG3 / KLF2','12':'MARCO / APOC1','15':'MKI67 / TOP2A'}
(d/'inputs/FEATURE_REGISTRY.json').write_text(json.dumps({'marker_groups':marker_groups,'dynamic_groups':dynamic_groups,'clusters':clusters,'labels':labels,'marker_patient_n':dict(zip(['C'+c for c in clusters],sizes))},indent=2))
# External paired tissue data: every patient kept, independently ordered by measured MTA within cohort.
cohorts=['RC18','RC20','CPTAC','CPTAC_val'];bulk=[]
for c in cohorts:
 f=pd.read_csv(r/'phase162_MTA_polyamine_part1_reassessment_20260915_v1'/f'PUBLIC_PATIENTS_{c}.tsv',sep='\t',index_col=0)
 assert f.index.is_unique
 f=f.sort_values('MET:5-methylthioadenosine (MTA)',kind='stable')
 f['cohort']=c;f['within_cohort_MTA_rank']=np.arange(len(f))+1;bulk.append(f)
bulk=pd.concat(bulk);assert len(bulk)==255
bulk.to_csv(d/'tables/TISSUE_MTA_ORDERED_PATIENTS.tsv',sep='\t')
# Dynamic display: exact pre-existing tumor SPP1 patient/partition eligibility.
tr=pd.read_csv(q/'07_trajectory/dynamics_v2/PARTITION_PATIENT_DYNAMIC_TRENDS.tsv.gz',sep='\t',dtype={'partition':str})
gate=tr[(tr.scope=='tumor_only')&(tr.pt_version=='pt_primary')&(tr.feature=='gene:SPP1')&np.isfinite(tr.depth_adjusted_rho)].copy()
pt=pd.read_csv(q/'07_trajectory/monocle3_v3/ALL_CELL_PSEUDOTIME.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str,'partition':str},low_memory=False)
grid=np.linspace(.025,.975,40);dynmats={};audit=[];comps={}
for mode in ['primary','without_stress_cycle']:
 records={};composition={}
 sub=pt[(pt['mode']==mode)&(pt.source_tissue.str.lower()=='tumor')]
 for z in gate[gate['mode']==mode].itertuples():
  ix=(sub.curated_patient_id==z.patient)&(sub.partition==str(z.partition))&(sub.study==z.study)&np.isfinite(sub.pt_primary)
  md=sub[ix];assert len(md)==z.cells,(mode,z.patient,z.partition,len(md),z.cells)
  t=(rankdata(md.pt_primary)-.5)/len(md);w=np.exp(-.5*((grid[:,None]-t[None,:])/.075)**2)
  w/=w.sum(axis=1,keepdims=True)
  vals=E[m.index.get_indexer(md.index)][:,[needed.index(g) for g in dg]]
  sm=w@vals
  records.setdefault(z.patient,[]).append(sm)
  cm=w@np.column_stack([(md.candidate=='4').to_numpy(float),md.candidate.isin(['1','Monocyte']).to_numpy(float)])
  composition.setdefault(z.patient,[]).append(cm)
  audit.append({'mode':mode,'patient':z.patient,'study':z.study,'partition':z.partition,'cells':len(md),'C4_cells':int((md.candidate=='4').sum())})
 dynmats[mode]=np.mean([np.mean(v,axis=0) for v in records.values()],axis=0).T
 comps[mode]=np.mean([np.mean(v,axis=0) for v in composition.values()],axis=0).T
 assert len(records)==(14 if mode=='primary' else 24),len(records)
 pd.DataFrame(dynmats[mode],index=dg,columns=np.round(grid,3)).to_csv(d/'tables'/f'DYNAMIC_{mode}_MEAN_LOGCP10K.tsv',sep='\t')
np.savez_compressed(d/'inputs/DYNAMIC_DISPLAY.npz',genes=np.array(dg),grid=grid,primary=dynmats['primary'],sensitivity=dynmats['without_stress_cycle'],composition_primary=comps['primary'],composition_sensitivity=comps['without_stress_cycle'])
pd.DataFrame(audit).to_csv(d/'tables/DYNAMIC_DISPLAY_ELIGIBILITY.tsv',sep='\t',index=False)
(d/'inputs/DISPLAY_STATUS.json').write_text(json.dumps({'mature_macrophages':len(ma),'patients':int(ma.curated_patient_id.nunique()),'C4':2674,'tissue_patients':255,'dynamic_patients':[14,24],'missing_markers':missing,'programs':programs,'new_inferential_tests':False},indent=2))
print('PREPARED',len(ma),hm.shape,pmean.shape,len(bulk),[(k,v.shape) for k,v in dynmats.items()],flush=True)
