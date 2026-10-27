from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='4'
import h5py,json,datetime,hashlib
import numpy as np,pandas as pd
from scipy.sparse import csc_matrix
P=Path(__file__).resolve().parents[1]
with h5py.File(P/'downloads/tenx_cell_feature_matrix.h5') as f:
 g=f['matrix'];x=csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=g['shape'][:]).T.tocsr();names=np.array([v.decode() for v in g['features/name'][:]]);types=np.array([v.decode() for v in g['features/feature_type'][:]]);barcodes=np.array([v.decode() for v in g['barcodes'][:]])
assert len(np.unique(barcodes))==len(barcodes) and (x.data>=0).all()
cells=pd.read_csv(P/'downloads/tenx_cells.csv.gz').set_index('cell_id');groups=pd.read_csv(P/'downloads/tenx_cell_groups.csv').set_index('cell_id');assert cells.index.is_unique and groups.index.is_unique
assert set(barcodes)==set(cells.index) and set(groups.index)<=set(cells.index)
cells=cells.loc[barcodes].copy();cells['author_group']=groups['group'].reindex(cells.index).fillna('Unannotated')
rna=types=='Gene Expression';protein=types=='Protein Expression';print('FEATURE_TYPES',pd.Series(types).value_counts().to_dict(),flush=True)
def support(gs):
 jj=np.flatnonzero(rna&np.isin(names,gs));return np.asarray((x[:,jj]>0).sum(1)).ravel()
panels={'T':['CD3D','CD3E','TRAC'],'CD8':['CD8A','CD8B'],'Fibroblast':['DPT','COL5A2','PDGFRA','PCOLCE','MFAP5','OGN','SFRP2','THBS2','TNC','FBLN1'],'Mural':['MYH11','RERGL','CNN1','DES'],'Endothelial':['PECAM1','VWF','EMCN','KDR'],'Myeloid':['AIF1','CD68','MRC1','FCGR1A','MPEG1'],'B':['MS4A1','CD79A'],'Epithelial':['EPCAM','KRT7','CDH1','CA9'],'Target':['TYMP','SPP1','COL1A1','COL3A1','POSTN','FAP','ACTA2','CD44']}
tumor=['Hypoxic Tumor','Differentiated Proximal-Like Tumor','Growth Signaling-Enriched Tumor']
tgroups=['Exhausted Cytotoxic T Cells','Cycling Exhausted T Cells','Vascular-Localized Cytotoxic T Cells','Memory T Cells','Helper T Cells','Proliferating T Cells','Lymphatic-Associated T Cells','Regulatory T Cells (Tregs)']
proto={'time':datetime.datetime.now().astimezone().isoformat(),'sample_unit':'one donor, one tissue section; no P/FDR','panels':panels,'tumor_reference':tumor,'T_author_groups':tgroups,'qualification':'author-cell identity and gene marker eligibility reviewed before spatial effects','QC':'RNA total>=10; finite x/y; positive cell area; nucleus_count==1; RNA-only positive markers; no protein raw>0 positivity claims','CD8':'author T group AND >=1CD8 RNA AND >=1 T RNA marker; exclude >=2myeloid/epithelial/endothelial markers and >=1B marker','Fibroblast':'author Unannotated AND >=2 fixed fibro markers; exclude >=2T/myeloid/epithelial/endothelial/mural or >=1B; call fibroblast candidates, not all proven CAF','Tumor':'three explicitly author-labeled tumor groups; exclude >=2T/myeloid/endo or >=1B; group marker audit before freeze; sensitivity require CA9>=1','preserve_unknown':True,'excluded_genes':'TYMP/SPP1 not measured, no imputation','no_spatial_results_in_this_stage':True}
(P/'resources/IDENTITY_PROTOCOL_initial.json').write_text(json.dumps(proto,indent=2))
sup={k:support(gs) for k,gs in panels.items()};tot=np.asarray(x[:,rna].sum(1)).ravel();q=(tot>=10)&np.isfinite(cells[['x_centroid','y_centroid']]).all(axis=1).to_numpy()&(cells.cell_area.to_numpy()>0)&cells.nucleus_count.eq(1).to_numpy()
cd8=q&cells.author_group.isin(tgroups).to_numpy()&(sup['CD8']>=1)&(sup['T']>=1)&(sup['Myeloid']<2)&(sup['Epithelial']<2)&(sup['Endothelial']<2)&(sup['B']<1)
fib=q&cells.author_group.eq('Unannotated').to_numpy()&(sup['Fibroblast']>=2)&(sup['T']<2)&(sup['Myeloid']<2)&(sup['Epithelial']<2)&(sup['Endothelial']<2)&(sup['Mural']<2)&(sup['B']<1)
tum=q&cells.author_group.isin(tumor).to_numpy()&(sup['T']<2)&(sup['Myeloid']<2)&(sup['Endothelial']<2)&(sup['B']<1)
cells['QC']=q;cells['RNA_counts']=tot;cells['CD8']=cd8;cells['Fibroblast_candidate']=fib;cells['Tumor']=tum;cells['Tumor_CA9']=tum&(support(['CA9'])>=1)
for k,v in sup.items():cells[k+'_marker_support']=v
rows=[];coverage=[]
for group,mask in [('all_QC',q),('CD8',cd8),('Fibroblast_candidate',fib),('Tumor',tum),*[(g,q&cells.author_group.eq(g).to_numpy()) for g in tumor]]:
 for modality in ['Gene Expression','Protein Expression']:
  for gene in dict.fromkeys(sum(panels.values(),[])+['CD8A','CD3E','PanCK','CD68','CD31','alphaSMA']):
   jj=np.flatnonzero((types==modality)&(names==gene))
   if not len(jj):continue
   assert len(jj)==1
   v=x[mask,jj[0]].toarray().ravel();rows.append({'group':group,'modality':modality,'gene':gene,'cells':len(v),'detection_fraction':float(np.mean(v>0)) if len(v) else None,'mean_raw_counts':float(np.mean(v)) if len(v) else None,'median_raw_counts':float(np.median(v)) if len(v) else None})
for panel,gs in panels.items():
 for gene in gs:coverage.append({'panel':panel,'gene':gene,'RNA_measured':bool(np.any(rna&(names==gene))),'protein_measured':bool(np.any((types=='Protein Expression')&(names==gene)))})
pd.DataFrame(rows).to_csv(P/'tables/IDENTITY_MARKER_AUDIT.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(P/'tables/PANEL_COVERAGE.tsv',sep='\t',index=False)
cells.to_csv(P/'tables/CELL_IDENTITY_CANDIDATES.tsv.gz',sep='\t')
counts={c:int(cells[c].sum()) for c in ['QC','CD8','Fibroblast_candidate','Tumor','Tumor_CA9']};counts['all_input']=len(cells)
(P/'qc/IDENTITY_AUDIT_COMPLETE.json').write_text(json.dumps({'status':'READY_FOR_REVIEW','counts':counts,'RNA_features':int(rna.sum()),'protein_features':int((types=='Protein Expression').sum()),'time':datetime.datetime.now().astimezone().isoformat()},indent=2))
print('COUNTS',counts,flush=True)
z=pd.DataFrame(rows);print(z[(z.group.isin(['CD8','Fibroblast_candidate','Tumor']))&(z.modality=='Gene Expression')&(z.detection_fraction>.1)].to_string(index=False),flush=True)
