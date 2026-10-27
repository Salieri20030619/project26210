from pathlib import Path
import json,collections,hashlib,sys,platform
import numpy as np,pandas as pd,scipy
from scipy import sparse
P=Path(__file__).resolve().parents[1]
src=P/'downloads/pilot_tumor1_region6.csv.gz'
raw=json.loads((P/'qc/PILOT_INPUT_VALIDATION.json').read_text())
assert hashlib.md5(src.read_bytes()).hexdigest()==raw['md5']
genes=sorted(raw['features']);gmap={g:i for i,g in enumerate(genes)}
cell_map={};all_counts=collections.Counter();nuc_counts=collections.Counter();nuc_sum=collections.defaultdict(lambda:np.zeros(4))
audit=collections.Counter();gene_all=collections.Counter();gene_nuc=collections.Counter()
for df in pd.read_csv(src,chunksize=250000,dtype={'cell_id':str,'feature_name':str}):
 audit['all_transcripts']+=len(df)
 assigned=~df.cell_id.isin(['-1','UNASSIGNED','nan','None',''])
 control=df.feature_name.str.match(r'(?i)^(NegControl|NegativeControl|Blank|Unassigned|antisense)')
 quality=df.qv.ge(20)
 audit['unassigned_transcripts']+=int((~assigned).sum());audit['qv_lt20_transcripts']+=int((~quality).sum());audit['control_transcripts']+=int(control.sum())
 keep=df[assigned & quality & ~control].copy();audit['qualified_assigned_transcripts']+=len(keep)
 assert np.isfinite(keep[['x_location','y_location','z_location']].to_numpy()).all()
 for c in keep.cell_id.unique():
  if c not in cell_map:cell_map[c]=len(cell_map)
 cnt=keep.groupby(['cell_id','feature_name'],sort=False).size();all_counts.update(cnt.to_dict())
 nu=keep[keep.overlaps_nucleus.eq(1)];audit['qualified_nuclear_transcripts']+=len(nu);nuc_counts.update(nu.groupby(['cell_id','feature_name'],sort=False).size().to_dict())
 gene_all.update(keep.feature_name.value_counts().to_dict());gene_nuc.update(nu.feature_name.value_counts().to_dict())
 sm=nu.groupby('cell_id')[['x_location','y_location','z_location']].sum();nc=nu.groupby('cell_id').size()
 for c,row in sm.iterrows():nuc_sum[c]+=np.array([row.x_location,row.y_location,row.z_location,nc[c]])
def matrix(counter):
 rows=np.fromiter((cell_map[c] for c,g in counter),dtype=np.int64)
 cols=np.fromiter((gmap[g] for c,g in counter),dtype=np.int64)
 vals=np.fromiter(counter.values(),dtype=np.int64)
 return sparse.coo_matrix((vals,(rows,cols)),shape=(len(cell_map),len(genes))).tocsr()
A=matrix(all_counts);N=matrix(nuc_counts)
sparse.save_npz(P/'resources/pilot_assigned_counts_qv20.npz',A);sparse.save_npz(P/'resources/pilot_nuclear_counts_qv20.npz',N)
rows=[]
for c,i in cell_map.items():
 sums=nuc_sum[c];xyz=sums[:3]/sums[3] if sums[3]>0 else np.repeat(np.nan,3)
 rows.append({'cell_id':c,'matrix_row':i,'assigned_RNA':int(A[i].sum()),'assigned_genes':int(A[i].getnnz()),'nuclear_RNA':int(N[i].sum()),'nuclear_genes':int(N[i].getnnz()),'nucleus_x_proxy_um':xyz[0],'nucleus_y_proxy_um':xyz[1],'nucleus_z_proxy_um':xyz[2]})
cells=pd.DataFrame(rows);cells['technical_qc_pass']=(cells.assigned_RNA>=20)&(cells.assigned_genes>=5)&(cells.nuclear_RNA>=5);cells.to_csv(P/'tables/PILOT_CELL_METADATA.tsv',sep='\t',index=False)
gf=pd.DataFrame({'gene':genes,'matrix_col':np.arange(len(genes)),'all_transcripts_raw':[raw['features'][g] for g in genes],'assigned_qv20_transcripts':[gene_all[g] for g in genes],'nuclear_qv20_transcripts':[gene_nuc[g] for g in genes],'assigned_cells':np.asarray((A>0).sum(axis=0)).ravel(),'nuclear_cells':np.asarray((N>0).sum(axis=0)).ravel()});gf.to_csv(P/'tables/PILOT_FEATURE_AUDIT.tsv',sep='\t',index=False)
sets={'project_targets':['TYMP','SPP1'],'T_CD8':['CD3D','CD3E','CD3G','TRAC','CD8A','CD8B'],'myeloid':['LST1','TYROBP','CD68','C1QA','C1QB','C1QC','FCGR3A','CSF1R'],'fibroblast':['COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRA','FAP'],'mural':['RGS5','CSPG4','MCAM','MYH11','PDGFRB'],'endothelial':['PECAM1','VWF','EMCN','KDR'],'epithelial_RCC':['EPCAM','KRT8','KRT18','KRT19','CA9','PAX8']}
markers=[]
for family,gs in sets.items():
 for g in gs:
  markers.append({'family':family,'gene':g,'observed_in_complete_file':g in gmap,'assigned_cells':int((A[:,gmap[g]]>0).sum()) if g in gmap else None,'nuclear_cells':int((N[:,gmap[g]]>0).sum()) if g in gmap else None,'status':'observed feature' if g in gmap else 'not observed; official panel coverage unresolved'})
pd.DataFrame(markers).to_csv(P/'tables/PILOT_MARKER_COVERAGE.tsv',sep='\t',index=False)
assert int(A.sum())==audit['qualified_assigned_transcripts'];assert int(N.sum())==audit['qualified_nuclear_transcripts'];assert (A-N).min()>=0
summary={'status':'RECONSTRUCTED_PENDING_INDEPENDENT_VALIDATION','input':str(src),'region':'tumor1_region6','donor_mapping':'unresolved','rows_raw':audit['all_transcripts'],'cells_assigned':len(cells),'technical_qc_cells':int(cells.technical_qc_pass.sum()),'cells_no_nuclear_coords':int(cells.nuclear_RNA.eq(0).sum()),'features_observed':len(genes),'counts':dict(audit),'matrix_shape':list(A.shape),'targets':{g:next(x for x in markers if x['gene']==g) for g in ['TYMP','SPP1']},'no_effect_tests':True,'coordinates':'quality passing nuclear transcript centroid proxy; not true segmented nucleus outline or membrane contact','original_segmentation':'author platform nuclear expansion assignment; not Proseg output','software':{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__}}
(P/'qc/RECONSTRUCTION_COMPLETE.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)
