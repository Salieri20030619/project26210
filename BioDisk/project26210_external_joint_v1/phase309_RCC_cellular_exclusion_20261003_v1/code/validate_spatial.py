from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='4'
import hashlib,json,datetime
import numpy as np,pandas as pd,h5py
from scipy import sparse
from scipy.spatial.distance import cdist
from scipy.stats import rankdata
from sklearn.neighbors import BallTree
P=Path(__file__).resolve().parents[1];checks=[]
def check(name,condition,detail=''):
 checks.append({'name':name,'pass':bool(condition),'detail':str(detail)});assert condition,(name,detail)
record=json.loads((P/'qc/TENX_DOWNLOAD_VALIDATION.json').read_text())
for f in record['files']:
 path=P/'downloads'/f['name'];check('input_size_'+f['name'],path.stat().st_size==f['bytes']);check('input_SHA_'+f['name'],hashlib.sha256(path.read_bytes()).hexdigest()==f['sha256'])
for stage in ['acquire','identity_audit','analyze_spatial']:check('exit_'+stage,(P/'logs'/f'{stage}.exit').read_text().strip()=='0')
with h5py.File(P/'downloads/tenx_cell_feature_matrix.h5') as f:
 g=f['matrix'];x=sparse.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=g['shape'][:]).tocsr();names=np.array([z.decode() for z in g['features/name'][:]]);types=np.array([z.decode() for z in g['features/feature_type'][:]]);bc=np.array([z.decode() for z in g['barcodes'][:]])
rna=types=='Gene Expression';raw=pd.read_csv(P/'downloads/tenx_cells.csv.gz').set_index('cell_id').loc[bc];labels=pd.read_csv(P/'downloads/tenx_cell_groups.csv').set_index('cell_id')['group'].reindex(bc).fillna('Unannotated')
stored=pd.read_csv(P/'tables/CELL_IDENTITY_CANDIDATES.tsv.gz',sep='\t').set_index('cell_id').loc[bc];proto=json.loads((P/'resources/IDENTITY_FROZEN.json').read_text());panels=proto['panels'];support={}
for name,genes in panels.items():
 v=np.zeros(len(bc),dtype=np.int64)
 for gene in genes:
  jj=np.flatnonzero(rna&(names==gene))
  if len(jj):v+=x[jj[0]].toarray().ravel()>0
 support[name]=v;check('raw_marker_support_'+name,np.array_equal(v,stored[name+'_marker_support']))
total=np.asarray(x[rna].sum(axis=0)).ravel();xy=raw[['x_centroid','y_centroid']].to_numpy();q=(total>=10)&np.isfinite(xy).all(axis=1)&(raw.cell_area.to_numpy()>0)&raw.nucleus_count.eq(1).to_numpy()
check('rawRNA_totals',np.array_equal(total,stored.RNA_counts));check('source_csv_RNA_totals',np.array_equal(total,raw.transcript_counts));check('raw_coordinates',np.allclose(xy,stored[['x_centroid','y_centroid']],rtol=0,atol=1e-10));check('QC',np.array_equal(q,stored.QC))
cd8=q&labels.isin(proto['T_author_groups']).to_numpy()&(support['CD8']>=1)&(support['T']>=1)&(support['Myeloid']<2)&(support['Epithelial']<2)&(support['Endothelial']<2)&(support['B']==0)
fib=q&labels.eq('Unannotated').to_numpy()&(support['Fibroblast']>=2)&(support['T']<2)&(support['Myeloid']<2)&(support['Epithelial']<2)&(support['Endothelial']<2)&(support['Mural']<2)&(support['B']==0)
tum=q&labels.isin(proto['tumor_reference']).to_numpy()&(support['T']<2)&(support['Myeloid']<2)&(support['Endothelial']<2)&(support['B']==0)
ca9=tum&(x[np.flatnonzero(rna&(names=='CA9'))[0]].toarray().ravel()>0)
for name,mask in [('CD8',cd8),('Fibroblast_candidate',fib),('Tumor',tum),('Tumor_CA9',ca9)]:check('raw_identity_'+name,np.array_equal(mask,stored[name]))
check('disjoint_lineages',not np.any(cd8&fib|cd8&tum|fib&tum));check('TYMP_SPP1_absent',not np.any(rna&np.isin(names,['TYMP','SPP1'])))
endo=q&(support['Endothelial']>=2)&(support['T']<2)&(support['Myeloid']<2)&(support['Epithelial']<2)
mural=q&(support['Mural']>=2)&(support['Endothelial']<2)&(support['T']<2)&(support['Myeloid']<2)&(support['Epithelial']<2)
d=pd.read_csv(P/'tables/CD8_SPATIAL_READOUTS.tsv.gz',sep='\t').set_index('cell_id');check('CD8_ID_order',np.array_equal(d.index,bc[cd8]));pos=xy[cd8]
subset=np.unique(np.linspace(0,len(pos)-1,64).astype(int))
for name,mask in [('QC',q),('Fibroblast',fib),('Tumor',tum),('Tumor_CA9',ca9),('Endothelial',endo),('Mural',mural)]:
 tree=BallTree(xy[mask]);dist=cdist(pos[subset],xy[mask])
 for radius in ([50,100] if name=='Fibroblast' else [100]):
  col=f'{name}_count_{radius}um';v=tree.query_radius(pos,r=radius,count_only=True)-(1 if name=='QC' else 0);check('BallTree_all_'+col,np.array_equal(v,d[col]));v2=(dist<=radius).sum(axis=1)-(1 if name=='QC' else 0);check('bruteforce64_'+col,np.array_equal(v2,d[col].to_numpy()[subset]))
 if name in ['Tumor','Tumor_CA9']:
  v=tree.query(pos,k=1)[0].ravel();check('BallTree_all_distance_'+name,np.allclose(v,d[name+'_distance_um'],rtol=0,atol=1e-8));check('bruteforce64_distance_'+name,np.allclose(dist.min(axis=1),d[name+'_distance_um'].to_numpy()[subset],rtol=0,atol=1e-8));nearest=raw.loc[d[name+'_nearest_cell'],['x_centroid','y_centroid']].to_numpy();check('nearest_actual_coordinates_'+name,np.allclose(np.linalg.norm(pos-nearest,axis=1),v,rtol=0,atol=1e-8))
 print('RAW_GEOMETRY_VERIFIED',name,flush=True)
tile=np.floor((pos-xy[q].min(axis=0))/500).astype(int);check('fixed_grid',np.array_equal(tile,d[['tile_x','tile_y']].to_numpy()))
def independent_rho(g,exposure,target,controls):
 xx=rankdata(g[exposure]);yy=rankdata(g[target]);cols=[np.ones(len(g))]
 for c in controls:
  z=rankdata(g[c]);
  if z.std()>0:cols.append((z-z.mean())/z.std())
 design=np.column_stack(cols);rank=np.linalg.matrix_rank(design)
 if len(g)-rank<10 or rank!=design.shape[1] or np.ptp(xx)==0 or np.ptp(yy)==0:return np.nan,rank
 Q=np.linalg.qr(design,mode='reduced')[0];xx=xx-Q@(Q.T@xx);yy=yy-Q@(Q.T@yy)
 if np.linalg.norm(xx)<1e-9 or np.linalg.norm(yy)<1e-9:return np.nan,rank
 return float(np.dot(xx/np.linalg.norm(xx),yy/np.linalg.norm(yy))),rank
models=pd.read_csv(P/'tables/SPATIAL_ASSOCIATIONS_ALL.tsv',sep='\t');blocks=pd.read_csv(P/'tables/SPATIAL_BLOCK_EFFECTS_ALL.tsv',sep='\t');check('eight_models',len(models)==8);check('no_inferential_P',models[['P','FDR']].isna().all().all() and blocks[['P','FDR']].isna().all().all());check('one_donor',models.donors.eq(1).all())
for k,row in models.iterrows():
 g=d[d[row.target+'_distance_um']<=200];controls=[] if row.model=='unadjusted' else [row.target+'_count_100um','QC_count_100um','Endothelial_count_100um','Mural_count_100um','RNA_counts','x_centroid','y_centroid'];rho,rank=independent_rho(g,f'Fibroblast_count_{row.fibro_radius_um}um',row.target+'_distance_um',controls)
 check('model_'+str(k),np.isclose(rho,row.rho,atol=1e-10,equal_nan=True));check('model_n_'+str(k),len(g)==row.CD8_cells);check('model_rank_'+str(k),rank==row.design_rank)
for k,row in blocks.iterrows():
 g=d[(d[row.target+'_distance_um']<=200)&d.tile_500um.eq(row.tile)]
 controls=[row.target+'_count_100um','QC_count_100um','Endothelial_count_100um','Mural_count_100um','RNA_counts','x_centroid','y_centroid']
 rho=np.nan if len(g)<30 else independent_rho(g,f'Fibroblast_count_{row.fibro_radius_um}um',row.target+'_distance_um',controls)[0]
 check('block_'+str(k),np.isclose(rho,row.rho,atol=1e-10,equal_nan=True));check('block_n_'+str(k),len(g)==row.CD8_cells)
summ=pd.read_csv(P/'tables/SPATIAL_BLOCK_SUMMARY.tsv',sep='\t')
for k,row in summ.iterrows():
 v=blocks.loc[(blocks.target==row.target)&(blocks.fibro_radius_um==row.fibro_radius_um),'rho'].dropna();check('block_summary_'+str(k),len(v)==row.eligible_blocks and int((v>0).sum())==row.positive_blocks and int((v<0).sum())==row.negative_blocks and np.isclose(v.median(),row.median_rho,equal_nan=True))
contrasts=pd.read_csv(P/'tables/FIBRO_NEIGHBORHOOD_CONTRASTS.tsv',sep='\t')
for k,row in contrasts.iterrows():
 g=d[d[row.target+'_distance_um']<=200];x=g[f'Fibroblast_count_{row.fibro_radius_um}um'];lo,hi=np.quantile(x,[.25,.75]);g=g[(x<=lo) if row['group']=='low' else (x>=hi)] if hi>lo else g.iloc[:0]
 check('contrast_n_'+str(k),len(g)==row.CD8_cells);check('contrast_distance_'+str(k),np.isclose(g[row.target+'_distance_um'].median(),row.median_cancer_distance_um,equal_nan=True))
 for radius in [20,50]:check(f'contrast_fraction_{k}_{radius}',np.isclose((g[row.target+'_distance_um']<=radius).mean(),row[f'within{radius}um_fraction'],equal_nan=True))
completion=json.loads((P/'qc/SPATIAL_ANALYSIS_COMPLETE.json').read_text());check('protocol_hash_unchanged',hashlib.sha256((P/'resources/SPATIAL_PROTOCOL_FROZEN.json').read_bytes()).hexdigest()==completion['protocol_sha256']);check('identity_hash_unchanged',hashlib.sha256((P/'tables/CELL_IDENTITY_CANDIDATES.tsv.gz').read_bytes()).hexdigest()==completion['identity_sha256'])
(P/'qc/VALIDATION_CHECKS.json').write_text(json.dumps(checks,indent=2));result={'status':'PASS','checks':len(checks),'time':datetime.datetime.now().astimezone().isoformat(),'method':'Original H5 marker/QC identities; full independent BallTree geometry plus64 brute-force cdist cells; QR residual effects; all blocks/group contrasts/sample units; source SHA256'};(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(result,indent=2));print('VALIDATION_PASS',len(checks),flush=True)
