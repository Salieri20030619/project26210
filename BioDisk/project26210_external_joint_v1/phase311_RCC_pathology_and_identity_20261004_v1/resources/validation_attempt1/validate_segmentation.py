from pathlib import Path
import json,hashlib,gzip,collections,datetime
import numpy as np,pandas as pd,anndata as ad
P=Path(__file__).resolve().parents[1];q=P/'resources/proseg_region6_v1';old=P.parent/'phase310_RCC_mapped_spatial_inputs_20261004_v1'
checks=[]
def ck(n,b):
 checks.append({'check':n,'pass':bool(b)})
 if not b:raise AssertionError(n)
ck('proseg_exit0',(P/'logs/proseg.exit').read_text().strip()=='0')
ck('cluster_exit0',(P/'logs/cluster_identity.exit').read_text().strip()=='0')
M=pd.read_csv(q/'maxpost-counts.csv.gz');E=pd.read_csv(q/'expected-counts.csv.gz');cm=pd.read_csv(q/'cell-metadata.csv.gz',dtype={'original_cell_id':str})
ck('counts_metadata_dimensions',len(M)==len(E)==len(cm))
ck('point_and_expected_gene_order',list(M)==list(E))
ck('integer_nonnegative_point_counts',np.isfinite(M.to_numpy()).all() and (M.to_numpy()>=0).all() and np.equal(M.to_numpy(),np.floor(M.to_numpy())).all())
ck('finite_nonnegative_expected_counts',np.isfinite(E.to_numpy()).all() and (E.to_numpy()>=0).all())
ck('cell_integer_row_index',np.array_equal(cm.cell,np.arange(len(cm))))
ck('unique_original_cell_ids',cm.original_cell_id.is_unique)
raw=pd.read_csv(old/'downloads/pilot_tumor1_region6.csv.gz',usecols=['transcript_id','feature_name','x_location','y_location','qv'])
ck('raw_unique_transcript_ids',raw.transcript_id.is_unique)
raw=raw.set_index('transcript_id');gmap={g:i for i,g in enumerate(M.columns)};recount=np.zeros(M.shape,dtype=np.int64)
seen=[];max_xy_error=0.;output_rows=0;accepted_rows=0;headers=None
for df in pd.read_csv(q/'transcript-metadata.csv.gz',chunksize=250000):
 if headers is None:
  headers=list(df);ck('transcript_metadata_schema14',len(headers)==14)
 # Schema order independently checked in official output.rs: id, xyz, original xyz, gene,qv,fov,cell,prob,bg,confusion.
 ids=df.iloc[:,0].to_numpy();ck('source_id_block_'+str(output_rows),pd.Index(ids).isin(raw.index).all())
 ref=raw.loc[ids];ck('gene_identity_block_'+str(output_rows),np.array_equal(ref.feature_name.astype(str),df.iloc[:,7].astype(str)))
 err=max(float(np.max(np.abs(ref.x_location.to_numpy()-df.iloc[:,4].to_numpy()))),float(np.max(np.abs(ref.y_location.to_numpy()-df.iloc[:,5].to_numpy()))));max_xy_error=max(max_xy_error,err)
 ck('QV20_block_'+str(output_rows),(df.iloc[:,8]>=20).all())
 ck('probability_bounds_block_'+str(output_rows),df.iloc[:,11].between(0,1).all())
 ci=df.iloc[:,10].to_numpy(np.uint64);pr=df.iloc[:,11].to_numpy(np.float32);good=(ci<len(cm))&(pr>np.float32(0.1))
 gi=df.iloc[:,7].map(gmap);ck('gene_dictionary_block_'+str(output_rows),gi.notna().all())
 np.add.at(recount,(ci[good].astype(int),gi.to_numpy(int)[good]),1)
 accepted_rows+=int(good.sum());output_rows+=len(df);seen.append(ids)
seen=np.concatenate(seen);ck('unique_output_transcript_ids',len(np.unique(seen))==len(seen))
ck('all_point_counts_recomputed',np.array_equal(recount,M.to_numpy()))
ck('all_point_counts_conserved',int(M.to_numpy().sum())==accepted_rows)
ck('original_xy_float32_tolerance',max_xy_error<0.002)
oldcells=pd.read_csv(old/'tables/PILOT_CELL_METADATA.tsv',sep='\t',dtype={'cell_id':str});ck('original_cell_ids_exist',set(cm.original_cell_id).issubset(set(oldcells.cell_id)))
a=ad.read_h5ad(P/'resources/proseg_region6_identity.h5ad');meta=pd.read_csv(P/'tables/PROSEG_CELL_IDENTITY_INPUT.tsv',sep='\t',dtype={'original_cell_id':str})
sel=meta.technical_qc.to_numpy(bool)
ck('QC_rule_independent',np.array_equal(sel,(M.sum(axis=1)>=20)&((M>0).sum(axis=1)>=5)&np.isfinite(cm[['centroid_x','centroid_y']].to_numpy()).all(axis=1)))
ck('h5ad_rows',np.array_equal(a.obs_names.astype(str),cm.loc[sel,'original_cell_id'].astype(str)))
ck('h5ad_counts',np.array_equal(a.layers['maxpost_counts'].toarray(),M.loc[sel].to_numpy()))
ck('h5ad_coordinates',np.allclose(a.obsm['spatial'],cm.loc[sel,['centroid_x','centroid_y']]))
ck('clustering_all_QC',a.obs.identity_cluster.notna().all())
ck('target_excluded_record',set(a.uns['identity_parameters']['target_excluded'])=={'TYMP','SPP1'})
manifest={}
for f in q.iterdir():
 if f.is_file():manifest[f.name]={'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()}
(P/'qc/PROSEG_OUTPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
rec={'status':'PASS','time':datetime.datetime.now().astimezone().isoformat(),'named_checks':len(checks),'checks':checks,'original_input_transcripts':len(raw),'output_transcripts':output_rows,'point_estimate_transcripts':accepted_rows,'all_matrix_entries_compared':int(recount.size),'source_xy_max_abs_error_um':max_xy_error,'rows_cells':len(cm),'QC_cells':len(a),'genes':M.shape[1],'raw_input_sha256':hashlib.sha256((old/'downloads/pilot_tumor1_region6.csv.gz').read_bytes()).hexdigest(),'scope':'data identity and reconstruction, not biological truth or unique segmentation'}
(P/'qc/PROSEG_IDENTITY_VALIDATION.json').write_text(json.dumps(rec,indent=2));print(json.dumps({k:v for k,v in rec.items() if k!='checks'},indent=2),flush=True)
