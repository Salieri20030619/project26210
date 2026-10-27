from pathlib import Path
import json,hashlib,datetime,platform,importlib.metadata as im
import h5py,numpy as np,pandas as pd,anndata as ad
from anndata.io import read_elem,sparse_dataset
from scipy import sparse
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');O=R/'01_input_audit/run_v1';O.mkdir();Q=R/'02_myeloid_qc/run_v1';Q.mkdir()
SRC=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad')
def log(x):print(datetime.datetime.now().isoformat(),x,flush=True)
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as z:
  for b in iter(lambda:z.read(32*1024*1024),b''):h.update(b)
 return h.hexdigest()
log('input structure and identity')
before=SRC.stat()
with h5py.File(SRC,'r') as h:
 obs=read_elem(h['obs']);var=read_elem(h['var'])
 assert len(obs)==258120 and obs.sample_id.nunique()==75 and obs.curated_patient_id.nunique()==40
 assert set(obs.study.astype(str))=={'GSE156632','GSE159115','GSE178481','GSE207493','GSE242299'}
 mask=obs.cell_type_broad_v1.astype(str).isin(['Macrophage','Monocyte']);idx=np.flatnonzero(mask)
 meta=obs.loc[mask].copy();assert meta.index.is_unique and var.index.is_unique
 log('extract raw counts only for myeloid cells')
 x=sparse_dataset(h['layers/counts'])[idx,:].tocsr()
 oldrep=h['obsm/X_scVI'][idx,:] if 'X_scVI' in h['obsm'] else None
assert x.shape==(50761,37795)
assert np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.rint(x.data),atol=0,rtol=0)
rs=np.asarray(x.sum(axis=1)).ravel();assert np.allclose(rs,meta.total_counts_raw.to_numpy(),rtol=0,atol=1e-6)
assert meta.analysis_keep.astype(bool).all()
a=ad.AnnData(X=x,obs=meta,var=var)
if oldrep is not None:a.obsm['X_scVI_parent_reference_only']=oldrep
# Descriptive lineage checks only; SPP1/C1QC never define selection.
panels={'macrophage':['C1QA','C1QB','C1QC','APOE','CSF1R','CD68','LST1','TYROBP','FCER1G','CTSB'],'monocyte':['FCN1','S100A8','S100A9','VCAN','CTSS','LYZ'],'T_cell':['CD3D','CD3E','TRAC'],'epithelial':['EPCAM','KRT8','KRT18','KRT19'],'B_cell':['CD79A','MS4A1','CD37'],'endothelial':['PECAM1','VWF','KDR'],'mast':['TPSAB1','TPSB2','KIT']}
checks=meta[['study','sample_id','curated_patient_id','source_tissue','cell_type_broad_v1','annotation_confidence_v1','annotation_ambiguous_v1','pct_counts_mt_raw','n_genes_raw','total_counts_raw','scrublet_score','scrublet_predicted_doublet']].copy()
norm=x.astype(np.float64).multiply(10000/rs[:,None]).tocsr();norm.data=np.log1p(norm.data)
coverage=[]
for key,genes in panels.items():
 use=[g for g in genes if g in var.index];ii=var.index.get_indexer(use)
 checks[key+'_mean_log1pCP10K']=np.asarray(norm[:,ii].mean(axis=1)).ravel();checks[key+'_n_detected']=np.asarray((x[:,ii]>0).sum(axis=1)).ravel();coverage.append({'panel':key,'requested':';'.join(genes),'present':';'.join(use)})
checks['mixed_lineage_review_flag']=False
for key in ['T_cell','epithelial','B_cell','endothelial','mast']:
 flag=(checks[key+'_n_detected']>=2)&(checks[key+'_mean_log1pCP10K']>checks.macrophage_mean_log1pCP10K)&(checks[key+'_mean_log1pCP10K']>checks.monocyte_mean_log1pCP10K)
 checks.loc[flag,'mixed_lineage_review_flag']=True
checks.to_csv(Q/'cell_lineage_and_qc_audit.tsv.gz',sep='\t',compression='gzip',index_label='cell_id')
pd.DataFrame(coverage).to_csv(Q/'lineage_panel_coverage.tsv',sep='\t',index=False)
pd.crosstab(meta.cell_type_broad_v1,meta.study).to_csv(O/'lineage_by_study.tsv',sep='\t')
pd.crosstab(meta.cell_type_broad_v1,meta.curated_patient_id).to_csv(O/'lineage_by_patient.tsv',sep='\t')
pd.crosstab(meta.cell_type_broad_v1,meta.sample_id).to_csv(O/'lineage_by_library.tsv',sep='\t')
pd.crosstab(meta.cell_type_broad_v1,meta.source_tissue).to_csv(O/'lineage_by_tissue.tsv',sep='\t')
meta.reset_index(names='cell_id').to_csv(O/'myeloid_metadata.tsv.gz',sep='\t',index=False,compression='gzip')
log('write isolated raw-count myeloid checkpoint')
f=O/'myeloid_raw_counts_v1.h5ad';a.write_h5ad(str(f)+'.part',compression='gzip');Path(str(f)+'.part').replace(f)
summary=[]
for group,g in checks.groupby('cell_type_broad_v1',observed=True):
 summary.append({'lineage':str(group),'n_cells':len(g),'n_patients':g.curated_patient_id.nunique(),'n_libraries':g.sample_id.nunique(),'ambiguous':int(g.annotation_ambiguous_v1.astype(bool).sum()),'mixed_lineage_review_flags':int(g.mixed_lineage_review_flag.sum()),'predicted_doublets':int(g.scrublet_predicted_doublet.astype(bool).sum()),'confidence':g.annotation_confidence_v1.value_counts().to_dict()})
log('hash source and new checkpoint')
source_hash=sha(SRC);assert before.st_size==SRC.stat().st_size and before.st_mtime_ns==SRC.stat().st_mtime_ns
manifest=[{'role':'immutable_source','path':str(SRC),'bytes':before.st_size,'sha256':source_hash},{'role':'new_myeloid_checkpoint','path':str(f),'bytes':f.stat().st_size,'sha256':sha(f)}]
pd.DataFrame(manifest).to_csv(O/'INPUT_OUTPUT_SHA256.tsv',sep='\t',index=False)
versions={k:im.version(k) for k in ['anndata','scanpy','numpy','scipy','pandas','h5py','harmonypy','leidenalg','igraph']};(O/'SOFTWARE_VERSIONS.json').write_text(json.dumps(versions,indent=2))
status={'status':'INPUT_EXTRACTION_AND_QC_AUDIT_COMPLETE_PENDING_LINEAGE_REVIEW','source_cells':len(obs),'myeloid_cells':len(meta),'genes':len(var),'summary':summary,'source_sha256':source_hash,'legacy_labels_used_for_selection':False,'next':'Review flagged cells and freeze clustering input before stage2'}
(O/'STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2));print(json.dumps(status,ensure_ascii=False),flush=True)
