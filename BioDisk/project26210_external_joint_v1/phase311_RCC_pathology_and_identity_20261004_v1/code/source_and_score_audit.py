from pathlib import Path
import json,datetime
import numpy as np,pandas as pd,anndata as ad,scanpy as sc
from scipy import sparse
P=Path(__file__).resolve().parents[1];O=P.parent/'phase310_RCC_mapped_spatial_inputs_20261004_v1';q=P/'resources/proseg_region6_v1'
a=ad.read_h5ad(P/'resources/proseg_region6_identity.h5ad');frozen=pd.read_csv(P/'tables/CELL_IDENTITY_FROZEN.tsv.gz',sep='\t',dtype={'original_cell_id':str})
oldmeta=pd.read_csv(O/'tables/PILOT_CELL_METADATA.tsv',sep='\t',dtype={'cell_id':str}).set_index('cell_id').loc[a.obs_names]
fg=pd.read_csv(O/'tables/PILOT_FEATURE_AUDIT.tsv',sep='\t');s0=int(fg.loc[fg.gene.eq('SPP1'),'matrix_col'].iloc[0]);gn=list(a.var_names);sj=gn.index('SPP1')
A=sparse.load_npz(O/'resources/pilot_assigned_counts_qv20.npz')[oldmeta.matrix_row.to_numpy(int)]
N=sparse.load_npz(O/'resources/pilot_nuclear_counts_qv20.npz')[oldmeta.matrix_row.to_numpy(int)]
C=a.layers['maxpost_counts'].toarray()
assays={'Proseg':(C[:,sj],C.sum(axis=1)),'OriginalAssigned':(A[:,s0].toarray().ravel(),np.asarray(A.sum(axis=1)).ravel()),'OriginalNuclear':(N[:,s0].toarray().ravel(),np.asarray(N.sum(axis=1)).ravel())}
rows=[]
for assay,(sp,lib) in assays.items():
 for lineage in sorted(frozen.reviewed_cluster_lineage.unique()):
  ix=frozen.reviewed_cluster_lineage.eq(lineage).to_numpy()
  rows.append({'assay':assay,'population':lineage,'cell_count':int(ix.sum()),'SPP1_detected_cells':int((sp[ix]>0).sum()),'SPP1_detection_fraction':float((sp[ix]>0).mean()),'SPP1_counts':int(sp[ix].sum()),'fraction_of_total_SPP1_counts':float(sp[ix].sum()/sp.sum())})
 for role in ['mature_myeloid','strict_fibroblast','strict_CD8']:
  ix=frozen[role].to_numpy(bool)
  rows.append({'assay':assay,'population':role,'cell_count':int(ix.sum()),'SPP1_detected_cells':int((sp[ix]>0).sum()),'SPP1_detection_fraction':float((sp[ix]>0).mean()) if ix.any() else None,'SPP1_counts':int(sp[ix].sum()),'fraction_of_total_SPP1_counts':float(sp[ix].sum()/sp.sum())})
pd.DataFrame(rows).to_csv(P/'tables/SPP1_CELLULAR_SOURCE_AUDIT.tsv',sep='\t',index=False)
# Explicit numeric-boundary sensitivity, triggered by documented terminal-tick issue; no threshold search.
cm=pd.read_csv(q/'cell-metadata.csv.gz',dtype={'original_cell_id':str});M=pd.read_csv(q/'maxpost-counts.csv.gz');gm={g:i for i,g in enumerate(M)}
delta=np.zeros(M.shape,np.int32)
for d in pd.read_csv(q/'transcript-metadata.csv.gz',usecols=['assignment','probability','gene'],chunksize=250000):
 pr=d.probability.to_numpy(np.float32);ci=d.assignment.to_numpy(np.uint64);good=(ci<len(cm))&(pr>np.float32(.1))&(pr<=np.float32(.100011))
 gi=d.gene.map(gm).to_numpy(int);np.add.at(delta,(ci[good].astype(int),gi[good]),1)
assert int(delta.sum())==923
assert (M.to_numpy()-delta>=0).all()
idx=cm.set_index('original_cell_id').loc[a.obs_names,'cell'].to_numpy(int)
dc=delta[idx];newC=C-dc;keep=(newC.sum(axis=1)>=20)&((newC>0).sum(axis=1)>=5)
rec={'status':'PASS','reason':'conservative removal of all model scores within one terminal tick of .1 threshold; not claimed exact posterior correction','removed_transcripts':int(delta.sum()),'removed_fraction':float(delta.sum()/M.to_numpy().sum()),'QC_cells_with_any_count_change':int((dc.sum(axis=1)>0).sum()),'QC_cells_lost_after_conservative_change':int((~keep).sum()),'SPP1_removed_all':int(delta[:,sj].sum()),'SPP1_removed_in_strict_myeloid':int(dc[frozen.mature_myeloid,sj].sum()),'strict_myeloid_SPP1_detection_lost':int(((C[:,sj]>0)&(newC[:,sj]==0)&frozen.mature_myeloid.to_numpy()).sum()),'labels':'identity clusters held fixed; per-cell membership and technical QC reapplied'}
b=a[keep].copy();b.X=sparse.csr_matrix(newC[keep].astype(np.float32));b.layers['maxpost_counts']=b.X.copy();b.obs['counts']=newC[keep].sum(axis=1);b.obs['genes']=(newC[keep]>0).sum(axis=1)
sc.pp.normalize_total(b,target_sum=1e4);sc.pp.log1p(b);b.uns['numerical_boundary_sensitivity']='Scores(.1,.100011] excluded; original cluster identity frozen; no re-clustering'
b.write_h5ad(P/'resources/proseg_region6_identity_score_boundary.h5ad',compression='gzip')
pd.DataFrame(delta,columns=M.columns).to_csv(P/'resources/score_boundary_removed_counts.csv.gz',index=False,compression='gzip')
(P/'qc/NUMERICAL_BOUNDARY_AUDIT.json').write_text(json.dumps(rec,indent=2))
print(json.dumps(rec,indent=2));print(pd.DataFrame(rows).query("assay=='Proseg'").to_string(index=False),flush=True)
