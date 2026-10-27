from common import *
from scipy import sparse
import h5py
from sklearn.metrics import adjusted_rand_score
def main():
 gate('validation')
 audit=R/'phase286_CAF_identity_boundary_audit_20260929_v1'
 o=read(audit/'inputs/CELL_IDENTITY_AUDIT.tsv.gz')
 m=read(P/'tables/STROMAL_CLUSTER_MEMBERSHIP.tsv.gz',dtype={'cluster_primary':str,'cluster_sensitive':str,'cluster_original_scvi':str})
 wanted=(o.cell_type_broad_v1.isin(['Fibroblast','Pericyte/SMC'])|o.fibro_candidate)&o.qc_pass
 assert m.cell_id.tolist()==o.loc[wanted,'cell_id'].tolist()
 for c in ['cell_type_broad_v1','curated_patient_id','study','sample_id','source_tissue']:
  assert m[c].tolist()==o.loc[wanted,c].tolist()
 raw=sparse.load_npz(P/'inputs/STROMAL_COUNTS.npz');genes=np.load(P/'inputs/GENES.npy');ids={g:i for i,g in enumerate(genes)}
 assert raw.shape==(len(m),len(genes))
 assert np.array_equal(np.asarray(raw.sum(1)).ravel(),m.total_counts_raw)
 hv=read(P/'tables/HVG_SELECTION.tsv',index_col=0)
 selected=[g for g,v in hv.highly_variable.items() if v and g not in ['TYMP','SPP1'] and not g.startswith(('MT-','RPL','RPS')) and g not in ['HBA1','HBA2','HBB','HBD','HBE1','HBG1','HBG2','HBM','HBQ1','HBZ']]
 result=json.loads((P/'qc/STROMAL_CLUSTER_COMPLETE.json').read_text())
 assert result['hvg']==len(selected)
 rr=np.random.default_rng(287).choice(len(m),120,replace=False)
 with h5py.File(H,'r') as h:
  assert np.array_equal(h['obs/cell_id'].asstr()[:][m.original_row],m.cell_id)
  x=h['layers/counts']
  for i in rr:
   r=int(m.original_row.iloc[i]);lo,hi=x['indptr'][r:r+2]
   row=sparse.csr_matrix((x['data'][lo:hi],x['indices'][lo:hi],[0,hi-lo]),shape=(1,len(genes)))
   assert (row!=raw[i]).nnz==0
 norm=raw.astype(float).multiply((1e4/m.total_counts_raw.to_numpy())[:,None]).tocsr();norm.data=np.log1p(norm.data)
 markers=read(P/'tables/CLUSTER_LINEAGE_MARKERS.tsv',dtype={'cluster':str})
 err=0.
 for (mode,cluster),tabx in markers.groupby(['mode','cluster']):
  ii=np.flatnonzero(m[mode].eq(cluster).to_numpy())
  for r in tabx.itertuples():
   col=ids[r.gene];mean=float(norm[ii,col].mean());det=float((raw[ii,col]>0).mean())
   err=max(err,abs(mean-r.mean_logCP10K));assert abs(det-r.detection)<1e-10
 assert err<2e-4
 iq=read(P/'tables/IMAGE_MORPHOLOGY_QC.tsv')
 isp=read(P/'tables/IMAGE_ONLY_SPOT_CLASSES.tsv.gz')
 for r in iq.itertuples():
  d=isp[isp['sample']==r.sample];g=np.load(P/'inputs'/f'{r.sample}_IMAGE_REGIONS.npz')
  for size in [64,96]:
   expect=g[f'grid{size}'][(d.y_hires_px.to_numpy()//size).astype(int),(d.x_hires_px.to_numpy()//size).astype(int)]
   assert np.array_equal(expect,d[f'morphology_class_{size}'])
  good=d.morphology_class_64.ge(0)&d.morphology_class_96.ge(0)
  assert good.sum()==r.cross_scale_common_spots
  if good.sum()>=30:assert abs(adjusted_rand_score(d.loc[good,'morphology_class_64'],d.loc[good,'morphology_class_96'])-r.cross_scale_ARI)<1e-12
  assert len(d)==r.spots and d.morphology_class_64.ge(0).sum()==r.assigned_spots_64
 for n in ['stromal','image']:assert (P/'logs'/f'{n}_exit.txt').read_text().strip()=='0'
 assert not iq.RNA_used.any() and not iq.validated_tumor_stroma_boundary.any()
 st=H.stat();assert st.st_size==result['source_bytes'] and st.st_mtime_ns==result['source_mtime_ns']
 assert sha(P/'PROTOCOL_BEFORE_ANALYSIS.md')==json.loads((P/'FREEZE.json').read_text())['sha256']
 dump(P/'qc/INDEPENDENT_VALIDATION.json',dict(status='PASS',candidate_set_and_original_metadata_match=True,counts_totals_match=True,source_rows_checked=len(rr),lineage_mean_max_error=err,all_marker_detections_match=True,image_spot_classes_reconstructed=True,all_cross_scale_ARI_recomputed=True,source_unchanged=True,protocol_unchanged=True,images=len(iq),cells=len(m)))
 print('INDEPENDENT_VALIDATION_PASS',flush=True)
if __name__=='__main__':main()
