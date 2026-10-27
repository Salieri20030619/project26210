from common import *
import h5py
from anndata.io import read_elem
from scipy import sparse
def main():
 gate('independent_validation')
 z=np.load(P/'inputs/IDENTITY_PANEL_COUNTS.npz')
 raw=z['counts'];genes=z['genes'].tolist();ix={g:i for i,g in enumerate(genes)}
 a=read(P/'inputs/CELL_IDENTITY_AUDIT.tsv.gz');o=read(P/'inputs/ORIGINAL_CELL_METADATA.tsv.gz')
 assert a.cell_id.tolist()==z['cell_id'].tolist()==o.cell_id.tolist()
 pd.testing.assert_frame_equal(a[o.columns],o)
 panels=json.loads((P/'resources/IDENTITY_PANELS_FROZEN.json').read_text())
 for k,gs in panels.items():
  q=raw[:,[ix[g] for g in gs if g in ix]]
  assert np.array_equal((q>0).sum(1),a['det_'+k])
  assert np.allclose(np.log1p(q/z['total_counts'][:,None]*10000).mean(1),a['score_'+k],rtol=1e-12,atol=1e-12)
 for g,i in ix.items():assert np.allclose(np.log1p(raw[:,i]/z['total_counts']*10000),a['expr_'+g],rtol=1e-12,atol=1e-12)
 det=lambda gs:(raw[:,[ix[g] for g in gs]]>0).sum(1)
 anch=det(['DCN','LUM','PDGFRA']);coll=det(['COL1A1','COL1A2','COL3A1'])
 candidate=(anch+coll>=3)&(anch>0)&(coll>0)
 strict=(anch+coll>=4)&(anch>=2)&(coll>0)
 q=o[['analysis_keep','low_quality','scrublet_predicted_doublet']].astype(str).apply(lambda s:s.str.lower().map({'true':True,'false':False,'1':True,'0':False,'1.0':True,'0.0':False}))
 qc=q.notna().all(axis=1)&q.analysis_keep.fillna(False)&~q.low_quality.fillna(True)&~q.scrublet_predicted_doublet.fillna(True)&(o.n_genes_raw>=500)&(o.total_counts_raw>=1000)
 mural=(det(['RGS5','CSPG4','ABCC9','KCNJ8'])>=2)|((raw[:,ix['MYH11']]>0)&(raw[:,ix['CNN1']]>0))
 high=strict&qc&~mural
 for k in ['immune','endothelial','epithelial']:high &= det(panels[k])<2
 for k,b in [('fibro_candidate',candidate),('fibro_strict',strict),('qc_pass',qc),('strong_mural',mural),('high_priority_review',high)]:assert np.array_equal(a[k],b)
 review=read(P/'tables/CANDIDATE_CELL_REVIEW.tsv.gz')
 assert set(review.cell_id)==set(a.loc[candidate|a.cell_type_broad_v1.eq('Fibroblast'),'cell_id'])
 pt=read(P/'tables/PATIENT_CAF_REVIEW_CAPACITY.tsv')
 for r in pt.itertuples():
  g=a[(a.study==r.study)&(a.curated_patient_id==r.patient)&(a.source_tissue==r.tissue)]
  assert len(g)==r.cells and g.high_priority_review.sum()==r.high_priority_review and g.fibro_candidate.sum()==r.candidate and g.cell_type_broad_v1.eq('Fibroblast').sum()==r.original_fibro
 caps=read(P/'tables/REFERENCE_CAPACITY_DESCRIPTIVE.tsv')
 for r in caps.itertuples():
  g=pt[(pt.tissue==r.tissue)&(pt[r.screen]>=r.threshold)]
  assert len(g)==r.patients and g.study.nunique()==r.studies and g[r.screen].sum()==r.cells
 rng=np.random.default_rng(286)
 selected=set(rng.choice(len(a),200,replace=False).tolist())
 for _,g in a.groupby(['study','cell_type_broad_v1'],observed=True):
  selected.add(int(g.index[0]))
 selected=sorted(selected)
 with h5py.File(H,'r') as f:
  sourceids=f['obs']['cell_id'].asstr()[:];assert sourceids.tolist()==a.cell_id.tolist()
  for k in o.columns.drop('cell_id'):
   s=pd.Series(np.asarray(read_elem(f['obs'][k])))
   assert np.array_equal(s.fillna('__NA__').astype(str),o[k].fillna('__NA__').astype(str)),k
  gg=f['var']['gene'].asstr()[:].tolist();jj=[gg.index(g) for g in genes];x=f['layers/counts']
  for i in selected:
   lo,hi=x['indptr'][i:i+2]
   row=sparse.csr_matrix((x['data'][lo:hi],x['indices'][lo:hi],np.array([0,hi-lo])),shape=(1,len(gg)))
   assert np.array_equal(row[:,jj].toarray().ravel(),raw[i])
   assert row.sum()==z['total_counts'][i] and (row>0).sum()==z['n_genes_current'][i]
 src=json.loads((P/'qc/CELL_AUDIT_COMPLETE.json').read_text())['source'];st=H.stat()
 assert src['bytes']==st.st_size and src['mtime_ns']==st.st_mtime_ns
 ims=read(P/'tables/IMAGE_COORDINATE_AVAILABILITY.tsv')
 assert len(ims)==38 and ims.registered_dimensions_match.all() and ims.coordinate_finite.all() and (ims.spots==ims.spots_inside).all()
 poly=read(P/'tables/POLYGON_GEOMETRY_AND_PROVENANCE.tsv')
 assert poly.vertices_in_image.all() and poly.finite_vertices.all() and poly.image_dimensions_match.all() and (poly.area_px2>0).all()
 for path,g in poly.groupby('path'):assert sha(path)==g.sha256.iloc[0]
 for s,t,n in [('GLMF1',903,1998),('GLMF2',2481,488)]:
  g=poly[(poly['sample']==s)&poly.path.str.contains('wangmingyuan')]
  assert g.loc[g.label=='tumor','spots_inside'].sum()==t and g.loc[g.label=='NAT','spots_inside'].sum()==n
 assert not poly.eligible_as_independent_gold_standard.any()
 inventory=read(P/'tables/BOUNDARY_CANDIDATE_INVENTORY.tsv')
 decisions=[]
 for r in inventory.itertuples():
  f=Path(r.path);kind=r.kind
  if f.name=='annotation_manifest.json':
   j=json.loads(f.read_text());assert set(j['samples'])=={'GLMF1','GLMF2'};decision='coordinate_and_HTML_asset_manifest_no_new_boundary'
  elif kind=='TABLE_PROVENANCE_REQUIRED':
   d=pd.read_csv(f)
   if f.name=='pathology_annotation_summary.csv':
    assert set(d['sample'])=={'GLMF1','GLMF2'};decision='summary_of_same_two_AI_drafts'
   else:
    assert set(d.pathology_label.dropna()) <= {'tumor','NAT','unassigned'};decision='spot_labels_and_distance_derived_from_same_AI_draft'
  else:decision=kind
  decisions.append(dict(path=r.path,classification=decision,independent_tumor_stroma_boundary=False))
 tab(decisions,'BOUNDARY_ASSET_ADJUDICATION.tsv')
 assert sha(P/'PROTOCOL_BEFORE_AUDIT.md')==json.loads((P/'FREEZE.json').read_text())['protocol_sha256']
 for name in ['cell_identity','boundary']:assert (P/'logs'/f'{name}_exit.txt').read_text().strip()=='0'
 report=dict(status='PASS',all_cell_metadata_matches_source=True,all_panel_scores_and_screens_recomputed=True,all_patient_capacity_tables_recomputed=True,independent_source_count_rows=len(selected),candidate_review_complete=True,source_unchanged=True,images=38,coordinate_spots=int(ims.spots.sum()),all_spots_inside=True,polygon_geometry_and_label_counts_match=True,inventory_hits_adjudicated=len(decisions),protocol_hash_verified=True,new_biological_inference=False)
 dump(P/'qc/INDEPENDENT_VALIDATION.json',report)
 print(json.dumps(report),flush=True)
if __name__=='__main__':main()
