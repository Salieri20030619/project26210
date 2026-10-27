from common import *
import h5py
from anndata.io import read_elem
from scipy import sparse,stats,linalg
from scipy.optimize import lsq_linear
from sklearn.neighbors import BallTree
def main():
 gate('validate_spatial')
 pb=read(P/'tables/REFERENCE_PSEUDOBULK_REGISTER.tsv');vv=read(P/'inputs/REFERENCE_PSEUDOBULK_CPM.tsv.gz');genes=vv.columns.tolist();v=vv.to_numpy();classes=['Fibroblast_like','Mural','Endothelial','Epithelial','Macrophage','Tcell']
 qual=read(P/'tables/SPATIAL_REFERENCE_QUALIFICATION.tsv');eff=read(P/'tables/NEIGHBORHOOD_SECTION_EFFECTS.tsv');maxrho=0.;maxcoef=0.;maxmean=0.;models=0
 for r in qual.itertuples():
  if not r.qualified:continue
  z=np.load(P/'inputs'/f'{r.sample}_SPATIAL_PANEL.npz');gs=z['genes'].tolist();jj=[genes.index(g) for g in gs]
  cols=[]
  for cl in classes:
   ii=pb.label.eq(cl);cols.append(np.mean([v[ii&pb.study.eq(st)][:,jj].mean(0) for st in pb.loc[ii,'study'].unique()],axis=0))
  A=np.array(cols).T;s=np.sqrt(np.maximum(A.mean(1),1.));A=A/s[:,None]
  for i in np.linspace(0,len(z['CPM'])-1,7,dtype=int):
   fit=lsq_linear(A,z['CPM'][i]/s,bounds=(0,np.inf),method='bvls',tol=1e-12)
   assert fit.success
   err=float(np.max(np.abs(fit.x-z['coefficients'][i])));maxcoef=max(maxcoef,err);assert err<1e-6
  d=read(R/'phase234_part4_replan_execution_20260921_v1/inputs'/f'{r.sample}_INPUTS.tsv.gz');assert d.spot_id.tolist()==z['spot_id'].tolist()
  d['fibro_coefficient']=z['coefficients'][:,0]
  # Independently verify spot counts extraction in three representative source formats.
  if r.sample in ['GLMF1','H46_T','GSM5924030']:
   register=read(R/'phase234_part4_replan_execution_20260921_v1/tables/INPUT_SECTION_REGISTER.tsv')
   f=register.loc[register['sample']==r.sample,'input_path'].iloc[0]
   with h5py.File(f,'r') as h:
    names=np.asarray(read_elem(h['var/gene_symbol_original'])) if 'gene_symbol_original' in h['var'] else np.asarray(read_elem(h['var'][h['var'].attrs['_index']]))
    ids=np.asarray(read_elem(h['obs'][h['obs'].attrs['_index']]))
    node=h['layers/counts'] if 'layers/counts' in h else h['X']
    for i in [0,len(d)//2,len(d)-1]:
     k=int(np.flatnonzero(ids==d.spot_id.iloc[i])[0]);lo,hi=node['indptr'][k:k+2];row=sparse.csr_matrix((node['data'][lo:hi],node['indices'][lo:hi],[0,hi-lo]),shape=(1,len(names)))
     expect=row[:,[np.flatnonzero(names==g)[0] for g in gs]].toarray().ravel()/z['totals'][i]*1e6
     assert np.allclose(expect,z['CPM'][i],rtol=1e-12,atol=1e-12)
  tree=BallTree(d[['physical_x_um','physical_y_um']].to_numpy())
  for rad in [200,400]:
   n=read(P/'inputs'/f'{r.sample}_NEIGHBORHOOD_{rad}.tsv.gz')
   assert n.spot_id.tolist()==d.spot_id.tolist()
   cols=['log1p_UMI','total_macro_fraction','malignant_fraction','endothelial_fraction','pericyte_fraction','Tcell_fraction','SPP1_RNA','fibro_coefficient']
   # Independent neighbor search on deterministic representative spots, no self.
   for k in np.linspace(0,len(d)-1,20,dtype=int):
    neighbors=tree.query_radius(d[['physical_x_um','physical_y_um']].iloc[[k]].to_numpy(),r=rad)[0];neighbors=neighbors[neighbors!=k]
    assert len(neighbors)==n.neighbors.iloc[k]
    if len(neighbors):
     expect=d.iloc[neighbors][cols].mean().to_numpy();got=n.iloc[k][['near_'+x for x in cols]].to_numpy(float)
     error=np.nanmax(np.abs(expect-got));maxmean=max(maxmean,float(error));assert np.allclose(expect,got,atol=1e-10,rtol=1e-10,equal_nan=True)
   joint=(stats.rankdata(n.near_fibro_coefficient)/(len(n)+1))*(stats.rankdata(n.near_SPP1_RNA)/(len(n)+1));assert np.allclose(joint,n.joint_neighborhood_rank)
   for e in eff[(eff['sample']==r.sample)&(eff.radius==rad)].itertuples():
    core=['log1p_UMI','total_macro_fraction','malignant_fraction'];vc=['endothelial_fraction','pericyte_fraction'] if e.model=='vascular_sensitivity' else []
    X=np.column_stack([d[core+vc].to_numpy(),n[['near_'+c for c in core+vc]].to_numpy()])
    if e.family=='immune_distribution':X=np.column_stack([X,n.near_Tcell_fraction.to_numpy()]);x=joint;y=d.CD8_expression.to_numpy()
    elif e.endpoint=='SPP1_near_fibro':x=d.SPP1_RNA.to_numpy();y=n.near_fibro_coefficient.to_numpy()
    else:x=d.TYMP_RNA.to_numpy();y=n.near_fibro_coefficient.to_numpy()
    good=n.eligible.to_numpy()&np.isfinite(x)&np.isfinite(y)&np.isfinite(X).all(1)
    assert int(good.sum())==e.spots
    if not np.isfinite(e.effect):continue
    qr=stats.rankdata(X[good],axis=0);qr=qr[:,qr.std(0)>1e-9];qr=(qr-qr.mean(0))/qr.std(0);design=np.column_stack([np.ones(good.sum()),qr])
    pair=np.column_stack([stats.rankdata(x[good]),stats.rankdata(y[good])]);coef=linalg.lstsq(design,pair,lapack_driver='gelsy',cond=max(design.shape)*np.finfo(float).eps)[0];resid=pair-design@coef;rho=float(np.corrcoef(resid.T)[0,1]);maxrho=max(maxrho,abs(rho-e.effect));assert abs(rho-e.effect)<1e-9,(r.sample,rad,e.model,e.endpoint,rho,e.effect,np.linalg.matrix_rank(design),design.shape);models+=1
  print('VALIDATED',r.sample,flush=True)
 pt=read(P/'tables/NEIGHBORHOOD_PATIENT_EFFECTS.tsv');keys=['dataset','radius','model','family','endpoint']
 for row in pt.to_dict('records'):
  mask=np.ones(len(eff),bool)
  for k in keys+['patient']:mask &= eff[k].eq(row[k])
  vals=eff.loc[mask,'effect'].dropna().to_numpy();expect=np.tanh(np.mean(np.arctanh(np.clip(vals,-.999999,.999999)))) if len(vals) else np.nan
  assert np.isclose(expect,row['effect'],equal_nan=True)
 su=read(P/'tables/NEIGHBORHOOD_COHORT_SUMMARY.tsv')
 for ix,row in su.iterrows():
  mask=np.ones(len(pt),bool)
  for k in keys:mask &= pt[k].eq(row[k])
  vals=pt.loc[mask,'effect'].dropna().to_numpy();assert len(vals)==row.patients
  assert np.isclose(np.median(vals),row.median_effect,equal_nan=True)
  p=stats.wilcoxon(vals).pvalue if len(vals)>=10 and np.any(vals!=0) else (1. if len(vals)>=10 else np.nan)
  assert np.isclose(p,row.P,equal_nan=True,atol=1e-12)
 for _,g in su.groupby(['dataset','radius','model','family']):
  g=g[g.P.notna()]
  if len(g):assert np.allclose(stats.false_discovery_control(g.P.to_numpy(),method='bh'),g.FDR.to_numpy(),atol=1e-12)
 for n,freeze in [('REFERENCE_PROTOCOL_BEFORE_TEST.md','REFERENCE_PROTOCOL_FREEZE.json'),('SPATIAL_PROTOCOL_BEFORE_TEST.md','SPATIAL_PROTOCOL_FREEZE.json')]:
  assert sha(P/n)==json.loads((P/'qc'/freeze).read_text())['sha256']
 for r in json.loads((P/'qc/SPATIAL_INPUT_MANIFEST.json').read_text()):
  st=Path(r['path']).stat();assert st.st_size==r['bytes'] and st.st_mtime_ns==r['mtime_ns']
 assert (P/'logs/spatial_pilot_exit.txt').read_text().strip()=='0'
 dump(P/'qc/SPATIAL_INDEPENDENT_VALIDATION.json',dict(status='PASS',sections=len(qual),NNLS_spot_fits_checked=len(qual)*7,NNLS_max_coefficient_difference=maxcoef,neighborhood_mean_max_difference=maxmean,partial_correlations_checked=models,partial_rho_max_difference=maxrho,all_patient_aggregates_P_FDR_match=True,source_spatial_counts_formats_checked=3,all_source_fingerprints_unchanged=True,protocols_unchanged=True))
 print('SPATIAL_VALIDATION_PASS',flush=True)
if __name__=='__main__':main()
