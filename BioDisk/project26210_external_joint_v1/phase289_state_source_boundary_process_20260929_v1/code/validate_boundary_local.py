from common import *
from scipy import sparse,stats,linalg
from scipy.sparse.csgraph import connected_components
from sklearn.neighbors import BallTree
def independent_masks(d,G,patch):
 coords={(int(r.x)//64,int(r.y)//64):[r.H_mean,r.E_mean,r.gray_mean,r.sat_mean] for r in patch.itertuples()}
 X=np.array([coords.get((int(x)//64,int(y)//64),[np.nan]*4) for x,y in d[['x_hires_px','y_hires_px']].to_numpy()])
 good=np.isfinite(X).all(1)&np.isfinite(d[['p_tumor_anchor','p_epithelial_anchor','p_normal_kidney']]).all(1).to_numpy()
 if good.sum()<100:return good,{}
 X[good]=(X[good]-X[good].mean(0))/np.maximum(X[good].std(0),1e-8)
 rows,cols=G.nonzero();valid=good[rows]&good[cols];rows=rows[valid];cols=cols[valid]
 dis=np.sum((X[rows]-X[cols])**2,axis=1);scale=np.median(dis[dis>0]) if (dis>0).any() else 1.
 w=np.exp(-dis/(2*max(scale,1e-8)));contrast=(d.p_tumor_anchor-d.p_normal_kidney).to_numpy()
 den=np.bincount(rows,weights=w,minlength=len(d));num=np.bincount(rows,weights=w*np.nan_to_num(contrast[cols]),minlength=len(d))
 sm=.5*contrast+.5*num/np.maximum(den,1e-12)
 seed=good&d.p_tumor_anchor.ge(d.loc[good,'p_tumor_anchor'].quantile(.75)).to_numpy()&d.p_epithelial_anchor.ge(d.loc[good,'p_epithelial_anchor'].quantile(.5)).to_numpy()&d.d_tumor_anchor.ge(2).to_numpy()&(contrast>0)
 epi=d.p_epithelial_anchor.ge(d.loc[good,'p_epithelial_anchor'].quantile(.25)).to_numpy()
 out={}
 for name,score in [('morphology_supported',sm),('expression_sensitivity',contrast)]:
  cand=good&epi&(score>0);visited=np.zeros(len(d),bool);mask=visited.copy()
  for s in np.flatnonzero(seed&cand):
   if visited[s]:continue
   stack=[int(s)];visited[s]=True;comp=[]
   while stack:
    i=stack.pop();comp.append(i)
    for j in G.indices[G.indptr[i]:G.indptr[i+1]]:
     if cand[j] and not visited[j]:visited[j]=True;stack.append(int(j))
   if len(comp)>=10:mask[comp]=True
  out[name]=mask
 return good,out
def main():
 gate('validate_boundary_local')
 q=read(P/'tables/BOUNDARY_QUALIFICATION.tsv');ef=read(P/'tables/LOCAL_SECTION_EFFECTS.tsv')
 panels=json.loads((P/'PANELS_FROZEN.json').read_text());maxavg=0.;maxrho=0.;checked=0;diagnostics=[]
 for sample,qq in q.groupby('sample',sort=False):
  d=read(P/'inputs'/f'{sample}_SPOT_PROGRAMS.tsv.gz');z=np.load(P/'inputs'/f'{sample}_PANELS.npz')
  assert d.spot_id.tolist()==z['spot_id'].tolist()
  log=np.log1p(z['counts']/z['total_counts'][:,None]*1e4);idx={g:i for i,g in enumerate(z['genes'].tolist())}
  for name,gs in panels.items():
   inds=[idx[g] for g in gs if g in idx]
   if len(inds)>=.8*len(gs):
    assert np.allclose(log[:,inds].mean(1),d['p_'+name],equal_nan=True,atol=1e-12)
    assert np.array_equal((z['counts'][:,inds]>0).sum(1),d['d_'+name])
  G=sparse.load_npz(R/'phase234_part4_replan_execution_20260921_v1/graphs'/f'{sample}_ring.npz').tocsr()
  patch=read(R/'phase287_stromal_refinement_image_morphology_20260929_v1/tables'/f'{sample}_IMAGE_PATCHES_64.tsv.gz')
  good,masks=independent_masks(d,G,patch);assert np.array_equal(good,d.boundary_input_valid)
  tree=BallTree(d[['physical_x_um','physical_y_um']].to_numpy());graph=sparse.diags(good.astype(float))@G@sparse.diags(good.astype(float))
  _,components=connected_components(graph,directed=False)
  for definition,inside in masks.items():
   assert np.array_equal(inside,d['inside_'+definition])
   outside=good&~inside
   row=qq[qq.definition==definition].iloc[0]
   rr,cc=graph.nonzero();cross=(rr<cc)&(inside[rr]!=inside[cc]);assert cross.sum()==row.crossing_edges
   assert inside.sum()==row.inside and outside.sum()==row.outside
   for col in ['malignant_fraction','p_heldout_epithelial','p_tumor_anchor','p_normal_kidney']:
    diagnostics.append(dict(sample=sample,definition=definition,readout=col,inside_mean=float(d.loc[inside,col].mean()),outside_mean=float(d.loc[outside,col].mean()),interpretation='descriptive_consistency_not_independent_pathology'))
   if not cross.any():continue
   memberships=json.loads((P/'inputs'/f'{sample}_{definition}_MEMBERSHIPS.json').read_text())
   md={(r['radius'],r['edge']):r for r in memberships}
   macro=d.total_macro_fraction.ge(d.total_macro_fraction.quantile(.75))&d.malignant_fraction.le(d.malignant_fraction.quantile(.5))
   cols=['g_TYMP','g_SPP1','p_CAF_activation','p_transition_RCC','p_TGF_response','p_ECM_organization','p_ECM_remodeling','p_source_identity','CD8_expression','Tcell_fraction','log1p_UMI','total_macro_fraction','malignant_fraction','endothelial_fraction','pericyte_fraction']
   for radius in [200,400]:
    e=read(P/'inputs'/f'{sample}_{definition}_{radius}_EDGES.tsv.gz')
    for k in np.linspace(0,len(e)-1,min(8,len(e)),dtype=int):
     r=e.iloc[k];nodein=int(r.inner_index);nodeout=int(r.outer_index)
     assert inside[nodein] and outside[nodeout] and G[nodein,nodeout]>0
     neighbor=tree.query_radius(np.array([[r.x,r.y]]),r=radius)[0]
     neighbor=neighbor[good[neighbor]&(components[neighbor]==components[nodein])]
     ii=np.sort(neighbor[inside[neighbor]]);oo=np.sort(neighbor[~inside[neighbor]]);mo=oo[macro.iloc[oo].to_numpy()]
     for side,inds in [('in',ii),('out',oo),('macro_out',mo)]:
      rec=md[(radius,int(r.edge))];stored=rec['inside' if side=='in' else ('outside' if side=='out' else 'macro_out')]
      assert np.array_equal(inds,np.sort(stored))
      if len(inds)>=3:
       expected=d.iloc[inds][cols].mean().to_numpy();got=r[[side+'_'+c for c in cols]].to_numpy(float)
       err=np.abs(expected-got);maxavg=max(maxavg,float(np.nanmax(err)));assert np.allclose(expected,got,atol=1e-10,rtol=1e-10,equal_nan=True)
    for rr in ef[(ef['sample']==sample)&(ef.definition==definition)&(ef.radius==radius)].itertuples():
     core=['in_log1p_UMI','out_log1p_UMI','in_total_macro_fraction','out_total_macro_fraction','in_malignant_fraction','out_malignant_fraction','log_side_count_ratio']
     if rr.model=='vascular_sensitivity':core+=['in_endothelial_fraction','out_endothelial_fraction','in_pericyte_fraction','out_pericyte_fraction']
     xs={'tumor_side_TYMP_CAF':'in_g_TYMP','macro_region_TYMP_CAF':'macro_out_g_TYMP','macro_region_SPP1_CAF':'macro_out_g_SPP1','joint_CAF_SPP1_inner_CD8':'joint_rank'}
     y='in_CD8_expression' if rr.family=='immune_distribution' else 'out_p_CAF_activation'
     if rr.family=='immune_distribution':core+=['out_CD8_expression','out_Tcell_fraction']
     x=e[xs[rr.endpoint]].to_numpy();yy=e[y].to_numpy();cov=e[core].to_numpy()
     use=e.geometry_eligible.to_numpy()&bool(row.eligible)&np.isfinite(x)&np.isfinite(yy)&np.isfinite(cov).all(1)
     assert int(use.sum())==rr.edges
     if not np.isfinite(rr.effect):continue
     C=stats.rankdata(cov[use],axis=0);C=C[:,C.std(0)>1e-9];C=(C-C.mean(0))/C.std(0);D=np.column_stack([np.ones(use.sum()),C])
     pair=np.column_stack([stats.rankdata(x[use]),stats.rankdata(yy[use])]);coef=linalg.lstsq(D,pair,cond=max(D.shape)*np.finfo(float).eps,lapack_driver='gelsy')[0];res=pair-D@coef;rho=float(np.corrcoef(res.T)[0,1])
     maxrho=max(maxrho,abs(rho-rr.effect));assert abs(rho-rr.effect)<1e-9,(sample,rr.endpoint,rho,rr.effect);checked+=1
  print('BOUNDARY_VALIDATED',sample,flush=True)
 tab(diagnostics,'BOUNDARY_EXTERNAL_CONSISTENCY_DESCRIPTIVE.tsv')
 pt=read(P/'tables/LOCAL_PATIENT_EFFECTS.tsv');su=read(P/'tables/LOCAL_COHORT_SUMMARY.tsv');keys=['dataset','definition','radius','model','family','endpoint']
 for r in pt.to_dict('records'):
  keep=np.ones(len(ef),bool)
  for k in keys+['patient']:keep &= ef[k].eq(r[k])
  vals=ef.loc[keep,'effect'].dropna().to_numpy();v=np.tanh(np.mean(np.arctanh(np.clip(vals,-.999999,.999999)))) if len(vals) else np.nan
  assert np.isclose(v,r['effect'],equal_nan=True)
 for r in su.to_dict('records'):
  keep=np.ones(len(pt),bool)
  for k in keys:keep &= pt[k].eq(r[k])
  vals=pt.loc[keep,'effect'].dropna().to_numpy();assert len(vals)==r['patients']
  p=stats.wilcoxon(vals).pvalue if len(vals)>=10 and np.any(vals!=0) else (1. if len(vals)>=10 else np.nan)
  assert np.isclose(p,r['P'],equal_nan=True)
 for _,g in su.groupby(['dataset','definition','radius','model','family']):
  z=g[g.P.notna()]
  if len(z):assert np.allclose(stats.false_discovery_control(z.P.to_numpy()),z.FDR)
 for r in json.loads((P/'qc/SPATIAL_INPUT_MANIFEST.json').read_text()):
  st=Path(r['path']).stat();assert st.st_size==r['bytes'] and st.st_mtime_ns==r['mtime_ns']
 for n in ['boundary','local_tests']:assert (P/'logs'/f'{n}_exit.txt').read_text().strip()=='0'
 dump(P/'qc/BOUNDARY_LOCAL_INDEPENDENT_VALIDATION.json',dict(status='PASS',all_boundary_masks_independently_reconstructed=True,source_score_arrays_match=True,neighbor_average_max_error=maxavg,local_partial_correlations_checked=checked,local_rho_max_error=maxrho,all_patient_aggregates_P_FDR_match=True,source_inputs_unchanged=True))
 print('BOUNDARY_LOCAL_VALIDATION_PASS',flush=True)
if __name__=='__main__':main()
