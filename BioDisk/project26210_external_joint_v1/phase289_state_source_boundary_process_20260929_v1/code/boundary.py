from common import *
import h5py
from anndata.io import read_elem
from scipy import sparse,stats
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
def extract(r,panels):
 d=read(R/'phase234_part4_replan_execution_20260921_v1/inputs'/f'{r.sample}_INPUTS.tsv.gz')
 with h5py.File(r.input_path,'r') as h:
  genes=np.asarray(read_elem(h['var/gene_symbol_original'])) if 'gene_symbol_original' in h['var'] else np.asarray(read_elem(h['var'][h['var'].attrs['_index']]))
  ids=np.asarray(read_elem(h['obs'][h['obs'].attrs['_index']])).astype(str);inds=pd.Index(ids).get_indexer(d.spot_id);assert (inds>=0).all()
  raw=read_elem(h['layers/counts'] if 'layers/counts' in h else h['X']).tocsr();assert (raw.data>=0).all() and np.allclose(raw.data,np.round(raw.data))
  total=np.asarray(read_elem(h['obs/total_counts'])) if 'total_counts' in h['obs'] else np.asarray(raw.sum(1)).ravel()
  union=sorted(set(sum(panels.values(),[])));present=[g for g in union if g in set(genes)];jj=[int(np.flatnonzero(genes==g)[0]) for g in present]
  x=raw[:,jj].toarray()[inds].astype(float);tot=total[inds];log=np.log1p(x/tot[:,None]*1e4);ix={g:i for i,g in enumerate(present)}
 cov=[]
 for name,gs in panels.items():
  kept=[g for g in gs if g in ix];ok=len(kept)>=.8*len(gs)
  d['p_'+name]=log[:,[ix[g] for g in kept]].mean(1) if ok else np.nan
  d['d_'+name]=(x[:,[ix[g] for g in kept]]>0).sum(1) if ok else np.nan
  cov.append(dict(sample=r.sample,program=name,present=len(kept),requested=len(gs),missing=';'.join(g for g in gs if g not in ix),eligible=ok))
 for g in ['TYMP','SPP1','TGFB1','CTSB']:d['g_'+g]=log[:,ix[g]] if g in ix else np.nan
 np.savez_compressed(P/'inputs'/f'{r.sample}_PANELS.npz',spot_id=d.spot_id.to_numpy(str),genes=np.array(present),counts=x,total_counts=tot)
 return d,cov
def define(d,G,patch):
 coords={(int(r.x)//64,int(r.y)//64):[r.H_mean,r.E_mean,r.gray_mean,r.sat_mean] for r in patch.itertuples()}
 features=np.array([coords.get((int(x)//64,int(y)//64),[np.nan]*4) for x,y in d[['x_hires_px','y_hires_px']].to_numpy()])
 valid=np.isfinite(features).all(1)&np.isfinite(d[['p_tumor_anchor','p_epithelial_anchor','p_normal_kidney']]).all(1).to_numpy()
 out={}
 if valid.sum()<100:return valid,out,dict(reason='fewer_than_100_valid_spots')
 F=features.copy();F[valid]=(F[valid]-F[valid].mean(0))/np.maximum(F[valid].std(0),1e-8)
 rr,cc=G.nonzero();keep=valid[rr]&valid[cc];rr=rr[keep];cc=cc[keep]
 dist=np.sum((F[rr]-F[cc])**2,axis=1);positive=dist[dist>0];sigma=np.median(positive) if len(positive) else 1.
 weights=np.exp(-dist/(2*max(sigma,1e-8)))
 W=sparse.csr_matrix((weights,(rr,cc)),shape=G.shape);degree=np.asarray(W.sum(1)).ravel();W=sparse.diags(1/np.maximum(degree,1e-12))@W
 contrast=(d.p_tumor_anchor-d.p_normal_kidney).to_numpy()
 seed=valid&(d.p_tumor_anchor>=np.quantile(d.loc[valid,'p_tumor_anchor'],.75))&(d.p_epithelial_anchor>=np.quantile(d.loc[valid,'p_epithelial_anchor'],.5))&d.d_tumor_anchor.ge(2).to_numpy()&(contrast>0)
 epi=d.p_epithelial_anchor>=np.quantile(d.loc[valid,'p_epithelial_anchor'],.25)
 smooth=.5*contrast+.5*(W@np.nan_to_num(contrast))
 graph=sparse.csr_matrix((np.ones(len(rr)),(rr,cc)),shape=G.shape)
 for name,score in [('morphology_supported',smooth),('expression_sensitivity',contrast)]:
  candidate=valid&epi.to_numpy()&(score>0);ii=np.flatnonzero(candidate);inside=np.zeros(len(d),bool)
  if len(ii):
   nc,labels=connected_components(graph[ii][:,ii],directed=False)
   for cl in range(nc):
    ids=ii[labels==cl]
    if len(ids)>=10 and seed[ids].any():inside[ids]=True
  out[name]=inside
 denom=out['morphology_supported'].sum()+out['expression_sensitivity'].sum();dice=2*np.sum(out['morphology_supported']&out['expression_sensitivity'])/denom if denom else 0.
 return valid,out,dict(reason='COMPUTED',seeds=int(seed.sum()),valid_spots=int(valid.sum()),inside_dice=float(dice),H_E_edge_scale=float(sigma))
def main():
 gate('spatial_boundary_extract')
 panels=json.loads((P/'PANELS_FROZEN.json').read_text());reg=read(R/'phase234_part4_replan_execution_20260921_v1/tables/INPUT_SECTION_REGISTER.tsv')
 quals=[];coverage=[];manifest=[];topology=[]
 for r in reg.itertuples():
  d,cov=extract(r,panels);coverage+=cov
  G=sparse.load_npz(R/'phase234_part4_replan_execution_20260921_v1/graphs'/f'{r.sample}_ring.npz').tocsr()
  patch=read(R/'phase287_stromal_refinement_image_morphology_20260929_v1/tables'/f'{r.sample}_IMAGE_PATCHES_64.tsv.gz')
  valid,masks,q=define(d,G,patch);d['boundary_input_valid']=valid
  for key,mask in masks.items():d['inside_'+key]=mask
  d.to_csv(P/'inputs'/f'{r.sample}_SPOT_PROGRAMS.tsv.gz',sep='\t',index=False)
  f=Path(r.input_path);st=f.stat();manifest.append(dict(sample=r.sample,path=str(f),bytes=st.st_size,mtime_ns=st.st_mtime_ns))
  if not masks:
   quals.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,definition='none',eligible=False,**q));continue
  xy=d[['physical_x_um','physical_y_um']].to_numpy();tree=cKDTree(xy)
  gg=sparse.diags(valid.astype(float))@G@sparse.diags(valid.astype(float));nc,components=connected_components(gg,directed=False)
  rr,cc=gg.nonzero()
  for definition,inside in masks.items():
   outside=valid&~inside;take=(rr<cc)&(inside[rr]!=inside[cc]);a=rr[take];b=cc[take]
   inner=np.where(inside[a],a,b);outer=np.where(inside[a],b,a)
   ok=int(inside.sum())>=50 and int(outside.sum())>=50 and len(a)>=20 and len(set(inner))>=10 and len(set(outer))>=10 and q['inside_dice']>=.6 and bool(r.physical_eligible)
   quals.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,definition=definition,inside=int(inside.sum()),outside=int(outside.sum()),crossing_edges=len(a),inner_endpoints=len(set(inner)),outer_endpoints=len(set(outer)),eligible=bool(ok),**q))
   if not len(a):continue
   mid=(xy[a]+xy[b])/2
   macro=(d.total_macro_fraction>=d.total_macro_fraction.quantile(.75))&(d.malignant_fraction<=d.malignant_fraction.quantile(.5))
   rows=[];neighbor_records=[]
   for radius in [200,400]:
    edge=[]
    for k,lst in enumerate(tree.query_ball_point(mid,radius)):
     ids=np.array([i for i in lst if valid[i] and components[i]==components[inner[k]]],int)
     ii=ids[inside[ids]];oo=ids[~inside[ids]];mo=oo[macro.iloc[oo].to_numpy()]
     rec=dict(edge=k,inner_index=int(inner[k]),outer_index=int(outer[k]),x=float(mid[k,0]),y=float(mid[k,1]),inside_spots=len(ii),outside_spots=len(oo),macro_enriched_outside_spots=len(mo),geometry_eligible=bool(ok and len(ii)>=3 and len(oo)>=3))
     for side,inds in [('in',ii),('out',oo),('macro_out',mo)]:
      for col in ['g_TYMP','g_SPP1','p_CAF_activation','p_transition_RCC','p_TGF_response','p_ECM_organization','p_ECM_remodeling','p_source_identity','CD8_expression','Tcell_fraction','log1p_UMI','total_macro_fraction','malignant_fraction','endothelial_fraction','pericyte_fraction']:
       rec[side+'_'+col]=float(d[col].iloc[inds].mean()) if len(inds)>=3 else np.nan
     rec['log_side_count_ratio']=float(np.log((len(ii)+1)/(len(oo)+1)))
     edge.append(rec);neighbor_records.append(dict(definition=definition,radius=radius,edge=k,inside=ii.tolist(),outside=oo.tolist(),macro_out=mo.tolist()))
    ed=pd.DataFrame(edge)
    x=ed.macro_out_g_SPP1.to_numpy();y=ed.out_p_CAF_activation.to_numpy();finite=np.isfinite(x)&np.isfinite(y)
    joint=np.full(len(ed),np.nan)
    if finite.any():joint[finite]=(stats.rankdata(x[finite])/(finite.sum()+1))*(stats.rankdata(y[finite])/(finite.sum()+1))
    ed['joint_rank']=joint
    high=np.zeros(len(ed),bool)
    if finite.any():high[finite]=(x[finite]>=np.quantile(x[finite],.75))&(y[finite]>=np.quantile(y[finite],.75))&(x[finite]>0)&(y[finite]>0)
    ed['joint_high']=high
    # Connected high crossing midpoints, <=101um; descriptive graph continuity only.
    good=finite&ed.geometry_eligible.to_numpy();hs=np.flatnonzero(high&good);longest=0
    if len(hs):
     pairs=list(cKDTree(mid[hs]).query_pairs(101.));ar=np.array(pairs,dtype=int).reshape(-1,2);B=sparse.csr_matrix((np.ones(2*len(ar)),(np.r_[ar[:,0],ar[:,1]],np.r_[ar[:,1],ar[:,0]])),shape=(len(hs),len(hs)))
     _,lab=connected_components(B,directed=False);longest=int(np.bincount(lab).max())
    topology.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,definition=definition,radius=radius,eligible=bool(ok),evaluable_edges=int(good.sum()),joint_high_edges=int((high&good).sum()),joint_coverage_fraction=float((high&good).sum()/good.sum()) if good.sum() else np.nan,largest_connected_high_edge_count=longest,interpretation='graph_continuity_proxy_not_cellular_wall'))
    ed.to_csv(P/'inputs'/f'{r.sample}_{definition}_{radius}_EDGES.tsv.gz',sep='\t',index=False)
   dump(P/'inputs'/f'{r.sample}_{definition}_MEMBERSHIPS.json',neighbor_records)
  print('BOUNDARY_DONE',r.sample,[(k,int(v.sum())) for k,v in masks.items()],flush=True)
 tab(quals,'BOUNDARY_QUALIFICATION.tsv');tab(coverage,'SPATIAL_PROGRAM_COVERAGE.tsv');tab(topology,'BOUNDARY_CONTINUITY_DESCRIPTIVE.tsv')
 dump(P/'qc/SPATIAL_INPUT_MANIFEST.json',manifest);dump(P/'qc/BOUNDARY_EXTRACTION_COMPLETE.json',dict(status='COMPUTED',sections=len(reg),definitions=len(quals),eligible_definitions=sum(r['eligible'] for r in quals),uses_TYMP_SPP1_CAF_CD8_for_boundary=False,actual_CNV_inference=False,new_figures=0))
 print('BOUNDARY_EXTRACTION_COMPLETE',flush=True)
if __name__=='__main__':main()
