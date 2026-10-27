from pathlib import Path
import json,datetime,hashlib,traceback,ast,sys
import numpy as np,pandas as pd,h5py
from scipy import sparse,stats
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import connected_components
from anndata.io import read_elem
P=Path(__file__).resolve().parents[1];R=P.parent;B=R.parent
# Reuse the already frozen statistical functions without executing the previous pipeline.
source=R/'phase222_spatial_S2_S5_execution_20260920_v1/scripts/run_S2_S5.py';tree=ast.parse(source.read_text());funcs=[x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name in ['partial','infer','average_graph']];exec(compile(ast.Module(body=funcs,type_ignores=[]),str(source),'exec'))
def save(n,v):(P/n).write_text(json.dumps(v,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
def tab(n,d):d.to_csv(P/n,sep='\t',index=False)
def run():
 items=json.loads((R/'phase217_spatial_asset_progress_inventory_20260920_v1/LEGACY14_EXPRESSION_ASSETS.json').read_text());records=[];effects=[];aliases=[];coords_audit=[]
 for item in items:
  sid=item['sample'];cohort='GLMF_GSE210041' if item['cohort']=='GLMF' else 'TJRCC';sub='glmf12' if item['cohort']=='GLMF' else 'tjrcc12';base=B/'project26210_part4_tymp_reset_phase17_20260812_v1'/(sub+'_cell2location_tymp_continuous_v1')
  run=json.loads((base/'run_metadata_tymp_continuous_v1.json').read_text());assert run['reference_sha256']=='9f7aabc593422af2af10b7db1ecd728f272a8b7f5b2c550ebfc44ddaa6e744d4'
  rawpath=Path(item['input_h5ad']);post=base/sid/(sid+'_tymp_continuous_v1_cell2location.h5ad');abpath=base/sid/'cell_abundance_q05_tymp_continuous_v1.tsv.gz'
  with h5py.File(rawpath) as h:
   obs=read_elem(h['obs']);var=read_elem(h['var']);X=read_elem(h['X']);xy=read_elem(h['obsm']['spatial']);spatial=h['uns']['spatial'];images=[]
   for lib in spatial:
    if 'images' in spatial[lib]:
     for name in spatial[lib]['images']:images.append({'library':lib,'image':name,'shape':list(spatial[lib]['images'][name].shape)})
  with h5py.File(post) as h:postobs=read_elem(h['obs']);assert postobs.index.equals(obs.index)
  a=pd.read_csv(abpath,sep='\t',index_col=0);assert a.index.equals(obs.index);a.columns=[c.removeprefix('q05cell_abundance_w_sf_') for c in a.columns];assert len(a.columns)==13 and np.isfinite(a.to_numpy()).all() and (a.to_numpy()>=0).all()
  q=a.div(a.sum(axis=1),axis=0);totalmacro=q[['SPP1-state Macrophage','C1QC-state Macrophage','Other Macrophage']].sum(axis=1).to_numpy();share=q['SPP1-state Macrophage'].to_numpy()/totalmacro
  assert np.isfinite(X.data).all() and (X.data>=0).all() and np.allclose(X.data,np.round(X.data));depth=np.log1p(np.asarray(X.sum(axis=1)).ravel());mal=q.Malignant.to_numpy()
  # Independent lattice graph, with physical cutoff to avoid spurious distant edges.
  ar=obs.array_row.to_numpy(int);ac=obs.array_col.to_numpy(int);loc={(int(i),int(j)):k for k,(i,j) in enumerate(zip(ar,ac))};assert len(loc)==len(obs)
  ds=cKDTree(xy).query(xy,k=2)[0][:,1];spacing=float(np.median(ds));rr=[];cc=[]
  for i,(row,col) in enumerate(zip(ar,ac)):
   for dr,dc in [(0,2),(0,-2),(1,1),(1,-1),(-1,1),(-1,-1)]:
    j=loc.get((int(row+dr),int(col+dc)))
    if j is not None and np.linalg.norm(xy[i]-xy[j])<=1.3*spacing+1e-6:rr.append(i);cc.append(j)
  G=sparse.csr_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(obs),len(obs)));W,deg=average_graph(G);assert(G!=G.T).nnz==0 and deg.max()<=6;nc,comp=connected_components(G,directed=False);mask=(totalmacro>=.1)&(deg>=3);nm=W@totalmacro;nd=W@depth
  patient={'GLMF1':'GSE210041:Patient1','GLMF2':'GSE210041:Patient7'}.get(sid,'TJRCC:'+sid.removesuffix('_T'));unit=dict(sample=sid,dataset=cohort,patient=patient)
  frozen=pd.DataFrame(dict(spot_id=obs.index,source_eligible=mask,SPP1_share=share,total_macro_fraction=totalmacro,malignant_fraction=mal,log1p_UMI=depth,neighbor_macro=nm,neighbor_depth=nd,component=comp,x=xy[:,0],y=xy[:,1]))
  for label,key in [('Endothelial','endothelial'),('SMC/Pericyte','pericyte')]:
   v=q[label].to_numpy();y=W@v;cov=np.column_stack([totalmacro,v,nm,mal,depth,nd]);reason=[]
   if (v>=.05).sum()<20:reason.append('TARGET_COVERAGE_LT20')
   if len(np.unique(share[mask]))<10 or len(np.unique(y[mask]))<10:reason.append('INSUFFICIENT_VARIATION')
   eff,info=partial(share,y,cov,mask)
   if info['reason']:reason.append(info['reason'])
   if reason:eff=np.nan
   effects.append(dict(**unit,endpoint='S2_'+key,effect=eff,effect_kind='correlation',eligible=np.isfinite(eff),n_source=int(mask.sum()),n_target=int((v>=.05).sum()),condition=info['condition'],reason=';'.join(reason)))
   frozen[key+'_fraction']=v;frozen['neighbor_'+key]=y;np.savez_compressed(P/'inputs'/(sid+'_'+key+'_DESIGN.npz'),x=share,y=y,cov=cov,mask=mask)
  tumor=mal>=.2;tn=G@tumor.astype(float);core=mask&tumor&(tn==deg);boundary=mask&((tumor&(tn<deg))|((~tumor)&(tn>0)));ok=core.sum()>=20 and boundary.sum()>=20
  effects.append(dict(**unit,endpoint='S3_expression_boundary',effect=float(np.median(share[boundary])-np.median(share[core])) if ok else np.nan,effect_kind='difference',eligible=ok,n_core=int(core.sum()),n_interface=int(boundary.sum()),reason='' if ok else 'REGION_SOURCE_LT20'))
  gi=var.index.get_indexer(['COL1A1','COL1A2','DCN','LUM']);lo=hi=np.zeros(len(obs),bool);ok=False
  if (gi>=0).all():
   countsum=np.asarray(X.sum(axis=1)).ravel();ecm=np.log1p(X[:,gi].toarray()/countsum[:,None]*1e4).mean(axis=1);q1,q3=np.quantile(ecm,[.25,.75]);lo=mask&(ecm<=q1);hi=mask&(ecm>=q3);ok=q3>q1 and lo.sum()>=20 and hi.sum()>=20
  effects.append(dict(**unit,endpoint='S3_ECM_region',effect=float(np.median(share[hi])-np.median(share[lo])) if ok else np.nan,effect_kind='difference',eligible=ok,n_low=int(lo.sum()),n_high=int(hi.sum()),reason='' if ok else 'ECM_GENE_OR_REGION_COVERAGE'))
  frozen['core_source']=core;frozen['boundary_source']=boundary;frozen['ECM_low']=lo;frozen['ECM_high']=hi;tab('inputs/'+sid+'_FROZEN_TRANSFER_INPUTS.tsv.gz',frozen);sparse.save_npz(P/'inputs'/(sid+'_adjacency.npz'),G)
  rec=dict(**unit,spots=len(obs),genes=len(var),mapping_genes=run.get('gene_set',run.get('reference_shape')),images=images,components=int(nc),source_spots=int(mask.sum()),reference_signature_sha256=run['reference_sha256'],N_cells_per_location=run['parameters']['N_cells_per_location'],input_path=str(rawpath),posterior_path=str(post),abundance_path=str(abpath),TLS_endpoint='NOT_RUN_NO_AUDITED_TLS_ANNOTATION',role='MAIN_WHOLE_SECTION_CASE' if cohort=='GLMF_GSE210041' else 'EXTERNAL_COHORT')
  if cohort=='GLMF_GSE210041':
   accession='GSM6415705' if sid=='GLMF1' else 'GSM6415706';rp=Path('/mnt/BioDisk/ludingyu/project26210/external/GSE210041/unpacked')/(accession+'_'+sid+'_filtered_feature_bc_matrix.h5')
   with h5py.File(rp) as h:
    mat=h['matrix'];shape=mat['shape'][:];raw=sparse.csr_matrix((mat['data'][:],mat['indices'][:],mat['indptr'][:]),shape=(int(shape[1]),int(shape[0])));bars=mat['barcodes'].asstr()[:];names=mat['features']['name'].asstr()[:]
   # Gene order is retained in these17,943-gene QC inputs; do not use name-only approximations.
   assert len(names)==len(var) and np.array_equal(names,var.index.astype(str).to_numpy()),'raw/QC gene order differs'
   ri=pd.Index(bars).get_indexer(obs.index);assert(ri>=0).all();diff=raw[ri].astype(float)-X.astype(float);diff.eliminate_zeros();assert diff.nnz==0
   rec.update(GEO_sample=accession,raw_source=str(rp),raw_to_QC_counts_exact=True,original_filtered_spots=len(bars))
   rawsp=rp.parent/sid/'spatial';scale=json.loads((rawsp/'scalefactors_json.json').read_text());factor=scale['tissue_hires_scalef'];image=[im for im in images if im['image']=='hires'];assert len(image)==1;height,width=image[0]['shape'][:2];points=xy*factor;inside=(points[:,0]>=0)&(points[:,0]<width)&(points[:,1]>=0)&(points[:,1]<height);assert inside.all()
   aliases.append(dict(sample=sid,dataset='GSE210041',GEO_sample=accession,patient=patient,counts_exact=True,qc_spots=len(obs),raw_filtered_spots=len(bars),raw_file=str(rp),raw_sha256=hashlib.sha256(rp.read_bytes()).hexdigest(),prior_treatment='Nivolumab+Ipilimumab' if sid=='GLMF1' else 'Naive'))
   coords_audit.append(dict(sample=sid,hires_shape=image[0]['shape'],all_QC_spots_inside_image=bool(inside.all()),scalefactor=factor,pathology_status='AI_DRAFT_NOT_PATHOLOGIST_REVIEWED',full_cut_claim='complete deposited Visium tissue/image support; not proof of whole-organ coverage'))
  records.append(rec);print('TRANSFERRED',sid,flush=True)
 ef=pd.DataFrame(effects);tab('results/ALL_14_SECTION_EFFECTS.tsv',ef);save('ASSET_ROLE_MANIFEST.json',records);save('GSE210041_GLMF_ALIAS_PROOF.json',aliases);save('GLMF_IMAGE_COORDINATE_CHECK.json',coords_audit)
 summaries=[]
 for j,((dataset,endpoint),g) in enumerate(ef.groupby(['dataset','endpoint'])):
  values=g.effect.dropna().to_numpy();row=infer(values,262223+j);row.update(dataset=dataset,endpoint=endpoint,patients=','.join(g.loc[g.effect.notna(),'patient']),FDR=np.nan)
  if dataset=='GLMF_GSE210041':row.update(P=np.nan,CI_low=np.nan,CI_high=np.nan)
  summaries.append(row)
 ss=pd.DataFrame(summaries)
 for dataset,g in ss[ss.endpoint.str.startswith('S2_')].groupby('dataset'):
  ix=g.index;pp=g.P.fillna(1).to_numpy();q=stats.false_discovery_control(pp);ss.loc[ix,'FDR']=np.where(g.P.notna(),q,np.nan)
 tab('results/COHORT_SPECIFIC_SUMMARY.tsv',ss)
 old=json.loads((P/'audit/RDS_BEFORE.json').read_text())
 for a in old:
  st=Path(a['path']).stat();assert(st.st_size,st.st_mtime_ns,st.st_ino)==(a['bytes'],a['mtime_ns'],a['inode'])
 save('audit/RDS_PRESERVATION.json',dict(tracked=len(old),changed=0));save('STATUS.json',dict(status='TRANSFER_COMPUTE_COMPLETE',added_sections=14,GLMF_equals_GSE210041=True,new_dataset_downloads=0,new_maps=0,new_figures=0));(P/'EXIT_STATUS').write_text('0\n')
try:run()
except Exception:
 (P/'ERROR.txt').write_text(traceback.format_exc());(P/'EXIT_STATUS').write_text('1\n');raise
