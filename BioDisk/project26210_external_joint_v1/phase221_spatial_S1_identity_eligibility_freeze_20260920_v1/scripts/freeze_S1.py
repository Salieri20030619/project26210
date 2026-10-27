from pathlib import Path
import json,hashlib,datetime,traceback,shutil,platform
import numpy as np,pandas as pd,h5py,scipy
from scipy import sparse,stats
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from anndata.io import read_elem
P=Path(__file__).resolve().parents[1];R=P.parent;B=R.parent
rules=json.loads((P/'S1_RULES_BEFORE_ELIGIBILITY.json').read_text());sources=[]
def record(f):
 f=Path(f);sources.append({'path':str(f),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest() if f.stat().st_size<30_000_000 else None})
def write(name,obj):(P/name).write_text(json.dumps(obj,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
def tab(name,df):df.to_csv(P/name,sep='\t',index=False)
def run():
 mp=R/'phase168_reference_source_audit_20260916_v1/SPATIAL_26_SECTION_PATIENT_MAP.tsv';record(mp);pm=pd.read_csv(mp,sep='\t');pm=pm[pm.dataset=='GSE175540'].copy().set_index('section');assert len(pm)==24 and pm.patient_id.nunique()==23
 sec54=B/'project26210_external_joint_v1_phase54_gse175540_spatial_v1/run_v1/sections_v1';sec78=B/'project26210_phase78_fine_state_spatial_v1/run_v1/sections_v1';extract=sec54.parents[1]/'input_extract_v1';old=R/'phase167_A_24section_endpoint_check_20260916_v1'
 eligibility=[];technical=[];optional=[];checks=[]
 for sid,m in pm.iterrows():
  hpath=sec54/sid/'processed_section_v1.h5ad';apath=sec78/sid/'fine_state_abundance_v1.tsv.gz';gpath=old/(sid+'_adjacency.npz');tpath=old/(sid+'_SPOT_TECHNICAL_QC.tsv.gz')
  for f in [hpath,apath,gpath,tpath]:record(f)
  with h5py.File(hpath) as h:obs=read_elem(h['obs']);xy=read_elem(h['obsm']['spatial']);genes=set(read_elem(h['var']).index.astype(str))
  a=pd.read_csv(apath,sep='\t',index_col='spot_id');t=pd.read_csv(tpath,sep='\t',index_col='spot_id');assert a.index.equals(obs.index) and t.index.equals(obs.index) and obs.index.is_unique
  cols=[x for x in a.columns if x.startswith('q05__')];assert len(cols)==13
  q=a[cols].to_numpy(float);assert np.isfinite(q).all() and (q>=0).all() and (q.sum(axis=1)>0).all();f=q/q.sum(axis=1)[:,None];lookup={x.removeprefix('q05__'):j for j,x in enumerate(cols)}
  get=lambda name:f[:,lookup[name]]
  totalmacro=get('SPP1-state Macrophage')+get('C1QC-state Macrophage')+get('Other Macrophage');share=get('SPP1-state Macrophage')/totalmacro
  G=sparse.load_npz(gpath).tocsr();deg=np.asarray(G.sum(axis=1)).ravel();assert G.shape==(len(obs),len(obs)) and (G!=G.T).nnz==0 and not G.diagonal().any() and np.all(G.data==1) and deg.max()<=6
  assert np.array_equal(deg,t.degree.to_numpy()) and np.isfinite(xy).all() and len(np.unique(xy,axis=0))==len(xy)
  nc,comp=connected_components(G,directed=False);nearest=cKDTree(xy).query(xy,k=2)[0][:,1];spacing=np.median(nearest);gr,gc=G.nonzero();dist=np.linalg.norm(xy[gr]-xy[gc],axis=1);assert (dist<=1.3*spacing+1e-6).all()
  # Existing first-ring graph independently checked against Visium lattice edges.
  dr=np.abs(obs.array_row.to_numpy()[gr]-obs.array_row.to_numpy()[gc]);dc=np.abs(obs.array_col.to_numpy()[gr]-obs.array_col.to_numpy()[gc]);assert (((dr==0)&(dc==2))|((dr==1)&(dc==1))).all()
  W=sparse.diags(np.divide(1.,deg,out=np.zeros_like(deg),where=deg>0))@G
  src=(totalmacro>=.1)&(deg>=3);depth=np.log1p(obs.total_counts.to_numpy(float));nmacro=np.asarray(W@totalmacro).ravel();ndepth=np.asarray(W@depth).ravel()
  spot=pd.DataFrame({'spot_id':obs.index,'component':comp,'degree':deg,'x_pixels':xy[:,0],'y_pixels':xy[:,1],'source_eligible':src,'SPP1_share_within_mature_macro':share,'total_macro_fraction':totalmacro,'malignant_fraction':get('Malignant'),'log1p_UMI':depth,'neighbor_macro_fraction':nmacro,'neighbor_log1p_UMI':ndepth})
  for target in ['Endothelial','SMC/Pericyte']:
   key='endothelial' if target=='Endothelial' else 'pericyte';v=get(target);y=np.asarray(W@v).ravel();spot[key+'_fraction']=v;spot['neighbor_'+key+'_fraction']=y
   cov=np.column_stack([totalmacro,v,nmacro,get('Malignant'),depth,ndepth])[src];reasons=[];rank=-1;condition=None;vx=None;vy=None
   if src.sum()<50:reasons.append('FEWER_THAN_50_MACRO_SOURCE_SPOTS')
   if (v>=.05).sum()<20:reasons.append('FEWER_THAN_20_TARGET_SUPPORTED_SPOTS')
   if len(np.unique(share[src]))<10:reasons.append('SOURCE_SHARE_TOO_FEW_DISTINCT_VALUES')
   if len(np.unique(y[src]))<10:reasons.append('NEIGHBOR_TARGET_TOO_FEW_DISTINCT_VALUES')
   if src.sum()>=10:
    rk=np.column_stack([stats.rankdata(cov[:,j]) for j in range(cov.shape[1])]);sd=rk.std(axis=0);std=(rk-rk.mean(axis=0))/np.where(sd==0,1,sd);D=np.column_stack([np.ones(src.sum()),std]);rank=int(np.linalg.matrix_rank(D));condition=float(np.linalg.cond(D))
    if rank<D.shape[1]:reasons.append('COVARIATE_DESIGN_RANK_DEFICIENT')
    if not np.isfinite(condition) or condition>=30:reasons.append('DESIGN_CONDITION_GE30')
    rx=stats.rankdata(share[src]);ry=stats.rankdata(y[src]);rx=rx-D@np.linalg.lstsq(D,rx,rcond=None)[0];ry=ry-D@np.linalg.lstsq(D,ry,rcond=None)[0];vx=float(np.var(rx));vy=float(np.var(ry))
    if not np.isfinite([vx,vy]).all() or min(vx,vy)<=1e-12:reasons.append('NO_RESIDUAL_VARIATION')
    # Deliberately no covariance, correlation or effect-size computation in S1.
   eligibility.append(dict(section=sid,patient=m.patient_id,preparation=m.preparation,cohort=m.cohort,endpoint='S2_'+key,source_spots=int(src.sum()),target_supported_spots=int((v>=.05).sum()),source_with_target_neighbor=int((src&(np.asarray(G@(v>=.05).astype(float)).ravel()>0)).sum()),design_rank=rank,design_columns=7,condition=condition,source_residual_variance=vx,target_residual_variance=vy,eligible=not reasons,exclusion_reason=';'.join(reasons)))
  tab('frozen_inputs/'+sid+'_SPOT_INPUTS.tsv.gz',spot);shutil.copy2(gpath,P/'frozen_inputs'/(sid+'_adjacency.npz'))
  tlsfiles=list(extract.glob(sid+'*TLS_annotation.csv.gz'));kind='MISSING_FILE';positive=negative=unknown=0;complete=False
  if tlsfiles:
   assert len(tlsfiles)==1;record(tlsfiles[0]);tt=pd.read_csv(tlsfiles[0],dtype=str).set_index('Barcode');assert tt.index.is_unique;complete=obs.index.isin(tt.index).all()
   if 'TLS_2_cat' in tt:
    kind='TLS_VS_NO_TLS';labels=tt['TLS_2_cat'].reindex(obs.index);assert set(labels.dropna().unique())<= {'TLS','NO_TLS'};positive=int((labels=='TLS').sum());negative=int((labels=='NO_TLS').sum());unknown=len(obs)-positive-negative
   elif 'TLS' in tt:
    kind='T_AGGREGATE_ONLY';labels=tt['TLS'].reindex(obs.index);assert set(labels.dropna().unique())<= {'T_agg'};positive=int((labels=='T_agg').sum());negative=0;unknown=len(obs)-positive
   else:kind='UNRECOGNIZED_SCHEMA'
  optional.append(dict(section=sid,patient=m.patient_id,preparation=m.preparation,S3_malignant_identity_present=True,ECM_genes_present=';'.join(sorted(genes&{'COL1A1','COL1A2','DCN','LUM'})),S3_histologic_boundary_frozen=False,S3_fibroblast_identity_available=False,TLS_annotation_type=kind,annotated_positive_spots=positive,explicit_negative_spots=negative,unknown_spots=unknown if tlsfiles else len(obs),all_QC_barcodes_in_annotation=complete,S4_TLS_contrast_candidate=(kind=='TLS_VS_NO_TLS' and positive>=20 and negative>=50 and complete)))
  technical.append(dict(section=sid,patient=m.patient_id,preparation=m.preparation,spots=len(obs),components=nc,isolated_spots=int((deg==0).sum()),source_spots=int(src.sum()),spacing_pixels=float(spacing),units='section_spot_spacing',barcodes_exact=True,graph_lattice_verified=True,q05_valid=True))
  print(sid,'S1 checked',flush=True)
 e=pd.DataFrame(eligibility);o=pd.DataFrame(optional);tech=pd.DataFrame(technical);tab('S2_SECTION_ELIGIBILITY.tsv',e);tab('S3_S4_AVAILABILITY.tsv',o);tab('SECTION_GRAPH_IDENTITY_QC.tsv',tech)
 patient=[]
 for (endpoint,pat),g in e.groupby(['endpoint','patient']):
  ok=g[g.eligible];patient.append(dict(endpoint=endpoint,patient=pat,eligible=bool(len(ok)),eligible_sections=';'.join(ok.section),n_eligible_sections=len(ok),all_sections=';'.join(g.section),aggregation_rule='Fisher_z_equal_section_mean_if_two_eligible;single_section_if_only_one;NA_if_none'))
 tab('S2_PATIENT_ELIGIBILITY.tsv',pd.DataFrame(patient));tab('FROZEN_SECTION_PATIENT_MAP.tsv',pm.reset_index())
 write('INPUT_MANIFEST.json',sources)
 summary={'status':'S1_COMPLETE','sections':24,'patients':23,'spots':int(tech.spots.sum()),'vascular_endpoint_eligibility':{},'TLS_types':o.TLS_annotation_type.value_counts().to_dict(),'TLS_contrast_candidate_sections':int(o.S4_TLS_contrast_candidate.sum()),'TLS_candidate_patients':int(o.loc[o.S4_TLS_contrast_candidate,'patient'].nunique()),'S3_boundary_ready':False,'S4_distance_endpoint_frozen':False,'association_computed':False,'new_mapping_run':False}
 for endpoint,g in e.groupby('endpoint'):
  ok=g[g.eligible];summary['vascular_endpoint_eligibility'][endpoint]={'sections':len(ok),'patients':ok.patient.nunique(),'by_preparation':ok.preparation.value_counts().to_dict(),'excluded':g.loc[~g.eligible,['section','exclusion_reason']].to_dict('records')}
 before=json.loads((P/'audit/RDS_BEFORE.json').read_text());changed=[]
 for x in before:
  st=Path(x['path']).stat()
  if (st.st_size,st.st_mtime_ns,st.st_ino)!=(x['bytes'],x['mtime_ns'],x['inode']):changed.append(x['path'])
 assert not changed
 write('audit/RDS_PRESERVATION.json',dict(tracked=len(before),changed=changed));write('S1_SUMMARY.json',summary)
 frozen={str(f.relative_to(P)):hashlib.sha256(f.read_bytes()).hexdigest() for f in P.rglob('*') if f.is_file() and f.suffix in ['.tsv','.gz','.npz']};write('FREEZE_SHA256.json',frozen)
 write('SOFTWARE.json',dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,pandas=pd.__version__,h5py=h5py.__version__))
 (P/'EXIT_STATUS').write_text('0\n');print(json.dumps(summary),flush=True)
try:run()
except Exception:
 (P/'ERROR.txt').write_text(traceback.format_exc());(P/'EXIT_STATUS').write_text('1\n');raise
