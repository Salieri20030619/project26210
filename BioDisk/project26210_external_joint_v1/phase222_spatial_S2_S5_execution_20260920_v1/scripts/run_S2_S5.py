from pathlib import Path
import json,datetime,traceback,hashlib,sys,os,platform
import numpy as np,pandas as pd,h5py,scipy
from scipy import stats,sparse
from scipy.sparse.csgraph import dijkstra
from anndata.io import read_elem
P=Path(__file__).resolve().parents[1];R=P.parent;B=R.parent;F=R/'phase221_spatial_S1_identity_eligibility_freeze_20260920_v1'
def log(s):print(datetime.datetime.now().isoformat(),s,flush=True)
def save(n,v):(P/n).write_text(json.dumps(v,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
def tab(n,v):v.to_csv(P/n,sep='\t',index=False)
def partial(x,y,cov,mask):
 n=int(mask.sum());reason=[]
 if n<50:return np.nan,dict(n=n,reason='FEWER_THAN_50_SOURCE',rank=0,condition=np.nan)
 a=np.column_stack([x[mask],y[mask],cov[mask]])
 if not np.isfinite(a).all():return np.nan,dict(n=n,reason='NONFINITE_INPUT',rank=0,condition=np.nan)
 rr=np.column_stack([stats.rankdata(a[:,j]) for j in range(a.shape[1])]);c=rr[:,2:];sd=c.std(axis=0);c=(c-c.mean(axis=0))/np.where(sd==0,1,sd);D=np.column_stack([np.ones(n),c]);rank=np.linalg.matrix_rank(D);condition=np.linalg.cond(D)
 if rank<D.shape[1] or not np.isfinite(condition) or condition>=30:return np.nan,dict(n=n,reason='RANK_OR_CONDITION',rank=int(rank),condition=float(condition))
 E=rr[:,:2]-D@np.linalg.lstsq(D,rr[:,:2],rcond=None)[0]
 if min(E.var(axis=0))<=1e-12:return np.nan,dict(n=n,reason='ZERO_RESIDUAL_VARIANCE',rank=int(rank),condition=float(condition))
 return float(np.corrcoef(E.T)[0,1]),dict(n=n,reason='',rank=int(rank),condition=float(condition))
def average_graph(G):
 deg=np.asarray(G.sum(axis=1)).ravel().astype(float)
 return sparse.diags(np.divide(1.,deg,out=np.zeros_like(deg),where=deg>0))@G,deg

def summarize_patients(sec):
 rows=[]
 for keys,g in sec.groupby(['endpoint','model','patient'],dropna=False):
  a=g.effect.dropna().to_numpy();kind=g.effect_kind.iloc[0]
  value=float(np.tanh(np.arctanh(np.clip(a,-.999999,.999999)).mean())) if len(a) and kind=='correlation' else (float(a.mean()) if len(a) else np.nan)
  rows.append(dict(endpoint=keys[0],model=keys[1],patient=keys[2],effect=value,effect_kind=kind,n_eligible_sections=len(a),sections=';'.join(g.loc[g.effect.notna(),'section']),cohort=g.cohort.iloc[0]))
 return pd.DataFrame(rows)
def infer(a,seed):
 a=np.asarray(a,float);a=a[np.isfinite(a)];n=len(a);pos=int((a>0).sum());neg=int((a<0).sum());result=dict(n_patients=n,positive=pos,negative=neg,zero=n-pos-neg,median_effect=float(np.median(a)) if n else np.nan,P=np.nan,CI_low=np.nan,CI_high=np.nan)
 if n>=10:
  result['P']=float(stats.binomtest(pos,pos+neg,.5).pvalue) if pos+neg else 1.
  rng=np.random.default_rng(seed);med=np.median(a[rng.integers(0,n,size=(5000,n))],axis=1);result['CI_low'],result['CI_high']=np.quantile(med,[.025,.975]).tolist()
 return result

def run():
 expected=json.loads((F/'FREEZE_SHA256.json').read_text());assert all(hashlib.sha256((F/k).read_bytes()).hexdigest()==v for k,v in expected.items())
 meta=pd.read_csv(F/'FROZEN_SECTION_PATIENT_MAP.tsv',sep='\t').set_index('section');optional=pd.read_csv(F/'S3_S4_AVAILABILITY.tsv',sep='\t').set_index('section')
 root54=B/'project26210_external_joint_v1_phase54_gse175540_spatial_v1/run_v1/sections_v1';root78=B/'project26210_phase78_fine_state_spatial_v1/run_v1/sections_v1';extract=root54.parents[1]/'input_extract_v1'
 shifts=[(dq,dr) for dq in range(-8,9) for dr in range(-8,9) if 3<=max(abs(dq),abs(dr),abs(dq+dr))<=8];rng=np.random.default_rng(262222);shifts=[shifts[i] for i in rng.permutation(len(shifts))[:99]];save('S2/FROZEN_TRANSLATION_VECTORS.json',shifts)
 sections=[];translation=[];regionqc=[];tlsc=[];spotassign=[]
 for sid,m in meta.iterrows():
  z=pd.read_csv(F/'frozen_inputs'/(sid+'_SPOT_INPUTS.tsv.gz'),sep='\t',index_col='spot_id');G=sparse.load_npz(F/'frozen_inputs'/(sid+'_adjacency.npz')).tocsr();W,deg=average_graph(G);mask=z.source_eligible.to_numpy(bool);n=len(z);comp=z.component.to_numpy();depth=z.log1p_UMI.to_numpy();ndepth=z.neighbor_log1p_UMI.to_numpy();source=z.SPP1_share_within_mature_macro.to_numpy();macro=z.total_macro_fraction.to_numpy();nm=z.neighbor_macro_fraction.to_numpy();mal=z.malignant_fraction.to_numpy()
  with h5py.File(root54/sid/'processed_section_v1.h5ad') as h:
   obs=read_elem(h['obs']);var=read_elem(h['var']);X=read_elem(h['X']);assert obs.index.equals(z.index)
  ab=pd.read_csv(root78/sid/'fine_state_abundance_v1.tsv.gz',sep='\t',index_col='spot_id');me=ab.filter(regex='^mean__').copy();me=me.div(me.sum(axis=1),axis=0);mm=(me['mean__SPP1-state Macrophage']+me['mean__C1QC-state Macrophage']+me['mean__Other Macrophage']).to_numpy();sxmean=me['mean__SPP1-state Macrophage'].to_numpy()/mm
  # Exact graph distance2, excluding distance0 and distance1.
  G2=(G@G).tocsr();G2.data[:]=1;G2.setdiag(0);G2=G2-G2.multiply(G);G2.eliminate_zeros();W2,deg2=average_graph(G2)
  xy=z[['x_pixels','y_pixels']].to_numpy();xy=(xy-xy.mean(axis=0))/xy.std(axis=0);coords=np.column_stack([xy,xy[:,0]**2,xy[:,0]*xy[:,1],xy[:,1]**2]);present=np.unique(comp[mask]);dummies=np.column_stack([comp==c for c in present[1:]]) if len(present)>1 else np.empty((n,0));coord_cov=np.column_stack([coords,dummies])
  largest=int(np.argmax(np.bincount(comp)));base=dict(section=sid,patient=m.patient_id,preparation=m.preparation,cohort=m.cohort)
  def add(endpoint,model,value,kind,info):sections.append(dict(**base,endpoint=endpoint,model=model,effect=value,effect_kind=kind,**info))
  # Build translations independently of expression outcomes.
  ar=obs.array_row.to_numpy(int);ac=obs.array_col.to_numpy(int);axq=(ac-ar)//2;axr=ar;loc={(int(q),int(r)):i for i,(q,r) in enumerate(zip(axq,axr))};translations=[]
  for dq,dr in shifts:
   ix=np.array([loc.get((int(q+dq),int(r+dr)),-1) for q,r in zip(axq,axr)]);valid=ix>=0;valid &= comp[np.maximum(ix,0)]==comp
   valid_neighbours=np.asarray(G@valid.astype(float)).ravel();matched=mask&valid&(valid_neighbours==deg)
   if matched.sum()>=max(50,int(np.ceil(mask.sum()*.5))):translations.append((dq,dr,ix,valid,matched))
  for label,key in [('Endothelial','endothelial'),('SMC/Pericyte','pericyte')]:
   v=z[key+'_fraction'].to_numpy();y=z['neighbor_'+key+'_fraction'].to_numpy();cov=np.column_stack([macro,v,nm,mal,depth,ndepth]);effect,info=partial(source,y,cov,mask);assert np.isfinite(effect);add('S2_'+key,'primary',effect,'correlation',info)
   np.savez_compressed(P/'cache'/(sid+'_'+key+'_PRIMARY.npz'),x=source,y=y,cov=cov,mask=mask)
   vm=me['mean__'+label].to_numpy();ym=np.asarray(W@vm).ravel();mc=np.column_stack([mm,vm,np.asarray(W@mm).ravel(),me['mean__Malignant'].to_numpy(),depth,ndepth])
   for model,xx,yy,cv,ma in [('posterior_mean',sxmean,ym,mc,mask),('coordinate_component',source,y,np.column_stack([cov,coord_cov]),mask),('largest_component',source,y,cov,mask&(comp==largest)),('second_ring',source,np.asarray(W2@v).ravel(),cov,mask&(deg2>0))]:
    val,inf=partial(xx,yy,cv,ma);add('S2_'+key,model,val,'correlation',inf)
   local=[]
   for dq,dr,ix,valid,matched in translations:
    shifted=np.zeros(n);shifted[valid]=v[ix[valid]];shift_y=np.asarray(W@shifted).ravel();covshift=cov.copy();covshift[:,1]=shifted
    observed,io=partial(source,y,cov,matched);null,inull=partial(source,shift_y,covshift,matched)
    if np.isfinite(observed) and np.isfinite(null):local.append(dict(**base,endpoint='S2_'+key,dq=dq,dr=dr,n_matched=int(matched.sum()),matched_observed=observed,shifted=null,excess=observed-null))
   translation.extend(local)
   val=float(np.median([a['excess'] for a in local])) if len(local)>=30 else np.nan;add('S2_'+key,'local_translation_excess',val,'difference',dict(n=int(mask.sum()),reason='' if len(local)>=30 else 'FEWER_THAN_30_VALID_TRANSLATIONS',valid_shifts=len(local)))
  # Computational tumour-region boundary, not histology.
  tumor=mal>=.2;tn=np.asarray(G@tumor.astype(float)).ravel();core=tumor&(tn==deg);interface=(tumor&(tn<deg))|((~tumor)&(tn>0));ca=mask&core;ia=mask&interface
  good=ca.sum()>=20 and ia.sum()>=20;val=float(np.median(source[ia])-np.median(source[ca])) if good else np.nan;add('S3_expression_boundary','primary',val,'difference',dict(n=int(ca.sum()+ia.sum()),n_core=int(ca.sum()),n_interface=int(ia.sum()),reason='' if good else 'REGION_SOURCE_COVERAGE_LT20'))
  ecmgenes=['COL1A1','COL1A2','DCN','LUM'];idx=var.index.get_indexer(ecmgenes);assert (idx>=0).all();ecm=X[:,idx].toarray().astype(float).mean(axis=1);q1,q3=np.quantile(ecm,[.25,.75]);lo=mask&(ecm<=q1);hi=mask&(ecm>=q3);good=q3>q1 and lo.sum()>=20 and hi.sum()>=20
  val=float(np.median(source[hi])-np.median(source[lo])) if good else np.nan;add('S3_ECM_region','primary',val,'difference',dict(n=int(lo.sum()+hi.sum()),n_low=int(lo.sum()),n_high=int(hi.sum()),reason='' if good else 'ECM_VARIATION_OR_REGION_SOURCE_COVERAGE'))
  regionqc.append(dict(**base,n_tumour_spots=int(tumor.sum()),n_core_source=int(ca.sum()),n_interface_source=int(ia.sum()),ECM_q25=q1,ECM_q75=q3,n_ECM_high_source=int(hi.sum()),n_ECM_low_source=int(lo.sum())))
  # Supplemental annotation-constrained TLS geometry.
  near=np.zeros(n,bool);far=np.zeros(n,bool);dd=np.full(n,np.inf);TLS=np.zeros(n,bool);NO=np.zeros(n,bool);tlsreason='NOT_IN_S1_FULLY_MATCHED_TLS_CANDIDATES'
  if bool(optional.loc[sid,'S4_TLS_contrast_candidate']):
   path=list(extract.glob(sid+'*TLS_annotation.csv.gz'));assert len(path)==1;tt=pd.read_csv(path[0],dtype=str).set_index('Barcode');ls=tt.TLS_2_cat.reindex(z.index);TLS=(ls=='TLS').to_numpy();NO=(ls=='NO_TLS').to_numpy();known=TLS|NO;Kg=(sparse.diags(known.astype(float))@G@sparse.diags(known.astype(float))).tocsr();Kg.eliminate_zeros()
   dd=dijkstra(Kg,directed=False,indices=np.flatnonzero(TLS),min_only=True,unweighted=True);near=mask&NO&(dd>=1)&(dd<=3);far=mask&NO&np.isfinite(dd)&(dd>=6);tlsreason='' if near.sum()>=20 and far.sum()>=20 else 'TLS_NEAR_OR_FAR_SOURCE_LT20'
  val=float(np.median(source[near])-np.median(source[far])) if not tlsreason else np.nan;add('S4_TLS_periphery','primary',val,'difference',dict(n=int(near.sum()+far.sum()),n_peri=int(near.sum()),n_far=int(far.sum()),reason=tlsreason))
  tlsc.append(dict(**base,S1_candidate=bool(optional.loc[sid,'S4_TLS_contrast_candidate']),TLS_source=int((mask&TLS).sum()),peri_source=int(near.sum()),far_source=int(far.sum()),unreachable_or_unknown_source=int((mask&(~np.isfinite(dd))).sum()),eligible=not tlsreason,reason=tlsreason))
  spot=pd.DataFrame(dict(spot_id=z.index,tumor_expression_region=tumor,core_source=ca,interface_source=ia,ECM_score=ecm,ECM_low_source=lo,ECM_high_source=hi,TLS_positive=TLS,TLS_explicit_negative=NO,TLS_graph_distance_hops=dd,TLS_peri_source=near,TLS_far_source=far));tab('cache/'+sid+'_REGION_MEMBERS.tsv.gz',spot)
  log('S2/S3/S4 complete '+sid)
 sec=pd.DataFrame(sections);pat=summarize_patients(sec);tab('ALL_SECTION_EFFECTS.tsv',sec);tab('ALL_PATIENT_EFFECTS.tsv',pat);tab('S2/TRANSLATION_CONTROLS.tsv.gz',pd.DataFrame(translation));tab('S3/REGION_ELIGIBILITY.tsv',pd.DataFrame(regionqc));tab('S4/TLS_ELIGIBILITY.tsv',pd.DataFrame(tlsc))
 summaries=[]
 for j,((ep,model),g) in enumerate(pat.groupby(['endpoint','model'])):summaries.append(dict(endpoint=ep,model=model,**infer(g.effect.to_numpy(),262222+j)))
 sm=pd.DataFrame(summaries);sm['FDR']=np.nan;sm['family']='SINGLE_PRESET_ENDPOINT'
 for model,g in sm[sm.endpoint.str.startswith('S2')].groupby('model'):
  ii=g.index;assert len(ii)==2;ps=sm.loc[ii,'P'].fillna(1).to_numpy();qs=stats.false_discovery_control(ps);sm.loc[ii,'FDR']=np.where(sm.loc[ii,'P'].notna(),qs,np.nan);sm.loc[ii,'family']='VASCULAR_TWO_ENDPOINTS__'+model
 tab('PATIENT_RESULT_SUMMARY.tsv',sm)
 # Fixed S5 summaries: no spot or repeated-section inferential sample inflation.
 strata=[];loo=[]
 for ep in sec.endpoint.unique():
  sp=sec[(sec.endpoint==ep)&(sec.model=='primary')];pp=pat[(pat.endpoint==ep)&(pat.model=='primary')]
  for column in ['preparation','cohort']:
   for name,g in sp.groupby(column):
    zz=summarize_patients(g);result=infer(zz.effect.to_numpy(),262223);result['P']=np.nan;result['CI_low']=np.nan;result['CI_high']=np.nan;strata.append(dict(endpoint=ep,stratum_type=column,stratum=name,**result))
  for patient in pp.patient:
   a=pp.loc[pp.patient!=patient,'effect'].dropna().to_numpy();loo.append(dict(endpoint=ep,drop_type='patient',dropped=patient,n=len(a),median=float(np.median(a)) if len(a) else np.nan,positive=int((a>0).sum()),negative=int((a<0).sum())))
  for cohort in pp.cohort.unique():
   a=pp.loc[pp.cohort!=cohort,'effect'].dropna().to_numpy();loo.append(dict(endpoint=ep,drop_type='cohort',dropped=cohort,n=len(a),median=float(np.median(a)) if len(a) else np.nan,positive=int((a>0).sum()),negative=int((a<0).sum())))
 tab('S5/PREPARATION_COHORT_DIRECTIONS.tsv',pd.DataFrame(strata));tab('S5/LEAVE_ONE_PATIENT_COHORT.tsv',pd.DataFrame(loo))
 assert all(hashlib.sha256((F/k).read_bytes()).hexdigest()==v for k,v in expected.items())
 old=json.loads((P/'audit/RDS_BEFORE.json').read_text())
 for a in old:
  st=Path(a['path']).stat();assert(st.st_size,st.st_mtime_ns,st.st_ino)==(a['bytes'],a['mtime_ns'],a['inode'])
 save('audit/RDS_PRESERVATION.json',dict(tracked=len(old),changed=0))
 save('SOFTWARE.json',dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,h5py=h5py.__version__))
 save('STATUS.json',dict(status='COMPUTE_COMPLETE',S2='COMPLETE',S3='COMPLETE_EXPRESSION_DEFINED_REGIONS',S4='COMPLETE_QUALIFIED_TLS_ONLY',S5='COMPLETE',new_figures=0,new_mapping=False));(P/'EXIT_STATUS').write_text('0\n');log('COMPLETE_ALL')
try:run()
except Exception:
 (P/'ERROR.txt').write_text(traceback.format_exc());(P/'EXIT_STATUS').write_text('1\n');save('STATUS.json',dict(status='FAILED'));raise
