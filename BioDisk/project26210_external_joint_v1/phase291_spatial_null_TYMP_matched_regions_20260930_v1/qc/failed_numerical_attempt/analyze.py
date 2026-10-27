from pathlib import Path
import json,hashlib,subprocess,sys,platform
import numpy as np,pandas as pd,scipy
from scipy import stats,sparse,linalg
from scipy.spatial.distance import cdist
from scipy.sparse.csgraph import connected_components,dijkstra
from scipy.optimize import linear_sum_assignment
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
P=Path(__file__).resolve().parents[1]
A=R/'phase289_state_source_boundary_process_20260929_v1'
B=R/'phase290_TYMP_SPP1_stromal_chain_link_20260929_v1'
O=R/'phase234_part4_replan_execution_20260921_v1'
CORE=['in_log1p_UMI','out_log1p_UMI','in_total_macro_fraction','out_total_macro_fraction','in_malignant_fraction','out_malignant_fraction','log_side_count_ratio']
VASC=[a+'_'+b for a in ['in','out'] for b in ['endothelial_fraction','pericyte_fraction']]
TARGETS=[('organization','state_to_CAF','M_state','A_activation'),('organization','tumorTYMP_to_joint','in_g_TYMP','J_state_activation'),('organization','tumorTYMP_to_CAF','in_g_TYMP','A_activation'),('immune','joint_to_T_distribution','J_state_activation','T_delta'),('immune','joint_to_CD8_distribution','J_state_activation','CD8_delta')]
MATCH=[('organization','joint','J_state_activation'),('organization','CAF','A_activation'),('organization','state','M_state'),('immune','T_distribution','T_delta'),('immune','CD8_distribution','CD8_delta')]
def read(f):return pd.read_csv(f,sep='\t',float_precision='round_trip')
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def dump(f,x):Path(f).write_text(json.dumps(x,indent=2,default=lambda o:o.item() if hasattr(o,'item') else str(o)))
def tab(x,n):pd.DataFrame(x).to_csv(P/'tables'/n,sep='\t',index=False)
def thin(e,gap):
 xy=e[['x','y']].to_numpy();keep=[]
 for i in np.lexsort((e.edge.to_numpy(),xy[:,1],xy[:,0])):
  if not keep or np.linalg.norm(xy[keep]-xy[i],axis=1).min()>=gap-1e-8:keep.append(int(i))
 return np.array(keep,int)
def residual(v,D):return v-D@np.linalg.lstsq(D,v,rcond=None)[0]
def bh(v):
 v=np.asarray(v);order=np.argsort(v);out=np.empty(len(v));out[order]=np.minimum(1,np.minimum.accumulate((v[order]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1]);return out
def summarize(sec,name):
 keys=['method','dataset','radius','model','family','endpoint'];records=[]
 for key,g in sec.groupby(keys+['patient'],dropna=False):
  v=g.effect.dropna().to_numpy()
  effect=(float(np.tanh(np.arctanh(np.clip(v,-.999999999,.999999999)).mean())) if key[0]=='spectral' else float(v.mean())) if len(v) else np.nan
  records.append(dict(zip(keys+['patient'],key),effect=effect,sections=';'.join(g.loc[g.effect.notna(),'sample']),n_sections=len(v)))
 pa=pd.DataFrame(records);tab(pa,name+'_PATIENT.tsv');out=[]
 for key,g in pa.groupby(keys,dropna=False):
  v=g.effect.dropna().to_numpy();n=len(v);p=float(stats.wilcoxon(v,zero_method='wilcox',alternative='two-sided',method='auto').pvalue) if n>=10 and np.any(v!=0) else (1. if n>=10 else np.nan)
  out.append(dict(zip(keys,key),patients=n,positive=int((v>0).sum()),negative=int((v<0).sum()),median_effect=float(np.median(v)) if n else np.nan,p_patient=p,registered_family=3 if key[-2]=='organization' else 2))
 su=pd.DataFrame(out);su['FDR_patient']=np.nan;su['estimable_family']=0
 for _,ix in su.groupby(keys[:-1]).groups.items():
  ix=list(ix);ok=su.loc[ix,'p_patient'].notna();jj=list(su.loc[ix].index[ok]);su.loc[ix,'estimable_family']=len(jj)
  if jj:su.loc[jj,'FDR_patient']=bh(su.loc[jj,'p_patient'])
 tab(su,name+'_COHORT.tsv')
def spectral(e,G,radius,meta,model,effects,qcs):
 covs=CORE+(VASC if model=='vascular' else [])
 for fam,ep,xcol,ycol in TARGETS:
  cv=covs+(['T_outer_log'] if ycol=='T_delta' else ['out_CD8_expression','T_outer_log'] if ycol=='CD8_delta' else [])
  valid=np.isfinite(e[[xcol,ycol]+cv]).all(axis=1)
  z=e.loc[valid].copy().reset_index(drop=True)
  base=dict(**meta,method='spectral',model=model,family=fam,endpoint=ep,effect=np.nan,spatial_p=np.nan,regions=len(z),reason='fewer_than_20_regions')
  if len(z)<20:effects.append(base);continue
  ii=z.inner_index.to_numpy(int)
  geo=dijkstra(G,directed=False,indices=ii,limit=4*radius)[:,ii]
  W=np.exp(-geo**2/(2*(2*radius)**2));np.fill_diagonal(W,0)
  keep=W.sum(1)>0;W=W[keep][:,keep];z=z.loc[keep].reset_index(drop=True);n=len(z);base['regions']=n
  if n<20:effects.append(base);continue
  nc,cc=connected_components(sparse.csr_matrix(W),directed=False)
  X=np.column_stack([stats.rankdata(z[c]) for c in cv]+[z.x,z.y])
  X=(X-X.mean(0))/np.where(X.std(0)>1e-12,X.std(0),1)
  D=np.column_stack([np.ones(n),X]+[(cc==k).astype(float) for k in range(1,nc)])
  rank=np.linalg.matrix_rank(D)
  if n-rank<8:base['reason']='residual_dimension_lt8';effects.append(base);continue
  xr=residual(stats.rankdata(z[xcol]),D);yr=residual(stats.rankdata(z[ycol]),D)
  if min(xr.std(),yr.std())<1e-8:base['reason']='constant_residual';effects.append(base);continue
  # Full centered subspace including all zero eigenvalues; avoids lost spectral components.
  Q=linalg.null_space(np.ones((1,n)))
  ev,U=linalg.eigh(Q.T@W@Q);V=Q@U;a=V.T@xr
  seed=int(hashlib.sha256((str(meta)+model+ep).encode()).hexdigest()[:8],16)
  rng=np.random.default_rng(seed);signs=rng.choice([-1.,1.],size=(n-1,999))
  sur=V@(a[:,None]*signs)
  sr=residual(sur,D);den=np.linalg.norm(sr,axis=0)*np.linalg.norm(yr)
  assert (den>1e-10).all()
  null=yr@sr/den;obs=float(xr@yr/(np.linalg.norm(xr)*np.linalg.norm(yr)))
  sp=float((1+(np.abs(null)>=abs(obs)).sum())/1000)
  norm=np.sum(xr*xr);moran=float(n/W.sum()*(xr@W@xr)/norm)
  mi=n/W.sum()*np.sum(sur*(W@sur),axis=0)/np.sum(sur*sur,axis=0)
  mir=n/W.sum()*np.sum(sr*(W@sr),axis=0)/np.sum(sr*sr,axis=0)
  qc=dict(**meta,model=model,endpoint=ep,regions=n,residual_dimension=n-rank,components=nc,moran=moran,maximum_moran_error=float(abs(mi-moran).max()),maximum_norm_error=float(abs(np.sum(sur*sur,axis=0)/norm-1).max()),maximum_spectrum_error=float(abs((V.T@sur)**2-a[:,None]**2).max()),reprojected_moran_shift_median=float(np.median(mir-moran)),reprojected_moran_shift_max=float(abs(mir-moran).max()))
  dump(P/'qc/LAST_NUMERICAL_CHECK.json',dict(**qc,W_asym=float(abs(W-W.T).max()),xr_mean=float(xr.mean()),xr_reconstruction=float(abs(V@a-xr).max())))
  assert qc['maximum_moran_error']<1e-8 and qc['maximum_norm_error']<1e-8
  qcs.append(qc)
  base.update(effect=obs,spatial_p=sp,reason='OK',null_q025=float(np.quantile(null,.025)),null_q975=float(np.quantile(null,.975)))
  effects.append(base)
  np.savez_compressed(P/'inputs'/f"{meta['sample']}_{radius}_{model}_{ep}_SPECTRAL.npz",edge=z.edge.to_numpy(),D=D,xr=xr,yr=yr,W=W,V=V,a=a,signs=signs.astype(np.int8),null=null)
def matching(e,meta,model,pairs,balance,effects,qual):
 covs=CORE+(VASC if model=='vascular' else [])
 base=e[np.isfinite(e[['in_g_TYMP']+covs]).all(axis=1)].copy()
 if len(base):q1,q3=base.in_g_TYMP.quantile([.25,.75]).to_numpy()
 else:q1=q3=np.nan
 z=base.iloc[thin(base,2*meta['radius'])].copy().reset_index(drop=True) if len(base) else base
 selected=[]
 if len(z) and q3>q1:
  v=base[covs].to_numpy();sd=v.std(0);active=sd>1e-12;cv=(z[covs].to_numpy()[:,active]-v[:,active].mean(0))/sd[active]
  hi=np.flatnonzero(z.in_g_TYMP.ge(q3));lo=np.flatnonzero(z.in_g_TYMP.le(q1))
  if len(hi) and len(lo):
   delta=cv[hi,None,:]-cv[None,lo,:];cost=np.sqrt(np.mean(delta**2,axis=2))
   allow=(np.abs(delta).max(axis=2)<=1)&(cost<=.75)&(z.component.to_numpy()[hi,None]==z.component.to_numpy()[None,lo])
   ia,ib=linear_sum_assignment(np.where(allow,cost,1e6))
   selected=[(int(hi[a]),int(lo[b]),float(cost[a,b])) for a,b in zip(ia,ib) if allow[a,b]]
 else:active=np.zeros(len(covs),bool);sd=np.zeros(len(covs))
 used=[];maxsmd=np.nan
 if selected:
  hi=np.array([x[0] for x in selected]);lo=np.array([x[1] for x in selected])
  smd=(z[covs].iloc[hi].mean().to_numpy()-z[covs].iloc[lo].mean().to_numpy())[active]/sd[active]
  maxsmd=float(abs(smd).max()) if len(smd) else 0.
  for col,v in zip(np.array(covs)[active],smd):balance.append(dict(**meta,model=model,covariate=col,SMD=float(v)))
  for k,(i,j,cost) in enumerate(selected):
   row=dict(**meta,model=model,pair=k,high_edge=int(z.edge.iloc[i]),low_edge=int(z.edge.iloc[j]),cost=cost,high_TYMP=float(z.in_g_TYMP.iloc[i]),low_TYMP=float(z.in_g_TYMP.iloc[j]))
   for fam,ep,col in MATCH:
    row[ep+'_high']=float(z[col].iloc[i]);row[ep+'_low']=float(z[col].iloc[j]);row[ep+'_difference']=row[ep+'_high']-row[ep+'_low']
   pairs.append(row);used.append(row)
 ok=len(selected)>=3 and maxsmd<=.25
 qual.append(dict(**meta,model=model,candidates=len(base),nonoverlap_regions=len(z),TYMP_q25=q1,TYMP_q75=q3,pairs=len(selected),max_abs_SMD=maxsmd,eligible=ok,constant_covariates=';'.join(np.array(covs)[~active])))
 for fam,ep,col in MATCH:
  vals=np.array([r[ep+'_difference'] for r in used]);vals=vals[np.isfinite(vals)]
  effect=float(vals.mean()) if ok and len(vals)>=3 else np.nan
  sdout=base[col].std(ddof=0) if len(base) else np.nan
  effects.append(dict(**meta,method='matched',model=model,family=fam,endpoint=ep,effect=effect,standardized_effect=effect/sdout if sdout>0 else np.nan,pairs=len(vals),reason='OK' if np.isfinite(effect) else 'balance_or_pair_count_failed',raw_difference_before_qualification=float(vals.mean()) if len(vals) else np.nan))
 return z
def main():
 subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True)
 dump(P/'qc/SOFTWARE.json',dict(python=sys.version,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,platform=platform.platform()))
 q=read(A/'tables/BOUNDARY_QUALIFICATION.tsv');q=q[(q.definition=='morphology_supported')&q.eligible]
 manifest=[];spec=[];sqc=[];pairs=[];bal=[];mat=[];qual=[];geometry=[]
 for r in q.itertuples():
  paths=[O/'inputs'/f'{r.sample}_INPUTS.tsv.gz',A/'inputs'/f'{r.sample}_SPOT_PROGRAMS.tsv.gz',O/'graphs'/f'{r.sample}_ring.npz',A/'inputs'/f'{r.sample}_morphology_supported_MEMBERSHIPS.json']
  manifest += [dict(path=str(f),sha256=sha(f)) for f in paths]
  d=read(paths[0]);dp=read(paths[1]);assert d.spot_id.tolist()==dp.spot_id.tolist()
  G=sparse.load_npz(paths[2]).tocsr();valid=dp.boundary_input_valid.to_numpy(bool)
  co=G.tocoo();keep=valid[co.row]&valid[co.col];aa=co.row[keep];bb=co.col[keep]
  xy=d[['physical_x_um','physical_y_um']].to_numpy();dist=np.linalg.norm(xy[aa]-xy[bb],axis=1)
  G=sparse.csr_matrix((dist,(aa,bb)),shape=G.shape);G=G.maximum(G.T)
  _,cc=connected_components(G,directed=False)
  mem=json.loads(paths[3].read_text());lookup={(m['radius'],m['edge']):m for m in mem}
  for radius in [200,400]:
   f=B/'inputs'/f'{r.sample}_{radius}_SPP1_mapped_share_CHAIN.tsv.gz';manifest.append(dict(path=str(f),sha256=sha(f)))
   e=read(f);e=e[e.geometry_eligible].copy().reset_index(drop=True);e['component']=cc[e.inner_index.to_numpy(int)]
   for side in ['inside','outside']:
    vals=[]
    for edge in e.edge:
     ids=lookup[(radius,int(edge))][side];v=d.Tcell_q05.iloc[ids].to_numpy();v=v[np.isfinite(v)]
     vals.append(float(v.mean()) if len(v)>=3 else np.nan)
    e['T_'+side]=vals
   assert not (e[['T_inside','T_outside']]<0).any().any()
   e['T_outer_log']=np.log1p(e.T_outside)
   e['T_delta']=np.log1p(e.T_inside)-e.T_outer_log;e['CD8_delta']=e.in_CD8_expression-e.out_CD8_expression
   e.to_csv(P/'inputs'/f'{r.sample}_{radius}_REGIONS.tsv.gz',sep='\t',index=False)
   es=e.iloc[thin(e,radius)].copy().reset_index(drop=True)
   meta=dict(sample=r.sample,dataset=r.dataset,patient=r.patient,radius=radius)
   for model in ['core','vascular']:
    spectral(es,G,radius,meta,model,spec,sqc)
    zm=matching(e,meta,model,pairs,bal,mat,qual)
    used=set()
    for edge in zm.edge:
     m=lookup[(radius,int(edge))];ids=set(m['inside']+m['outside']);assert not used&ids;used|=ids
    geometry.append(dict(**meta,model=model,all_regions=len(e),spectral_regions=len(es),matched_nonoverlap_regions=len(zm),shared_spots=0))
  print('SECTION_COMPLETE',r.sample,flush=True)
 se=pd.DataFrame(spec);se['spatial_FDR']=np.nan
 for _,ix in se.groupby(['sample','radius','model','family']).groups.items():
  jj=[i for i in ix if pd.notna(se.at[i,'spatial_p'])]
  if jj:se.loc[jj,'spatial_FDR']=bh(se.loc[jj,'spatial_p'])
 tab(se,'SPECTRAL_SECTION.tsv');tab(sqc,'SPECTRAL_QC.tsv');tab(pairs,'MATCHED_PAIRS.tsv');tab(bal,'MATCH_BALANCE.tsv');tab(mat,'MATCH_SECTION.tsv');tab(qual,'MATCH_QUALIFICATION.tsv');tab(geometry,'GEOMETRY_QC.tsv')
 summarize(se,'SPECTRAL');summarize(pd.DataFrame(mat),'MATCH')
 manifest={x['path']:x for x in manifest};dump(P/'qc/INPUT_MANIFEST.json',list(manifest.values()))
 for x in manifest.values():assert sha(x['path'])==x['sha256']
 dump(P/'qc/COMPUTED.json',dict(status='COMPUTED_AWAITING_VALIDATION',sections=len(q),spectral_models=len(spec),valid_spectral=sum(x['reason']=='OK' for x in spec),match_models=len(mat),valid_match=sum(x['reason']=='OK' for x in mat),figures=0))
 print('COMPUTED_ALL',flush=True)
if __name__=='__main__':main()
