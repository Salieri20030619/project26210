from pathlib import Path
import json,datetime,hashlib,ast,gc,traceback,subprocess,os,sys
import numpy as np,pandas as pd
from scipy import sparse,stats
from scipy.optimize import linear_sum_assignment
from scipy.sparse.csgraph import connected_components
P=Path(__file__).resolve().parents[1];R=P.parent
OLD=R/'phase234_part4_replan_execution_20260921_v1';MARK=R/'phase246_spatial_marker_structure_20260924_v1'
def now():return datetime.datetime.now().astimezone().isoformat()
def dump(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
def tab(d,f):pd.DataFrame(d).to_csv(P/'tables'/f,sep='\t',index=False)
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def status(s,**kw):dump(P/'STATUS.json',dict(status=s,at=now(),**kw))
def gate(stage):
 q=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
 with open(P/'logs/memory_gate.log','a') as f:f.write(stage+'\n'+q.stdout+q.stderr)
 if q.returncode:raise RuntimeError('MEMORY_GATE')
tree=ast.parse((OLD/'scripts/common.py').read_text())
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['fit','infer']],type_ignores=[]),'inherited','exec'),globals())
KEYS=['dataset','family','endpoint','model','effect_type','quantile']
CORE=['depth_mean','depth_delta','macro_logmean','malignant_mean','malignant_delta','log_side_count_ratio']
EXT=CORE+['endothelial_mean','fibroblast_mean']
def summarize(rows):
 sec=pd.DataFrame(rows);tab(sec,'SECTION_EFFECTS.tsv');out=[]
 for k,g in sec.groupby(KEYS+['patient'],dropna=False,sort=True):
  v=g.effect.dropna().to_numpy();e=(float(np.tanh(np.arctanh(np.clip(v,-.999999999,.999999999)).mean())) if k[4]=='rho' else float(v.mean())) if len(v) else np.nan
  out.append(dict(zip(KEYS+['patient'],k),effect=e,n_sections=len(v),sections=';'.join(g.loc[g.effect.notna(),'sample'])))
 pat=pd.DataFrame(out);tab(pat,'PATIENT_EFFECTS.tsv');out=[]
 for i,(k,g) in enumerate(pat.groupby(KEYS,dropna=False,sort=True)):out.append(dict(zip(KEYS,k),**infer(g.effect,26210276+i)))
 su=pd.DataFrame(out)
 for _,g in su.groupby(['dataset','family']):
  for pc,qc in [('P','FDR'),('P_wilcoxon','FDR_wilcoxon')]:
   su.loc[g.index,qc]=np.where(g[pc].notna(),stats.false_discovery_control(g[pc].fillna(1).to_numpy()),np.nan)
  su.loc[g.index,'family_size_registered']=len(g);su.loc[g.index,'family_size_estimable']=g.P.notna().sum()
 su['evidence']=np.where(su.P.isna(),'DESCRIPTIVE',np.where(su.P<.05,np.where(su.FDR<.05,'FDR_SIGNIFICANT','NOMINAL_SIGNIFICANT'),'NOT_SUPPORTED'))
 tab(su,'COHORT_SUMMARY.tsv')
 loo=[];common=[]
 for k,g in pat[pat['quantile']==.75].groupby(KEYS):
  g=g[g.effect.notna()]
  if len(g)<10:continue
  for rr in g.itertuples():
   v=g.loc[g.patient!=rr.patient,'effect'].to_numpy();nz=np.count_nonzero(v)
   loo.append(dict(zip(KEYS,k),omitted_patient=rr.patient,n=len(v),median=float(np.median(v)),positive=int((v>0).sum()),negative=int((v<0).sum()),P=stats.binomtest(int((v>0).sum()),nz,.5).pvalue if len(v)>=10 and nz else np.nan))
 for k,g in pat.groupby(['dataset','endpoint','model','effect_type']):
  z=g.pivot(index='patient',columns='quantile',values='effect').dropna()
  for q in [.5,.75,.8]:
   v=z[q].to_numpy()
   common.append(dict(zip(['dataset','endpoint','model','effect_type'],k),quantile=q,n=len(v),median=float(np.median(v)) if len(v) else np.nan,positive=int((v>0).sum()),negative=int((v<0).sum())))
 tab(loo,'LEAVE_ONE_PATIENT_OUT.tsv');tab(common,'COMMON_PATIENT_THRESHOLD_SENSITIVITY.tsv')
 return sec,pat,su
def geometry(d,G):
 n=len(d);xy=d[['physical_x_um','physical_y_um']].to_numpy();inside=d.tumor_expression_region.to_numpy(bool)
 neighbors=[set(G.indices[G.indptr[i]:G.indptr[i+1]]) for i in range(n)]
 rr,cc=G.nonzero();assert np.allclose(np.linalg.norm(xy[rr]-xy[cc],axis=1),100)
 crossing=sorted((i,j) for i in range(n) for j in neighbors[i] if i<j and inside[i]!=inside[j])
 lookup={e:i for i,e in enumerate(crossing)};segments=[];triangles=0
 for i in range(n):
  for j in neighbors[i]:
   if j<=i:continue
   for k in neighbors[i].intersection(neighbors[j]):
    if k<=j:continue
    triangles+=1
    es=[tuple(sorted(e)) for e in [(i,j),(j,k),(k,i)] if inside[e[0]]!=inside[e[1]]]
    if es:assert len(es)==2;segments.append((lookup[es[0]],lookup[es[1]]))
 end=np.array(crossing,dtype=int).reshape(-1,2);mid=xy[end].mean(1)
 # Exclude isolated crossing edges from contour lengths, but record them.
 B=sparse.coo_matrix((np.ones(2*len(segments)),(np.array(segments+[(b,a) for a,b in segments])[:,0],np.array(segments+[(b,a) for a,b in segments])[:,1])),shape=(len(end),len(end))).tocsr()
 assert np.max(np.diff(B.indptr),initial=0)<=2 and np.max(B.data,initial=0)<=1
 ncomp,lab=connected_components(B,directed=False);dual=np.zeros(len(end));links=[]
 for a,b in segments:
  le=float(np.linalg.norm(mid[a]-mid[b]));assert 49.999<le<50.001
  dual[a]+=le/2;dual[b]+=le/2;links.append(dict(node_a=a,node_b=b,length_um=le))
 memberships=[];contours=[];segmeta=[];sid=0
 for comp in range(ncomp):
  ids=np.flatnonzero(lab==comp);active=[v for v in ids if B.indptr[v+1]>B.indptr[v]]
  if not active:continue
  endpoints=[v for v in active if B.indptr[v+1]-B.indptr[v]==1]
  key=lambda v:(mid[v,0],mid[v,1],v)
  start=min(endpoints if endpoints else active,key=key);path=[start];previous=-1;cur=start
  while True:
   ns=[int(v) for v in B.indices[B.indptr[cur]:B.indptr[cur+1]] if v!=previous]
   if not ns:break
   nxt=min(ns,key=key)
   if nxt==start:break
   assert nxt not in path
   path.append(nxt);previous,cur=cur,nxt
  assert len(path)==len(active)
  closed=len(endpoints)==0;pairs=list(zip(path[:-1],path[1:]))+([(path[-1],path[0])] if closed else [])
  pieces=[];arc=0.
  for a,b in pairs:
   le=float(np.linalg.norm(mid[a]-mid[b]))
   for node,left,right in [(a,arc,arc+le/2),(b,arc+le/2,arc+le)]:
    while right-left>1e-7:
     bn=int(np.floor((left+1e-7)/400));stop=min(right,(bn+1)*400.)
     pieces.append((bn,node,stop-left));left=stop
   arc+=le
  contours.append(dict(contour=comp,nodes=len(path),closed=closed,length_um=arc))
  pp=pd.DataFrame(pieces,columns=['bin','node','weight'])
  for bn,g in pp.groupby('bin',sort=True):
   length=float(g.weight.sum())
   segmeta.append(dict(segment=sid,contour=comp,bin=int(bn),length_um=length,geometry_eligible=length>=200-1e-6))
   for node,v in g.groupby('node').weight.sum().items():memberships.append(dict(segment=sid,node=int(node),length_weight=float(v)))
   sid+=1
 assert abs(sum(x['length_um'] for x in contours)-dual.sum())<1e-7
 return end,mid,B,lab,dual,pd.DataFrame(links),pd.DataFrame(contours),pd.DataFrame(segmeta),pd.DataFrame(memberships),triangles
def neighborhoods(end,D,n):
 out=[];distance=[];assign=np.full(n,-1,int);best=np.full(n,np.inf)
 for node,(a,b) in enumerate(end):
  near={int(a):50.,int(b):50.}
  for v in [a,b]:
   for j,val in zip(D.indices[D.indptr[v]:D.indptr[v+1]],D.data[D.indptr[v]:D.indptr[v+1]]):
    vv=float(val)+50
    if vv<=400+1e-6:near[int(j)]=min(near.get(int(j),np.inf),vv)
  ids=np.array(sorted(near));dis=np.array([near[j] for j in ids]);out.append(ids);distance.append(dis)
  update=dis<best[ids]-1e-6;best[ids[update]]=dis[update];assign[ids[update]]=node
 return out,distance,assign,best
def runs(high,B,dual):
 ids=np.flatnonzero(high&(dual>0))
 if not len(ids):return 0.,0
 nr,lab=connected_components(B[ids][:,ids],directed=False)
 return float(np.bincount(lab,weights=dual[ids]).max()),int(nr)

def main():
 gate('start');freeze=json.loads((P/'FREEZE.json').read_text());assert sha(P/'PROTOCOL_BEFORE_TESTING.md')==freeze['protocol_sha256']
 regs=pd.read_csv(OLD/'tables/INPUT_SECTION_REGISTER.tsv',sep='\t');assert len(regs)==38 and regs.patient.nunique()==37
 regs.to_csv(P/'tables/INPUT_SECTION_REGISTER.tsv',sep='\t',index=False)
 effects=[];topology=[];thresholds=[];geomqc=[];pairrows=[];matchqc=[];hashes=[];extension=[]
 for i,r in enumerate(regs.itertuples()):
  gate(r.sample);status('RUNNING',current=r.sample,completed=i,total=38)
  files=[OLD/'inputs'/f'{r.sample}_INPUTS.tsv.gz',OLD/'graphs'/f'{r.sample}_ring.npz',OLD/'graphs'/f'{r.sample}_distance800.npz',MARK/'inputs'/f'{r.sample}_MARKERS.tsv.gz']
  hashes.extend(dict(path=str(f),sha256=sha(f)) for f in files)
  d=pd.read_csv(files[0],sep='\t');m=pd.read_csv(files[3],sep='\t');assert np.array_equal(d.spot_id,m.spot_id)
  G=sparse.load_npz(files[1]).tocsr();D=sparse.load_npz(files[2]).tocsr();inside=d.tumor_expression_region.to_numpy(bool)
  assert np.array_equal(inside,d.malignant_fraction.to_numpy()>=.2) and (G-G.T).nnz==0
  end,mid,B,lab,dual,links,contours,segs,membership,triangles=geometry(d,G)
  neigh,ndist,assign,bdist=neighborhoods(end,D,len(d))
  rr,cc=B.nonzero();assert np.all(lab[rr]==lab[cc])
  nodes=pd.DataFrame(dict(node=np.arange(len(end)),spot_a=end[:,0],spot_b=end[:,1],x_um=mid[:,0],y_um=mid[:,1],contour=lab,dual_length_um=dual))
  macro=d.total_macro_q05.to_numpy();total=d.total_q05.to_numpy()
  abundance={'SPP1':d.V1_SPP1_q05.to_numpy(),'C1QC':(d.C1QC_fraction*d.total_q05).to_numpy()}
  spotshare={k:v/macro for k,v in abundance.items()}
  assert np.max(np.abs(spotshare['SPP1']-d.V3_SPP1_macro_share.to_numpy()))<1e-10
  macrofrac=macro/total;rich=macrofrac>=.05
  nrich=(macro[end].sum(1)/total[end].sum(1))>=.05
  nodes['macro_rich']=nrich
  for state,a in abundance.items():nodes[state+'_share']=a[end].sum(1)/macro[end].sum(1)
  ys={'T_abundance_log2':np.log2(d.Tcell_q05.to_numpy()) if (d.Tcell_q05>0).all() else np.full(len(d),np.nan),'CD8_RNA':m.CD8lineage.to_numpy()}
  depth=d.log1p_UMI.to_numpy();mal=d.malignant_fraction.to_numpy();endo=d.endothelial_fraction.to_numpy();fib=m.Fibroblast.to_numpy()
  byseg={int(sid):gg.node.to_numpy(int) for sid,gg in membership.groupby('segment')}
  primary_segment=np.full(len(end),-1,int);bestweight=np.full(len(end),-1.)
  for row in membership.sort_values(['segment','node']).itertuples():
   if row.length_weight>bestweight[row.node]+1e-6:primary_segment[row.node]=row.segment;bestweight[row.node]=row.length_weight
  spotsegment=np.full(len(d),-1,int);valid=assign>=0;spotsegment[valid]=primary_segment[assign[valid]]
  segrows=[]
  for sr in segs.itertuples():
   ns=byseg[sr.segment];ids=np.unique(np.concatenate([neigh[j] for j in ns]));ii=ids[inside[ids]];oo=ids[~inside[ids]]
   row=dict(segment=sr.segment,contour=sr.contour,length_um=sr.length_um,geometry_eligible=sr.geometry_eligible,n_inside=len(ii),n_outside=len(oo))
   good=bool(sr.geometry_eligible and min(len(ii),len(oo))>=10);row['eligible']=good
   mean=lambda v,j:float(np.mean(v[j])) if len(j) else np.nan
   row.update(depth_mean=mean(depth,ids),depth_delta=mean(depth,ii)-mean(depth,oo),macro_logmean=float(np.log1p(mean(macro,ids))),malignant_mean=mean(mal,ids),malignant_delta=mean(mal,ii)-mean(mal,oo),log_side_count_ratio=float(np.log(max(len(ii),1)/max(len(oo),1))),endothelial_mean=mean(endo,ids),fibroblast_mean=mean(fib,ids))
   for target,y in ys.items():
    row[target+'_inside']=mean(y,ii);row[target+'_outside']=mean(y,oo);row[target+'_delta']=mean(y,ii)-mean(y,oo)
   segrows.append(row)
  z=pd.DataFrame(segrows)
  M=sparse.coo_matrix((membership.length_weight,(membership.segment,membership.node)),shape=(len(segs),len(end))).tocsr()
  assert np.allclose(np.asarray(M.sum(1)).ravel(),z.length_um)
  for q in [.75,.5,.8]:
   for state in ['SPP1','C1QC']:
    vals=spotshare[state][rich];qualified=len(vals)>=100 and len(np.unique(vals))>=3
    cutoff=float(np.quantile(vals,q)) if qualified else np.nan
    high=nrich&(nodes[state+'_share'].to_numpy()>=cutoff) if qualified else np.zeros(len(end),bool)
    nodes[f'{state}_high_q{q}']=high
    cover=np.asarray(M@high.astype(float)).ravel()/z.length_um.to_numpy() if qualified else np.full(len(z),np.nan)
    z[f'{state}_coverage_q{q}']=cover
    longest,nruns=runs(high,B,dual);longgap,ngaps=runs(~high,B,dual)
    thresholds.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,state=state,quantile=q,threshold=cutoff,macro_rich_spots=int(rich.sum()),eligible=qualified))
    topology.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,state=state,quantile=q,threshold=cutoff,boundary_length_um=float(dual.sum()),coverage=float(np.dot(dual,high)/dual.sum()) if qualified and dual.sum()>0 else np.nan,longest_high_run_um=longest,high_runs=nruns,longest_low_run_um=longgap,low_runs=ngaps,closed_contours=int(contours.closed.sum()),contours=len(contours)))
    rawhigh=rich&(spotshare[state]>=cutoff)
    for sr in segs.itertuples():
     loc=(spotsegment==sr.segment)&rawhigh
     vi=bdist[loc&inside];vo=bdist[loc&~inside]
     extension.append(dict(sample=r.sample,segment=sr.segment,state=state,quantile=q,high_spots_inside=len(vi),high_spots_outside=len(vo),inside_extent_q90_um=float(np.quantile(vi,.9)) if len(vi)>=3 else np.nan,outside_extent_q90_um=float(np.quantile(vo,.9)) if len(vo)>=3 else np.nan,definition='assigned high-state spots within400um graph radius; not physical barrier thickness'))
    for target in ys:
     for model in ['unadjusted','core_adjusted','extended_adjusted']:
      cov=None if model=='unadjusted' else z[CORE if model=='core_adjusted' else EXT].to_numpy()
      e,n,why=fit(cover,z[target+'_delta'].to_numpy(),cov=cov,mask=z.eligible.to_numpy(bool),min_n=30 if model=='extended_adjusted' else 15)
      effects.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,family='coverage_association'+('_primary' if q==.75 else '_sensitivity'),endpoint=state+'__'+target,model=model,effect_type='rho',quantile=q,effect=e,eligible_segments=n,reason=why))
   for target in ys:
    base=z.eligible.to_numpy(bool)&np.isfinite(z[target+'_delta'])&np.isfinite(z[CORE]).all(1).to_numpy()
    cv=z[CORE].to_numpy();sd=cv[base].std(0) if base.any() else np.zeros(len(CORE));pairs=[];why=''
    highidx=np.flatnonzero(base&(z[f'SPP1_coverage_q{q}']>=.75));lowidx=np.flatnonzero(base&(z[f'SPP1_coverage_q{q}']<=.25))
    if not len(highidx) or not len(lowidx):why='NO_COVERED_OR_GAP_SEGMENTS'
    elif (sd==0).any():why='CONSTANT_MATCH_COVARIATE'
    else:
     cv=(cv-cv[base].mean(0))/sd
     cost=np.linalg.norm(cv[highidx,None,:]-cv[None,lowidx,:],axis=2)
     allowed=(z.contour.to_numpy()[highidx,None]==z.contour.to_numpy()[None,lowidx])&(cost<=2)
     aa,bb=linear_sum_assignment(np.where(allowed,cost,1e6))
     for a,b in zip(aa,bb):
      if not allowed[a,b]:continue
      h=int(highidx[a]);l=int(lowidx[b]);val=float(z[target+'_delta'].iloc[l]-z[target+'_delta'].iloc[h]);pairs.append(val)
      pairrows.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,quantile=q,target=target,covered_segment=h,gap_segment=l,contour=int(z.contour.iloc[h]),match_distance=float(cost[a,b]),covered_value=float(z[target+'_delta'].iloc[h]),gap_value=float(z[target+'_delta'].iloc[l]),gap_minus_covered=val))
    if len(pairs)<3:
     effect=np.nan
     if not why:why='FEWER_THAN_THREE_MATCHED_PAIRS'
    else:effect=float(np.mean(pairs))
    effects.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,family='matched_gaps'+('_primary' if q==.75 else '_sensitivity'),endpoint='SPP1_gap_minus_covered__'+target,model='same_contour_core6_caliper2',effect_type='difference',quantile=q,effect=effect,eligible_segments=2*len(pairs),reason=why))
    matchqc.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,quantile=q,target=target,covered_candidates=len(highidx),gap_candidates=len(lowidx),accepted_pairs=len(pairs),eligible=len(pairs)>=3,reason=why))
  for name,dd in [('nodes',nodes),('links',links),('contours',contours),('segments',z),('membership',membership)]:
   dd.to_csv(P/'inputs'/f'{r.sample}_{name}.tsv.gz',sep='\t',index=False)
  geomqc.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,spots=len(d),triangles=triangles,boundary_nodes=len(end),isolated_crossings=int((dual==0).sum()),boundary_links=len(links),max_degree=int(np.diff(B.indptr).max(initial=0)),boundary_length_um=float(dual.sum()),membership_length_um=float(membership.length_weight.sum()),contours=len(contours),closed_contours=int(contours.closed.sum()),all_segments=len(z),length_eligible_segments=int(z.geometry_eligible.sum()),both_sides_eligible_segments=int(z.eligible.sum()),outside_graph_bridges=0))
  print('SECTION_COMPLETE',i+1,r.sample,'segments',int(z.eligible.sum()),flush=True)
  del D,G,B,M;gc.collect()
 tab(topology,'BOUNDARY_TOPOLOGY.tsv');tab(thresholds,'FROZEN_STATE_THRESHOLDS.tsv');tab(geomqc,'GEOMETRY_QC.tsv');tab(pairrows,'MATCHED_PAIRS.tsv');tab(matchqc,'MATCHING_ELIGIBILITY.tsv');tab(extension,'TRANSVERSE_EXTENT_DESCRIPTIVE.tsv')
 sec,pat,su=summarize(effects)
 for row in hashes:assert sha(row['path'])==row['sha256']
 dump(P/'qc/READONLY_INPUT_SHA256.json',hashes)
 dump(P/'qc/COMPUTATION_COMPLETED.json',dict(status='COMPUTED_AWAITING_VALIDATION',at=now(),sections=len(regs),patients=regs.patient.nunique(),section_effects=len(sec),patient_effects=len(pat),cohort_endpoints=len(su),figures=0))
 status('COMPUTED_AWAITING_VALIDATION',completed=38,total=38);print('ALL_COMPUTATION_COMPLETED',len(su),flush=True)
if __name__=='__main__':
 try:main()
 except Exception:
  status('FAILED',traceback=traceback.format_exc());raise
