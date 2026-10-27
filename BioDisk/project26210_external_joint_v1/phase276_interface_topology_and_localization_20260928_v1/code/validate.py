from pathlib import Path
import json,math,hashlib,datetime
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[1]
read=lambda f:pd.read_csv(f,sep='\t',float_precision='round_trip')
s=read(P/'tables/SECTION_EFFECTS.tsv');p=read(P/'tables/PATIENT_EFFECTS.tsv');c=read(P/'tables/COHORT_SUMMARY.tsv')
keys=['dataset','family','endpoint','model','effect_type','quantile']
def eq(a,b,tol=1e-9):
 if pd.isna(a) and pd.isna(b):return
 assert np.isfinite(a) and np.isfinite(b) and abs(a-b)<tol,(a,b)
assert len(s)==38*42 and len(p)==37*42 and len(c)==126
assert not s.duplicated(['sample']+keys).any()
for r in p.itertuples():
 mask=s.patient.eq(r.patient)
 for k in keys:mask &= s[k].eq(getattr(r,k))
 v=s.loc[mask,'effect'].dropna().tolist()
 e=np.nan if not v else math.tanh(sum(math.atanh(max(-.999999999,min(.999999999,x))) for x in v)/len(v)) if r.effect_type=='rho' else sum(v)/len(v)
 eq(e,r.effect);assert len(v)==r.n_sections
def signP(v):
 z=[x for x in v if x!=0];n=len(z)
 if not n:return 1.
 k=min(sum(x>0 for x in z),sum(x<0 for x in z))
 return min(1.,2*sum(math.comb(n,j) for j in range(k+1))/2**n)
for r in c.itertuples():
 mask=np.ones(len(p),bool)
 for k in keys:mask &= p[k].eq(getattr(r,k)).to_numpy()
 v=p.loc[mask,'effect'].dropna().tolist();assert len(v)==r.n
 eq(r.P,signP(v) if len(v)>=10 else np.nan)
 eq(r.median,float(np.median(v)) if len(v) else np.nan)
for _,g in c.groupby(['dataset','family']):
 assert len(g) in [2,4,12,24]
 for pc,qc in [('P','FDR'),('P_wilcoxon','FDR_wilcoxon')]:
  ids=sorted(g.index,key=lambda j:g.at[j,pc] if pd.notna(g.at[j,pc]) else 1.)
  last=1.
  for rank,j in reversed(list(enumerate(ids,1))):
   val=g.at[j,pc] if pd.notna(g.at[j,pc]) else 1.;last=min(last,len(g)*val/rank,1.)
   eq(g.at[j,qc],last if pd.notna(g.at[j,pc]) else np.nan)
core=['depth_mean','depth_delta','macro_logmean','malignant_mean','malignant_delta','log_side_count_ratio'];ext=core+['endothelial_mean','fibroblast_mean']
pairs=read(P/'tables/MATCHED_PAIRS.tsv');top=read(P/'tables/BOUNDARY_TOPOLOGY.tsv');errors=[];lengtherrors=[]
for sample,g in s.groupby('sample'):
 nodes=read(P/'inputs'/f'{sample}_nodes.tsv.gz');z=read(P/'inputs'/f'{sample}_segments.tsv.gz');mem=read(P/'inputs'/f'{sample}_membership.tsv.gz');edges=read(P/'inputs'/f'{sample}_links.tsv.gz')
 xy=nodes[['x_um','y_um']].to_numpy()
 lens=np.linalg.norm(xy[edges.node_a]-xy[edges.node_b],axis=1)
 assert np.allclose(lens,50) and not edges.duplicated().any()
 eq(lens.sum(),nodes.dual_length_um.sum(),1e-7);eq(lens.sum(),mem.length_weight.sum(),1e-7)
 lengths=mem.groupby('segment').length_weight.sum();assert np.allclose(lengths.to_numpy(),z.length_um)
 adj=[[] for _ in range(len(nodes))]
 for row in edges.itertuples():adj[row.node_a].append(row.node_b);adj[row.node_b].append(row.node_a)
 assert max(map(len,adj))<=2
 for q in [.5,.75,.8]:
  for state in ['SPP1','C1QC']:
   high=nodes[f'{state}_high_q{q}'].to_numpy(bool);v=(mem.length_weight.to_numpy()*high[mem.node.to_numpy(int)])
   cov=pd.Series(v).groupby(mem.segment).sum().to_numpy()/z.length_um
   assert np.max(np.abs(cov-z[f'{state}_coverage_q{q}']))<1e-10
   seen=set();longest=0;nr=0
   for j in np.flatnonzero(high&(nodes.dual_length_um>0)):
    if int(j) in seen:continue
    todo=[int(j)];seen.add(int(j));length=0
    while todo:
     a=todo.pop();length+=nodes.dual_length_um.iloc[a]
     for b in adj[a]:
      if high[b] and b not in seen:seen.add(b);todo.append(b)
    nr+=1;longest=max(longest,length)
   tr=top[(top['sample']==sample)&(top.state==state)&(top['quantile']==q)].iloc[0]
   eq(tr.longest_high_run_um,longest,1e-7);assert tr.high_runs==nr
 for row in g[(g['quantile']==.75)&g.effect.notna()&g.effect_type.eq('rho')].itertuples():
  state,target=row.endpoint.split('__')
  cols=[f'{state}_coverage_q0.75',target+'_delta']+([] if row.model=='unadjusted' else core if row.model=='core_adjusted' else ext)
  a=z.loc[z.eligible,cols].dropna().rank(method='average').to_numpy()
  D=np.column_stack([np.ones(len(a)),a[:,2:]]);Q=np.linalg.qr(D,mode='reduced')[0];res=a[:,:2]-Q@(Q.T@a[:,:2])
  rho=float(np.dot(res[:,0],res[:,1])/np.sqrt(np.dot(res[:,0],res[:,0])*np.dot(res[:,1],res[:,1])))
  eq(rho,row.effect,1e-8);errors.append(abs(rho-row.effect))
 for (q,target),pp in pairs[pairs['sample']==sample].groupby(['quantile','target']):
  assert pp.covered_segment.is_unique and pp.gap_segment.is_unique
  assert not set(pp.covered_segment)&set(pp.gap_segment)
  base=z.eligible.to_numpy(bool)&np.isfinite(z[target+'_delta'])&np.isfinite(z[core]).all(1).to_numpy()
  cv=z[core].to_numpy();cv=(cv-cv[base].mean(0))/cv[base].std(0)
  vals=[]
  for row in pp.itertuples():
   h=int(row.covered_segment);l=int(row.gap_segment)
   assert z.contour.iloc[h]==z.contour.iloc[l] and z[f'SPP1_coverage_q{q}'].iloc[h]>=.75 and z[f'SPP1_coverage_q{q}'].iloc[l]<=.25
   dd=float(np.linalg.norm(cv[h]-cv[l]));eq(dd,row.match_distance);assert dd<=2
   val=float(z[target+'_delta'].iloc[l]-z[target+'_delta'].iloc[h]);eq(val,row.gap_minus_covered);vals.append(val)
  r=g[g['quantile'].eq(q)&g.endpoint.eq('SPP1_gap_minus_covered__'+target)].iloc[0]
  eq(r.effect,np.mean(vals) if len(vals)>=3 else np.nan)
assert (P/'logs/exit_status.txt').read_text()=='0'
freeze=json.loads((P/'FREEZE.json').read_text());assert hashlib.sha256((P/'PROTOCOL_BEFORE_TESTING.md').read_bytes()).hexdigest()==freeze['protocol_sha256']
assert not list(P.rglob('*.png')) and not list(P.rglob('*.pdf'))
result=dict(status='PASS',at=datetime.datetime.now().astimezone().isoformat(),sections=38,patients=37,section_effects=len(s),patient_effects=len(p),cohort_endpoints=len(c),FDR_families=c.groupby(['dataset','family']).ngroups,independent_rank_correlations=len(errors),max_correlation_error=max(errors) if errors else None,matched_pair_rows_checked=len(pairs),checks=['midpoint geometry and length conservation','degree<=2 and no artificial closure','direct weighted coverage','independent DFS high-state run lengths','patient aggregation','math exact binomial P','manual BH','QR partial rank correlation','match caliper/no-reuse/contour/raw differences'],figures=0,independent_pathology_or_multiplex_validation='NOT_AVAILABLE_IN_INSPECTED_ASSETS')
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
