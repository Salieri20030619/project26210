from common import *
from adjacency import CORE,OLD,H
from validate import ref
def main():
 assert (P/'logs/adjacency_exit.txt').read_text().strip()=='0'
 effects=read(P/'tables/ADJACENCY_SECTION_EFFECTS.tsv');mx=0.;nc=0;aggmax=0.;nagg=0
 for (sample,radius),g in effects.groupby(['sample','radius']):
  z=read(P/'inputs'/f'{sample}_ADJACENCY_{radius}.tsv.gz');d=read(OLD/'inputs'/f'{sample}_INPUTS.tsv.gz');D=sparse.load_npz(OLD/'graphs'/f'{sample}_distance800.npz').tocsr()
  for r in z.itertuples():
   ids=np.array([int(x) for x in str(r.neighbor_indices).split(';')]);i=int(r.spot_index);assert i not in ids and len(ids)>=3 and len(set(ids))==len(ids)
   dist=np.asarray(D[i,ids].toarray()).ravel();assert np.all(dist>0)&np.all(dist<=radius+1e-6)&np.all(d.malignant_fraction.iloc[ids]<.2)
   for name in ['POSTNCAF','ECM4']:
    value=float(sum(d[name].iloc[ids])/len(ids));er=abs(value-getattr(r,'receiver_'+name));assert er<1e-10;aggmax=max(aggmax,er);nagg+=1
  for r in g[g.effect.notna()].itertuples():
   cov=CORE if r.model=='composition_adjusted' else CORE+['source_'+r.target,'source_endo','receiver_endo'];v=z[['source_TYMP','receiver_'+r.target]+cov].dropna()
   value=ref(v.source_TYMP,v['receiver_'+r.target],v[cov].to_numpy());assert abs(value-r.effect)<1e-10;mx=max(mx,abs(value-r.effect));nc+=1
 su=read(P/'tables/ADJACENCY_COHORT_SUMMARY.tsv');pa=read(P/'tables/ADJACENCY_PATIENT_EFFECTS.tsv');keys=['dataset','family','radius','target','model']
 for row in pa.to_dict('records'):
  match=np.ones(len(effects),bool)
  for k in keys+['patient']:match &= effects[k].eq(row[k]).to_numpy()
  vv=effects.loc[match,'effect'].dropna().to_numpy()
  if len(vv):assert abs(float(np.tanh(np.mean(np.arctanh(np.clip(vv,-.99999999,.99999999)))))-row['effect'])<1e-12
  else:assert pd.isna(row['effect'])
 for row in su.to_dict('records'):
  match=np.ones(len(pa),bool)
  for k in keys:match &= pa[k].eq(row[k]).to_numpy()
  v=pa.loc[match,'effect'].dropna().to_numpy();assert len(v)==row['n'] and int((v>0).sum())==row['positive']
  if len(v)>=10:assert abs(stats.binomtest(int((v>0).sum()),int((v!=0).sum()),.5).pvalue-row['P'])<1e-12
  else:assert pd.isna(row['P'])
 for _,g in su.groupby(['dataset','family']):
  q=stats.false_discovery_control(g.P.fillna(1).to_numpy());ok=g.P.notna().to_numpy();assert np.allclose(q[ok],g.FDR.to_numpy()[ok],atol=1e-12)
 for rec in json.loads((P/'qc/ADJACENCY_INPUTS.json').read_text()):assert sha(rec['path'])==rec['sha256']
 dump(P/'qc/ADJACENCY_VALIDATION.json',dict(status='PASS',independent_rhos=nc,max_rho_error=mx,neighbor_aggregates=nagg,max_aggregate_error=aggmax,self_pairs=0,summary_endpoints=len(su),all_input_hashes_unchanged=True))
 print('ADJACENCY_VALIDATION_PASS',nc,nagg,mx,flush=True)
if __name__=='__main__':main()
