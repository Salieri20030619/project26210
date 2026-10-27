from common import *
from scipy import sparse,stats,linalg
PRE=R/'phase289_state_source_boundary_process_20260929_v1'
OLD=R/'phase234_part4_replan_execution_20260921_v1'
HARM=R/'phase277_metabolic_macrophage_barrier_cooperation_20260928_v1/harmonized_platform_v1'
CORE=['in_log1p_UMI','out_log1p_UMI','in_total_macro_fraction','out_total_macro_fraction','in_malignant_fraction','out_malignant_fraction','log_side_count_ratio']
VASC=['in_endothelial_fraction','out_endothelial_fraction','in_pericyte_fraction','out_pericyte_fraction']
def main():
 gate('chain_validation')
 q=read(P/'tables/INHERITED_BOUNDARY_QUALIFICATION.tsv');prior=read(PRE/'tables/BOUNDARY_QUALIFICATION.tsv');prior=prior[prior.definition.eq('morphology_supported')].reset_index(drop=True)
 pd.testing.assert_frame_equal(q,prior)
 ef=read(P/'tables/CHAIN_SECTION_EFFECTS.tsv');maxmean=0.;maxrho=0.;count=0;inputs=0
 for r in q.itertuples():
  mf=PRE/'inputs'/f'{r.sample}_morphology_supported_MEMBERSHIPS.json'
  if not mf.exists():continue
  d=read(OLD/'inputs'/f'{r.sample}_INPUTS.tsv.gz');h=read(HARM/'inputs'/f'{r.sample}_PROGRAMS.tsv.gz').set_index('spot_id').reindex(d.spot_id)
  finite=d.V3_SPP1_macro_share.dropna()
  assert finite.ge(-1e-10).all() and finite.le(1+1e-10).all()
  S=np.column_stack([d.V3_SPP1_macro_share,h.C4_platform15_noAPOEshared])
  mem=json.loads(mf.read_text())
  for radius in [200,400]:
   f=PRE/'inputs'/f'{r.sample}_morphology_supported_{radius}_EDGES.tsv.gz'
   if not f.exists():continue
   old=read(f);look={v['edge']:v['macro_out'] for v in mem if v['radius']==radius}
   rr=[];cc=[]
   for k,key in enumerate(old.edge):
    ids=look[int(key)];rr.extend([k]*len(ids));cc.extend(ids)
   G=sparse.csr_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(old),len(d)))
   denom=G@np.isfinite(S).astype(int);sums=G@np.nan_to_num(S)
   means=np.divide(sums,denom,out=np.full_like(sums,np.nan),where=denom>=3)
   for j,readout in enumerate(['SPP1_mapped_share','C4_disjoint_program']):
    e=read(P/'inputs'/f'{r.sample}_{radius}_{readout}_CHAIN.tsv.gz')
    pd.testing.assert_frame_equal(e[old.columns],old)
    assert np.array_equal(e.state_valid_spots,denom[:,j])
    assert np.allclose(means[:,j],e.M_state,equal_nan=True,atol=1e-12,rtol=1e-12)
    err=np.abs(means[:,j]-e.M_state.to_numpy())
    if np.isfinite(err).any():maxmean=max(maxmean,float(np.nanmax(err)))
    mask=e.geometry_eligible.to_numpy()&bool(r.eligible);m=e.M_state.to_numpy();a=e.A_activation.to_numpy()
    assert np.array_equal(a,e.out_p_CAF_activation.to_numpy(),equal_nan=True)
    common=mask&np.isfinite(m)&np.isfinite(a)
    ok=common.sum()>=20 and len(np.unique(m[common]))>=3 and len(np.unique(a[common]))>=3 and (m[common]>0).any() and (a[common]>0).any()
    joint=np.full(len(e),np.nan)
    if ok:joint[common]=(stats.rankdata(m[common])/(common.sum()+1))*(stats.rankdata(a[common])/(common.sum()+1))
    assert np.array_equal(joint,e.J_state_activation.to_numpy(),equal_nan=True)
    assert np.array_equal(e.joint_eligible.to_numpy(),common&ok)
    tests=ef[(ef['sample']==r.sample)&(ef.radius==radius)&(ef.readout==readout)]
    for t in tests.itertuples():
     covs=CORE+(VASC if t.model=='vascular_sensitivity' else [])
     interaction=t.family=='interaction'
     if t.endpoint.startswith('tumor'):x=e.in_g_TYMP.to_numpy()
     elif t.endpoint.startswith('macro'):x=e.macro_out_g_TYMP.to_numpy()
     else:x=joint
     if interaction:y=a
     elif t.endpoint=='joint_to_inner_CD8':y=e.in_CD8_expression.to_numpy();covs=covs+['out_CD8_expression','out_Tcell_fraction']
     elif t.endpoint=='joint_to_inner_T':y=e.in_Tcell_fraction.to_numpy();covs=covs+['out_Tcell_fraction']
     else:y=joint
     C=e[covs].to_numpy();use=mask&np.isfinite(x)&np.isfinite(y)&np.isfinite(C).all(1)
     if interaction:use &= np.isfinite(m)
     assert use.sum()==t.edges
     if not np.isfinite(t.effect):continue
     xx=stats.rankdata(x[use]);yy=stats.rankdata(y[use]);cov=stats.rankdata(C[use],axis=0);cov=cov[:,cov.std(0)>1e-9];cov=(cov-cov.mean(0))/cov.std(0) if cov.shape[1] else cov
     if interaction:
      mm=stats.rankdata(m[use]);xx=(xx-xx.mean())/xx.std();mm=(mm-mm.mean())/mm.std();cov=np.column_stack([cov,xx,mm]);xx=xx*mm
     D=np.column_stack([np.ones(use.sum()),cov]);pair=np.column_stack([xx,yy]);coef=linalg.lstsq(D,pair,cond=max(D.shape)*np.finfo(float).eps,lapack_driver='gelsy')[0]
     residual=pair-D@coef;rho=float(np.corrcoef(residual.T)[0,1]);maxrho=max(maxrho,abs(rho-t.effect));assert abs(rho-t.effect)<1e-9,(r.sample,t.endpoint,rho,t.effect);count+=1
    inputs+=1
  print('CHAIN_VALIDATED',r.sample,flush=True)
 keys=['dataset','radius','readout','model','family','endpoint'];pt=read(P/'tables/CHAIN_PATIENT_EFFECTS.tsv');su=read(P/'tables/CHAIN_COHORT_SUMMARY.tsv')
 for row in pt.to_dict('records'):
  mask=np.ones(len(ef),bool)
  for k in keys+['patient']:mask &= ef[k].eq(row[k])
  vals=ef.loc[mask,'effect'].dropna().to_numpy();value=np.tanh(np.arctanh(np.clip(vals,-.999999,.999999)).mean()) if len(vals) else np.nan
  assert np.isclose(value,row['effect'],equal_nan=True)
 for row in su.to_dict('records'):
  mask=np.ones(len(pt),bool)
  for k in keys:mask &= pt[k].eq(row[k])
  vals=pt.loc[mask,'effect'].dropna().to_numpy();n=len(vals);assert n==row['patients']
  p=stats.wilcoxon(vals).pvalue if n>=10 and (vals!=0).any() else (1. if n>=10 else np.nan)
  assert np.isclose(p,row['P'],equal_nan=True)
  assert int((vals>0).sum())==row['positive'] and int((vals<0).sum())==row['negative']
 for _,g in su.groupby(['dataset','radius','readout','model','family']):
  z=g[g.P.notna()]
  if len(z):assert np.allclose(stats.false_discovery_control(z.P),z.FDR,atol=1e-12)
 detail=read(P/'tables/SAME_PATIENT_CHAIN_DETAILS.tsv');summ=read(P/'tables/SAME_PATIENT_CHAIN_CONCORDANCE.tsv')
 for row in summ.to_dict('records'):
  mask=np.ones(len(detail),bool)
  for k in ['dataset','radius','readout','model','source','immune']:mask &= detail[k].eq(row[k])
  z=detail[mask];assert len(z)==row['paired_patients']
  assert ((z.TYMP_to_joint>0)&(z.joint_to_immune<0)).sum()==row['compatible']
  assert row['compatible']+row['positive_positive']+row['negative_negative']+row['negative_positive']+row['zero_any']==row['paired_patients']
  for r in z.itertuples():
   cand=pt[(pt.dataset==r.dataset)&(pt.radius==r.radius)&(pt.readout==r.readout)&(pt.model==r.model)&(pt.patient==r.patient)]
   up=cand.loc[cand.endpoint==r.source+'_TYMP_to_joint','effect'].iloc[0]
   down=cand.loc[cand.endpoint=='joint_to_inner_'+r.immune,'effect'].iloc[0]
   assert up==r.TYMP_to_joint and down==r.joint_to_immune
 for record in json.loads((P/'qc/INPUT_MANIFEST.json').read_text()):assert sha(record['path'])==record['sha256']
 assert sha(P/'PROTOCOL_BEFORE_TESTING.md')==json.loads((P/'FREEZE.json').read_text())['protocol_sha256']
 assert (P/'logs/chain_exit.txt').read_text().strip()=='0'
 dump(P/'qc/INDEPENDENT_VALIDATION.json',dict(status='PASS',edge_input_variants=inputs,source_state_sparse_mean_max_error=maxmean,all_J_and_eligibility_match=True,local_effects_checked=count,max_partial_rho_difference=maxrho,all_patient_aggregates_P_FDR_match=True,same_patient_chain_pairs_exact=True,all_sources_and_protocol_unchanged=True))
 print('CHAIN_VALIDATION_PASS',flush=True)
if __name__=='__main__':main()
