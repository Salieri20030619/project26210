from common import *
import json
m=pd.read_csv(D/'inputs/MACRO_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str})
dom=json.loads((D/'inputs/FROZEN_UNION.json').read_text())['dominant_C10_patient']
cache=pd.read_csv(R/'phase179_STAT3_external_activity_validation_20260916_v1/run_v2/SC_PRIMARY_SCORE_CACHE.tsv.gz',sep='\t')
defs={'U':['4','10','15'],'C4':['4'],'U_no10':['4','15'],'U_no15':['4','10']}
rows=[]
for resource,ss in cache.groupby('resource'):
 for pat,z in ss.groupby('patient'):
  z=z.set_index('cell_id');md=m.loc[z.index];v=z.score.to_numpy();cs=md.candidate
  for label,members in defs.items():
   hit=cs.isin(members).to_numpy()
   if min(hit.sum(),(~hit).sum())<20:continue
   rows.append(dict(resource=resource,group=label,reference='own_complement',patient=pat,study=md.study.iloc[0],target_score=v[hit].mean(),effect=v[hit].mean()-v[~hit].mean()))
  # Same non-U reference for two definitions.
  back=~cs.isin(defs['U']).to_numpy()
  if back.sum()<20:continue
  for label in ['U','C4']:
   hit=cs.isin(defs[label]).to_numpy()
   if hit.sum()>=20:rows.append(dict(resource=resource,group=label,reference='nonU',patient=pat,study=md.study.iloc[0],target_score=v[hit].mean(),effect=v[hit].mean()-v[back].mean()))
pt=pd.DataFrame(rows);pt.to_csv(D/'receiver/STAT3_ATTRIBUTION_PATIENTS.tsv',sep='\t',index=False)
out=[];delta=[];los=[]
for (resource,ref),z in pt.groupby(['resource','reference']):
 shared=set(z.loc[z.group=='U','patient'])&set(z.loc[z.group=='C4','patient'])
 for (label,),g in z.groupby(['group']):
  for pop,keep in [('own',g),('common_U_C4',g[g.patient.isin(shared)]),('without_C10_dominant',g[g.patient!=dom])]:
   out.append(dict(resource=resource,reference=ref,group=label,population=pop,**summary(keep.effect)))
  for st in g.study.unique():los.append(dict(resource=resource,reference=ref,group=label,omitted_study=st,**summary(g.loc[g.study!=st,'effect'])))
 w=z[z.group.isin(['U','C4'])].pivot(index='patient',columns='group',values='effect').dropna()
 v=(w.U-w.C4).to_numpy();rng=np.random.default_rng(262205);bs=np.median(v[rng.integers(0,len(v),(5000,len(v)))],axis=1)
 delta.append(dict(resource=resource,reference=ref,U_median=w.U.median(),C4_median=w.C4.median(),CI_low=np.quantile(bs,.025),CI_high=np.quantile(bs,.975),**summary(v)))
add_fdr(pd.DataFrame(out),['resource','reference','population']).to_csv(D/'receiver/STAT3_ATTRIBUTION_SUMMARY.tsv',sep='\t',index=False)
add_fdr(pd.DataFrame(delta),['reference']).to_csv(D/'receiver/STAT3_SHARED_PATIENT_CHANGE.tsv',sep='\t',index=False)
pd.DataFrame(los).to_csv(D/'receiver/STAT3_LEAVE_STUDY.tsv',sep='\t',index=False)
# Side-by-side attribution for the functional axes that actually changed.
f=pd.read_csv(D/'functional/ALL_PROGRAM_RESULTS.tsv.gz',sep='\t')
terms=['Glycolysis','E2F Targets','G2-M Checkpoint','Cholesterol metabolism','Lysosome']
sel=f[(f.version=='without_SPP1_TYMP')&f.term.isin(terms)&f.contrast.isin(['all_U','all_C4','all_U_no10','all_U_no15','all_U_without_C10_dominant_patient'])]
sel.to_csv(D/'functional/CHANGED_AXIS_COMPONENT_ATTRIBUTION.tsv',sep='\t',index=False)
print(pd.DataFrame(out).query("resource == 'collectri' and reference == 'own_complement'").to_string(index=False))
print('STAT3_SHARED',pd.DataFrame(delta).to_string(index=False))
(D/'qc/ATTRIBUTION_STATUS.json').write_text(json.dumps({'status':'PASS','definitions_fixed_before_primary_results':True,'same_background_compared':True,'do_not_select_best_component_subset':True},indent=2))
