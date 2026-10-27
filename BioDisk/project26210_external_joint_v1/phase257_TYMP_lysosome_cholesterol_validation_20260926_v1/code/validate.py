from common import *
from scipy import stats
def bh(p):
 p=np.asarray(p,float);ix=np.argsort(p,kind='stable');s=p[ix]*len(p)/np.arange(1,len(p)+1)
 v=np.minimum(1,np.minimum.accumulate(s[::-1])[::-1]);q=np.empty(len(p));q[ix]=v;return q
def main():
 gate('final_validation')
 assert (P/'analysis_EXIT_STATUS').read_text().strip()=='0'
 for name in ['SPATIAL_COMPLETE.json','PAIRED_COMPLETE.json']:assert json.loads((P/name).read_text())['status']=='PASS'
 manifests=[]
 for name in ['SPATIAL_SOURCE_MANIFEST.json','PAIRED_SOURCE_MANIFEST.json']:manifests+=json.loads((P/'qc'/name).read_text())
 unique={m['path']:m for m in manifests}
 for path,old in unique.items():assert fingerprint(path)==old,path
 jwrite('qc/INPUTS_AFTER.json',list(unique.values()))
 sec=pd.read_csv(P/'tables/SPATIAL_SECTION_EFFECTS.tsv',sep='\t')
 pat=pd.read_csv(P/'tables/SPATIAL_PATIENT_EFFECTS.tsv',sep='\t')
 sm=pd.read_csv(P/'tables/SPATIAL_COHORT_SUMMARY.tsv',sep='\t')
 pair=pd.read_csv(P/'tables/PAIRED_ASSOCIATIONS.tsv',sep='\t')
 raw=pd.read_csv(P/'scores/PAIRED_PATIENT_VALUES.tsv',sep='\t')
 keys=['dataset','family','x','y','scale','adjustment','version','role']
 assert raw.patient.is_unique and len(raw)==16 and pair.patient_n.eq(16).all()
 assert not pat.duplicated(keys+['patient']).any() and not sm.duplicated(keys).any()
 assert (sm.loc[sm.patient_n<10,['P','FDR']].isna().all().all())
 reconstructed=[]
 for key,g in sec.groupby(keys+['patient']):
  vals=g.rho.dropna().to_numpy()
  if not len(vals):continue
  rho=np.tanh(sum(np.arctanh(v) for v in np.clip(vals,-.999999999,.999999999))/len(vals))
  m=np.ones(len(pat),bool)
  for col,v in zip(keys+['patient'],key):m &= pat[col].eq(v)
  assert m.sum()==1;err=abs(float(pat.loc[m,'rho'].iloc[0])-rho);assert err<1e-12;reconstructed.append(err)
 for row in sm.itertuples():
  m=np.ones(len(pat),bool)
  for col in keys:m &= pat[col].eq(getattr(row,col))
  vals=pat.loc[m,'rho'].dropna().to_numpy()
  assert len(vals)==row.patient_n
  if len(vals):
   assert abs(np.median(vals)-row.median_rho)<1e-12
   pos=int((vals>1e-12).sum());neg=int((vals< -1e-12).sum());assert pos==row.positive and neg==row.negative
   if len(vals)>=10:
    pv=stats.binomtest(pos,pos+neg,.5).pvalue if pos+neg else 1
    assert abs(pv-row.P)<1e-12
 for frame,groups in [(sm,['dataset','family','role']),(pair,['family','role'])]:
  for k,g in frame.groupby(groups):
   expected=bh(g.P.fillna(1));ok=g.P.notna()
   np.testing.assert_allclose(expected[ok],g.loc[ok,'FDR'],rtol=1e-10,atol=1e-12)
 # Independently reconstruct adjusted correlation via explicit residual-maker rather than helper QR.
 checks=[]
 for row in pair.itertuples():
  x=raw[row.x].to_numpy();col=row.y+('__no_shared' if row.version=='no_shared' else '');y=raw[col].to_numpy()
  if not np.isfinite(row.rho):continue
  cv=[]
  if row.adjustment!='raw':
   D=pd.get_dummies(raw.study,drop_first=True,dtype=float)
   cv.extend(D[c].to_numpy() for c in D)
  if row.adjustment=='study_covariate':
   cv.append(np.log1p(raw.malignant_total_counts).to_numpy() if row.family=='tumor_intracellular' else raw.C4_fraction.to_numpy())
  cols=[]
  for v in cv:
   v=stats.rankdata(v)
   if v.std()>1e-10:cols.append((v-v.mean())/v.std())
  C=np.column_stack([np.ones(len(raw))]+cols)
  H=np.eye(len(raw))-C@np.linalg.pinv(C)
  ref=float(stats.pearsonr(H@stats.rankdata(x),H@stats.rankdata(y)).statistic)
  assert abs(ref-row.rho)<1e-10
  if row.adjustment=='raw':assert abs(stats.spearmanr(x,y).pvalue-row.P)<1e-12
  checks.append(dict(family=row.family,x=row.x,y=row.y,adjustment=row.adjustment,version=row.version,error=abs(ref-row.rho)))
 out('qc/INDEPENDENT_PAIRED_CORRELATION_CHECKS.tsv',pd.DataFrame(checks))
 coverage=pd.read_csv(P/'tables/PAIRED_PROGRAM_COVERAGE.tsv',sep='\t')
 for r in coverage.itertuples():
  assert not set(r.members.split(';'))&{'TYMP','SPP1','GPNMB','LGALS3','FABP5','CD9','LPL'}
 scorechecks=pd.read_csv(P/'qc/SPATIAL_SCORE_RECONSTRUCTION.tsv',sep='\t');assert len(scorechecks)==76 and scorechecks.error.max()<1e-10
 cohort=sm[(sm.role=='primary')];out('tables/SPATIAL_PRIMARY_RESULTS.tsv',cohort)
 out('tables/PAIRED_PRIMARY_RESULTS.tsv',pair[pair.role=='primary'])
 jwrite('qc/FINAL_VALIDATION.json',dict(status='PASS',time=stamp(),spatial_sections=sec['sample'].nunique(),spatial_patients=pat[['dataset','patient']].drop_duplicates().shape[0],paired_patients=16,paired_studies=raw.study.value_counts().to_dict(),spatial_correlations=len(sec),spatial_non_estimable=int(sec.rho.isna().sum()),spatial_patient_merges_checked=len(reconstructed),spatial_cohort_rows=len(sm),paired_associations=len(pair),paired_non_estimable=int(pair.rho.isna().sum()),paired_matrix_checks=len(checks),program_reconstructions=76,inputs_unchanged=len(unique),new_models_trained=0,new_figures=0,old_RDS_modified=False))
 print('VALIDATION_PASS',flush=True)
 print(cohort[['dataset','family','x','y','patient_n','positive','negative','median_rho','P','FDR']].to_string(index=False),flush=True)
 print(pair[pair.role=='primary'][['family','x','y','patient_n','rho','P','FDR']].to_string(index=False),flush=True)
if __name__=='__main__':main()
