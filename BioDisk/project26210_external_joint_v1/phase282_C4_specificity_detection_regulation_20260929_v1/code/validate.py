from common import *
def ref(x,y,cov,rank=True,minn=3):
 x=np.asarray(x,float);y=np.asarray(y,float);cv=np.asarray(cov,float)
 A=np.column_stack([stats.rankdata(x),stats.rankdata(y)]) if rank else np.column_stack([x,y])
 ranks=pd.DataFrame(cv).rank(method='average').to_numpy();sd=ranks.std(0);D=np.column_stack([np.ones(len(A)),(ranks-ranks.mean(0))/sd])
 U,sv,V=np.linalg.svd(D,full_matrices=False);res=A-U@(U.T@A)
 return float(stats.pearsonr(res[:,0],res[:,1]).statistic)
def main():
 for job in ['specificity','regulation','detection_extension']:assert (P/'logs'/f'{job}_exit.txt').read_text().strip()=='0'
 z=data();checks=[];mx=0.;nc=0
 for file in ['SCOPE_PATIENT_EFFECTS.tsv','STATE_GRID_PATIENT_EFFECTS.tsv']:
  d=read(P/'tables'/file,dtype={'cluster':str})
  for r in d[d.effect.notna()].itertuples():
   h=z[z.patient==r.patient]
   if file.startswith('STATE'):h=h[h.candidate==r.cluster];other=False;gene='SPP1'
   else:h=h[h.candidate.eq('4') if r.scope=='C4' else h.candidate.ne('4')];other=r.scope=='OtherMacro';gene=r.gene
   e=ref(h.TYMP,h[gene],covariates(h,r.model,other));assert abs(e-r.effect)<1e-10,(str(r),e,r.effect,abs(e-r.effect));mx=max(mx,abs(e-r.effect));nc+=1
 checks.append(dict(check='scope_and_state_rhos_SVD',n=nc,max_error=mx))
 d=read(P/'tables/DECOMPOSITION_PATIENT_EFFECTS.tsv');mx=0.;nc=0
 for r in d[d.effect.notna()].itertuples():
  h=z[z.patient==r.patient];h=h[h.candidate.eq('4') if r.scope=='C4' else h.candidate.ne('4')]
  binary=r.part=='co_detection'
  if not binary:h=h[h.det_TYMP&h.det_SPP1]
  e=ref(h.det_TYMP if binary else h.TYMP,h.det_SPP1 if binary else h.SPP1,covariates(h,r.model,r.scope=='OtherMacro'),rank=not binary)
  assert abs(e-r.effect)<1e-10,(str(r),e,r.effect,abs(e-r.effect));mx=max(mx,abs(e-r.effect));nc+=1
 checks.append(dict(check='detection_and_positive_intensity_SVD',n=nc,max_error=mx))
 tab2=read(P/'tables/DETECTION_2X2_COUNTS.tsv')
 for r in tab2.itertuples():
  h=z[z.patient==r.patient];h=h[h.candidate.eq('4') if r.scope=='C4' else h.candidate.ne('4')]
  counts=pd.crosstab(h.det_TYMP,h.det_SPP1).reindex(index=[False,True],columns=[False,True],fill_value=0).fillna(0)
  assert counts.at[True,True]==r.n11 and counts.at[True,False]==r.n10 and counts.at[False,True]==r.n01 and counts.at[False,False]==r.n00
 checks.append(dict(check='all_2x2_detection_tables',n=len(tab2)))
 d=read(P/'tables/TYMP_INTENSITY_SPP1_DETECTION_PATIENT_EFFECTS.tsv');mx=0.;nc=0
 for r in d[d.effect.notna()].itertuples():
  h=z[z.patient==r.patient];h=h[h.candidate.eq('4') if r.scope=='C4' else h.candidate.ne('4')]
  e=ref(stats.rankdata(h.TYMP),h.det_SPP1,covariates(h,r.model,r.scope=='OtherMacro'),rank=False);assert abs(e-r.effect)<1e-10,(str(r),e,r.effect,abs(e-r.effect));mx=max(mx,abs(e-r.effect));nc+=1
 checks.append(dict(check='continuous_to_detection_SVD',n=nc,max_error=mx))
 c=read(P/'tables/C4_OTHER_PAIRED_PATIENT_EFFECTS.tsv')
 assert not c[['patient','gene','model']].duplicated().any()
 assert np.allclose(c.effect,c.effect_C4-c.effect_Other,equal_nan=True)
 checks.append(dict(check='paired_same_patient_difference',n=len(c)))
 scores=read(P/'inputs/REGULON_CELL_SCORES.tsv.gz',float_precision='round_trip');zz=z.merge(scores,on=['cell_id','patient'],validate='one_to_one')
 rtab=read(P/'tables/REGULON_PATIENT_EFFECTS.tsv');mx=0.;nc=0
 for (patient,scope,model),d in rtab[rtab.effect.notna()].groupby(['patient','scope','model']):
  h=zz[zz.patient==patient];h=h[h.candidate.eq('4') if scope=='C4' else h.candidate.ne('4')];cov=covariates(h,model,scope=='OtherMacro')
  for r in d.itertuples():
   e=ref(h[r.regulon],h[r.gene],cov);assert abs(e-r.effect)<1e-10,(str(r),e,r.effect,abs(e-r.effect));mx=max(mx,abs(e-r.effect));nc+=1
 checks.append(dict(check='all_regulon_rhos_SVD',n=nc,max_error=mx))
 d=read(P/'tables/CONDITIONAL_REMAINING_PATIENT_EFFECTS.tsv');mx=0.;nc=0
 for (patient,scope),g in d[d.effect.notna()].groupby(['patient','scope']):
  h=zz[zz.patient==patient];h=h[h.candidate.eq('4') if scope=='C4' else h.candidate.ne('4')]
  for r in g.itertuples():
   e=ref(h.TYMP,h.SPP1,covariates(h,'context',scope=='OtherMacro',h[r.regulon]));assert abs(e-r.effect)<1e-10,(str(r),e,r.effect,abs(e-r.effect));mx=max(mx,abs(e-r.effect));nc+=1
 checks.append(dict(check='conditional_remaining_SVD',n=nc,max_error=mx))
 # Independently reconstruct every target-program score using pandas standardization.
 cfg=json.loads((P/'resources/REGULON_PROGRAMS_FROZEN.json').read_text());sets=cfg['programs'];excluded=set(cfg['excluded']);assert all(not(set(gs)&excluded) for gs in sets.values())
 meta=read(M/'inputs/MACRO_METADATA.tsv.gz',index_col=0,dtype={'candidate':str});genes=read(M/'inputs/GENES.tsv').gene.tolist();gi={g:i for i,g in enumerate(genes)};X=sparse.load_npz(M/'inputs/MACRO_COUNTS.npz').tocsr();total=np.asarray(X.sum(1)).ravel()
 union=sorted(set(sum(list(sets.values()),[])));lookup={g:i for i,g in enumerate(union)};mx=0.;nvals=0
 for patient,h in z.groupby('patient'):
  ids=meta.index.get_indexer(h.cell_id);raw=X[ids][:,[gi[g] for g in union]].toarray().astype(float);log=np.log1p(raw/total[ids,None]*1e4)
  saved=scores.set_index('cell_id').loc[h.cell_id]
  for reg,gs in sets.items():
   v=pd.DataFrame(log[:,[lookup[g] for g in gs]]);sd=v.std(ddof=0).replace(0,1);refscore=((v-v.mean())/sd).mean(axis=1).to_numpy();err=np.max(np.abs(refscore-saved[reg].to_numpy()));assert err<1e-10;mx=max(mx,float(err));nvals+=len(h)
 checks.append(dict(check='all_disjoint_program_scores_reconstructed',values=nvals,max_error=mx))
 # Full summary family, P and BH checks, including non-estimable members in family.
 pairs=[
 ('SCOPE_SUMMARY.tsv','SCOPE_PATIENT_EFFECTS.tsv',['scope','gene','model'],['scope']),
 ('C4_OTHER_PAIRED_SUMMARY.tsv','C4_OTHER_PAIRED_PATIENT_EFFECTS.tsv',['family','gene','model'],['family']),
 ('STATE_GRID_SUMMARY.tsv','STATE_GRID_PATIENT_EFFECTS.tsv',['family','cluster','model'],['family']),
 ('DECOMPOSITION_SUMMARY.tsv','DECOMPOSITION_PATIENT_EFFECTS.tsv',['family','scope','part','model'],['family']),
 ('REGULON_ASSOCIATION_SUMMARY.tsv','REGULON_PATIENT_EFFECTS.tsv',['scope','regulon','gene','model'],['scope']),
 ('CONDITIONAL_REMAINING_SUMMARY.tsv','CONDITIONAL_REMAINING_PATIENT_EFFECTS.tsv',['family','scope','regulon'],['family']),
 ('TYMP_INTENSITY_SPP1_DETECTION_SUMMARY.tsv','TYMP_INTENSITY_SPP1_DETECTION_PATIENT_EFFECTS.tsv',['family','scope','model'],['family'])]
 for sf,ef,keys,fam in pairs:
  s=read(P/'tables'/sf,dtype={'cluster':str});e=read(P/'tables'/ef,dtype={'cluster':str})
  for row in s.to_dict('records'):
   mask=np.ones(len(e),bool)
   for key in keys:mask &= e[key].eq(row[key]).to_numpy()
   v=e.loc[mask,'effect'].dropna().to_numpy();assert len(v)==row['n']
   assert int((v>0).sum())==row['positive'] and int((v<0).sum())==row['negative']
   if len(v):assert abs(np.median(v)-row['median'])<1e-12
   if len(v)>=10:
    pp=stats.binomtest(int((v>0).sum()),int((v!=0).sum()),.5).pvalue if np.any(v!=0) else 1.
    assert abs(pp-row['P'])<1e-12
   else:assert pd.isna(row['P'])
  for _,g in s.groupby(fam):
   qq=stats.false_discovery_control(g.P.fillna(1).to_numpy());ok=g.P.notna().to_numpy();assert np.allclose(qq[ok],g.FDR.to_numpy()[ok],atol=1e-12)
  checks.append(dict(check='sign_BH_'+sf,n=len(s)))
 # Phase280 C4 effects preserved exactly, despite new comparative families.
 old=read(PRE/'tables/WITHIN_MACRO_PATIENT_EFFECTS.tsv');new=read(P/'tables/SCOPE_PATIENT_EFFECTS.tsv')
 old=old[(old.group=='C4')&(old.model=='adjusted')&old.endpoint.isin(GENES)]
 j=old.merge(new[(new.scope=='C4')&(new.model=='core')],left_on=['patient','endpoint'],right_on=['patient','gene'],suffixes=('_old','_new'))
 assert np.allclose(j.effect_old,j.effect_new,equal_nan=True,atol=1e-10);checks.append(dict(check='Phase280_C4_effects_preserved',n=len(j)))
 for name in ['SPECIFICITY_INPUTS.json','REGULATION_INPUTS.json']:
  for rec in json.loads((P/'qc'/name).read_text()):assert sha(rec['path'])==rec['sha256']
 assert sha(P/'PROTOCOL_BEFORE_TESTING.md')==json.loads((P/'FREEZE.json').read_text())['protocol_sha256']
 dump(P/'qc/FINAL_VALIDATION.json',dict(status='PASS',checks=checks,all_input_hashes_unchanged=True,no_new_figures=True,software=dict(python=sys.version,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__)))
 print('FINAL_VALIDATION_PASS',checks,flush=True)
if __name__=='__main__':main()
