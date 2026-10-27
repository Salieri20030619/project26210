from pathlib import Path
import json,hashlib,datetime,traceback
import numpy as np,pandas as pd
from scipy.stats import rankdata
P=Path(__file__).resolve().parents[1];ROOT=P.parent;G=ROOT/'phase237_SpaMTP_registration_20260923_v1';S=ROOT/'phase235_SpaMTP_TJRCC_pilot_20260921_v1';F=ROOT/'phase234_part4_replan_execution_20260921_v1'
TARGETS=['TYMP_RNA','V1_SPP1_q05','V3_SPP1_macro_share','polyamine_mta_context','V4_SPP1_program'];PRIMARY=['TYMP_RNA','V1_SPP1_q05','V3_SPP1_macro_share']
def sha(q):return hashlib.sha256(q.read_bytes()).hexdigest()
def snapshot(q):return dict(path=str(q),bytes=q.stat().st_size,mtime_ns=q.stat().st_mtime_ns,sha256=sha(q))
def now():return datetime.datetime.now().astimezone().isoformat()
def rho(x,y):
 a=rankdata(x,method='average');b=rankdata(y,method='average');a-=a.mean();b-=b.mean();return float(a@b/np.sqrt((a@a)*(b@b)))
def direction(v):return 'positive' if v>1e-12 else 'negative' if v< -1e-12 else 'zero'
def main():
 assert json.loads((G/'FINAL_VALIDATION.json').read_text())['status']=='PASS'
 mapping=pd.read_csv(G/'transforms/CANDIDATE_ST_TO_SM_MAPS.tsv.gz',sep='\t');candidates=json.loads((G/'transforms/RETAINED_CANDIDATES.json').read_text());peaks=pd.read_csv(S/'tables/TARGETED_MASS_CANDIDATES_PASS5PPM.tsv',sep='\t')
 inp=[G/'transforms/CANDIDATE_ST_TO_SM_MAPS.tsv.gz',G/'transforms/RETAINED_CANDIDATES.json',S/'tables/TARGETED_MASS_CANDIDATES_PASS5PPM.tsv',F/'PROGRAM_DEFINITIONS.json']
 for sample in ['R114_T','Y12_T']:
  inp.append(F/'inputs'/(sample+'_INPUTS.tsv.gz'))
  for pol in ['pos','neg']:inp.append(S/'results/matrices'/(sample.replace('_','-')+'-'+pol+'.npz'))
 before=[snapshot(q) for q in inp];(P/'audit/INPUT_MANIFEST.json').write_text(json.dumps(before,indent=2))
 identities=[];effects=[];regions=[];covers=[];mapping_audits=[];checks=[];peakqc=[]
 for sample,patient in [('R114_T','TJRCC:R114'),('Y12_T','TJRCC:Y12')]:
  d=pd.read_csv(F/'inputs'/(sample+'_INPUTS.tsv.gz'),sep='\t');assert d.spot_id.is_unique and d['patient'].nunique()==1 and d.patient.iloc[0]==patient and d.dataset.eq('TJRCC').all();assert np.isfinite(d[TARGETS].to_numpy()).all();assert np.allclose(d.V4_SPP1_program,d.spp1_state_no_spp1)
  identities.append(dict(sample=sample,registered_patient=patient,dataset='TJRCC',tissue='tumor',ST_QC_spots=len(d),primary_region_labels=int(d.block_4.nunique()),SM_acquisitions=2,selection='Phase235 smallest complete paired-polarity archives selected before inspecting molecular associations'))
  for pol in ['neg','pos']:
   label=sample.replace('_','-')+'-'+pol;z=np.load(S/'results/matrices'/(label+'.npz'));X=z['intensity'];tic=z['tic'];xy=z['xy'];cs=sorted([c for c in candidates if c['label']==label],key=lambda c:c['rank']);assert cs and cs[0]['rank']==1
   frames={};sets=[]
   for c in cs:
    m=mapping[mapping.candidate_id==c['candidate_id']].copy();assert m.spot_id.is_unique
    frame=d.merge(m,on='spot_id',how='left',validate='one_to_one',suffixes=('','_mapping'));assert frame.candidate_id.notna().all() and frame.QC_pass.all();assert frame['sample'].eq(frame.sample_mapping).all()
    ids=frame.SM_pixel_index.to_numpy(dtype=int);assert np.all(ids>=0) and np.all(ids<len(X));assert np.allclose(xy[ids],frame[['SM_nearest_x_one_based','SM_nearest_y_one_based']].to_numpy());assert (tic[ids]>0).all();assert np.allclose(frame.x_hires_px,frame.ST_x_hires) and np.allclose(frame.y_hires_px,frame.ST_y_hires)
    included=frame.included.astype(bool);frames[c['candidate_id']]=frame;sets.append(set(frame.loc[included,'spot_id']))
    mapping_audits.append(dict(candidate_id=c['candidate_id'],sample=sample,polarity=pol,QC_barcodes=len(d),matched_barcodes=len(frame),positive_TIC_pixels_only=True,coordinates_and_pixel_indices_match=True))
   shared=set.intersection(*sets);assert shared
   for c in cs:
    frame=frames[c['candidate_id']];ids=frame.SM_pixel_index.to_numpy(dtype=int)
    for mode in ['candidate_specific','common_coverage']:
     mask=frame.included.astype(bool).to_numpy()
     if mode=='common_coverage':mask &= frame.spot_id.isin(shared).to_numpy()
     gframe=frame.loc[mask].copy();coverage=len(gframe)/len(d);covers.append(dict(sample=sample,patient=patient,polarity=pol,candidate_id=c['candidate_id'],rank=c['rank'],coverage_mode=mode,mapped_QC_spots=len(gframe),all_QC_spots=len(d),coverage=coverage,unique_SM_pixels=gframe.SM_pixel_index.nunique(),common_QC_spots=len(shared)))
     selected=peaks[(peaks['sample']==sample)&(peaks.polarity==pol)]
     for peak in selected.itertuples():
      j=int(peak.feature_index);assert abs(z['mz'][j]-peak.observed_mz)<1e-8
      values=X[:,j].astype(float)/np.where(tic>0,tic,1)*X.shape[1];assert np.isfinite(values).all() and (values>=0).all()
      work=gframe[['spot_id','block_4','SM_pixel_index']+TARGETS].copy();work['SM_quantity']=values[work.SM_pixel_index.to_numpy(dtype=int)];work['SM_detected']=(X[work.SM_pixel_index.to_numpy(dtype=int),j]>0).astype(float)
      grouped=work.groupby('block_4',sort=True)[['SM_quantity','SM_detected']+TARGETS].mean();n=work.groupby('block_4').size();u=work.groupby('block_4').SM_pixel_index.nunique();kept=n>=10;grouped=grouped.loc[kept]
      peakqc.append(dict(sample=sample,polarity=pol,candidate_id=c['candidate_id'],rank=c['rank'],coverage_mode=mode,molecule=peak.molecule,observed_mz=peak.observed_mz,mapped_spots=len(work),detected_spots=int(work.SM_detected.sum()),detected_spot_fraction=float(work.SM_detected.mean()),eligible_regions=len(grouped),regions_with_detected_peak=int((grouped.SM_detected>0).sum())))
      for block,r in grouped.iterrows():
       regions.append(dict(sample=sample,patient=patient,polarity=pol,candidate_id=c['candidate_id'],rank=c['rank'],coverage_mode=mode,molecule=peak.molecule,observed_mz=peak.observed_mz,block_4=int(block),mapped_spots=int(n.loc[block]),unique_SM_pixels=int(u.loc[block]),**r.to_dict()))
      for target in TARGETS:
       x=grouped.SM_quantity.to_numpy();y=grouped[target].to_numpy();reason=[]
       if coverage<.5:reason.append('coverage_below_frozen_rule')
       if len(work)<50:reason.append('too_few_mapped_spots')
       if len(grouped)<8:reason.append('too_few_eligible_regions')
       if len(np.unique(x))<3:reason.append('insufficient_peak_variation')
       if len(np.unique(y))<3:reason.append('insufficient_RNA_variation')
       r=rho(x,y) if not reason else np.nan
       if not reason:
        independent=pd.Series(x).rank(method='average').corr(pd.Series(y).rank(method='average'),method='pearson');checks.append(abs(r-independent));assert abs(r-independent)<1e-12
       effects.append(dict(sample=sample,patient=patient,polarity=pol,candidate_id=c['candidate_id'],rank=c['rank'],coverage_mode=mode,molecule=peak.molecule,observed_mz=peak.observed_mz,mass_identity='mass_compatible_candidate',RNA_target=target,analysis_role='primary' if peak.molecule in ['Spermidine','Spermine'] and target in PRIMARY else 'supplementary',n_regions=len(grouped),mapped_spots=len(work),QC_coverage=coverage,rho=r,direction=direction(r) if not reason else 'unevaluable',reason=';'.join(reason),biological_patients_in_effect=1,anatomically_validated_registration=False))
  print('CASE_COMPLETE',sample,flush=True)
 e=pd.DataFrame(effects);reg=pd.DataFrame(regions);cov=pd.DataFrame(covers);pd.DataFrame(identities).to_csv(P/'tables/CASE_IDENTITY.tsv',sep='\t',index=False);e.to_csv(P/'tables/ALL_REGIONAL_EFFECTS.tsv',sep='\t',index=False);reg.to_csv(P/'tables/REGIONAL_VALUES.tsv.gz',sep='\t',index=False);cov.to_csv(P/'tables/MAPPING_COVERAGE.tsv',sep='\t',index=False);pd.DataFrame(peakqc).to_csv(P/'tables/PEAK_COVERAGE.tsv',sep='\t',index=False);pd.DataFrame(mapping_audits).to_csv(P/'audit/MAPPING_INTEGRITY.tsv',sep='\t',index=False)
 summary=[]
 keys=['sample','patient','polarity','molecule','observed_mz','RNA_target','analysis_role','coverage_mode']
 for key,g in e.groupby(keys,sort=False):
  c=dict(zip(keys,key));g=g.sort_values('rank');main=g[g['rank']==1].iloc[0];v=g.rho.dropna().to_numpy();complete=len(v)==len(g);pos=int((v>1e-12).sum());neg=int((v< -1e-12).sum());zero=int((np.abs(v)<=1e-12).sum());n=len(g)
  stability='single_candidate_orientation_sensitivity_not_testable' if n==1 else 'incomplete_estimability' if not complete else 'direction_changes_across_candidates' if pos and neg else 'same_positive_direction' if pos==n else 'same_negative_direction' if neg==n else 'zero_or_mixed_with_zero'
  summary.append(dict(**c,primary_rho=main.rho,primary_direction=main.direction,primary_regions=main.n_regions,primary_mapped_spots=main.mapped_spots,candidates=n,estimable_candidates=len(v),rho_min=float(v.min()) if len(v) else np.nan,rho_max=float(v.max()) if len(v) else np.nan,positive_candidates=pos,negative_candidates=neg,zero_candidates=zero,orientation_stability=stability,range_is_confidence_interval=False))
 su=pd.DataFrame(summary);su.to_csv(P/'tables/DIRECTION_STABILITY.tsv',sep='\t',index=False);su[(su.coverage_mode=='candidate_specific')&(su.analysis_role=='primary')].to_csv(P/'tables/PRIMARY_FINDINGS.tsv',sep='\t',index=False)
 cross=[]
 for (mol,target,mode),g in su[su.analysis_role=='primary'].groupby(['molecule','RNA_target','coverage_mode']):
  assert len(g)==2;v=g.primary_rho.to_numpy();cross.append(dict(molecule=mol,RNA_target=target,coverage_mode=mode,biological_patients=2,R114_rho=float(g[g['sample']=='R114_T'].primary_rho.iloc[0]),Y12_rho=float(g[g['sample']=='Y12_T'].primary_rho.iloc[0]),primary_case_direction='both_positive' if (v>0).all() else 'both_negative' if (v<0).all() else 'opposite_or_zero',all_candidate_effects_same_positive=bool((g.positive_candidates==g.candidates).all()),all_candidate_effects_same_negative=bool((g.negative_candidates==g.candidates).all()),Y12_retained_orientation_count=int(g[g['sample']=='Y12_T'].candidates.iloc[0]),pooled_effect_computed=False))
 pd.DataFrame(cross).to_csv(P/'tables/CROSS_CASE_DESCRIPTIVE_COMPARISON.tsv',sep='\t',index=False)
 # Compare prior preferred candidates descriptively without transferring old coefficients.
 old=pd.read_csv(S/'tables/CONDITIONAL_MAPPING_EFFECTS.tsv',sep='\t');comparison=[]
 for r in su[su.coverage_mode=='candidate_specific'].itertuples():
  prev=old[(old['sample']==r.sample)&(old.polarity==r.polarity)&(old.molecule==r.molecule)&(old.RNA_target==r.RNA_target)].rho.dropna()
  comparison.append(dict(sample=r.sample,polarity=r.polarity,molecule=r.molecule,RNA_target=r.RNA_target,old_phase235_mapping_effect_min=float(prev.min()) if len(prev) else np.nan,old_phase235_mapping_effect_max=float(prev.max()) if len(prev) else np.nan,new_primary_rho=r.primary_rho,new_mapping_effect_min=r.rho_min,new_mapping_effect_max=r.rho_max,note='Definitions preserved; geometry and inclusion changed; ranges are not statistical confidence intervals'))
 pd.DataFrame(comparison).to_csv(P/'tables/PHASE235_CONTEXT_COMPARISON.tsv',sep='\t',index=False)
 for item in before:
  q=Path(item['path']);assert q.stat().st_size==item['bytes'] and q.stat().st_mtime_ns==item['mtime_ns'] and sha(q)==item['sha256']
 assert len(e)==210 and len(su)==100 and len(e[(e['rank']==1)&(e.coverage_mode=='candidate_specific')])==50
 assert not any(c.lower() in ['p','pvalue','p_value','fdr','qvalue'] for c in e.columns)
 qa=dict(status='PASS',finished_at=now(),patients=2,patient_ids=['TJRCC:R114','TJRCC:Y12'],candidate_poses=8,coverage_modes=2,all_effect_records=len(e),estimable_effects=int(e.rho.notna().sum()),unestimable_effects=int(e.rho.isna().sum()),preferred_pose_effects=50,main_polyamine_TYMP_SPP1_effects=12,regional_value_rows=len(reg),independent_rank_checks=len(checks),max_rank_correlation_discrepancy=max(checks),input_files_unchanged=len(before),old_RDS_modified=False,new_P_values=0,new_FDR_values=0,new_figures=0,validated_registration=False)
 (P/'ANALYSIS_QA.json').write_text(json.dumps(qa,indent=2));(P/'ANALYSIS_COMPLETE').write_text(now()+'\n');print(su[(su.coverage_mode=='candidate_specific')&(su.analysis_role=='primary')].to_string(index=False));print('QA',json.dumps(qa))
code=1
if __name__=='__main__':
 try:main();code=0
 except BaseException:(P/'logs/ERROR.txt').write_text(traceback.format_exc());print(traceback.format_exc(),flush=True)
 finally:(P/'EXIT_STATUS').write_text(str(code)+'\n')
 raise SystemExit(code)
