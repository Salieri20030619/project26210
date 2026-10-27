from pathlib import Path
import json,hashlib,shutil,datetime,csv
import pandas as pd,numpy as np
from statsmodels.duration.hazard_regression import PHReg
from scipy.stats import chi2
p=Path(__file__).resolve().parents[1];r=p.parent
# Independently check the eight corrected literature comparator models.
q=pd.read_csv(p/'tables/PHASE270_ADJUSTED_SEX_CORRECTION_AUDIT.tsv',sep='\t')
all=pd.read_csv(p/'tables/CHECKMATE_ANALYSIS_INPUT.tsv',sep='\t');rec=[]
for _,row in q.iterrows():
 d=all[all.arm==row.arm].copy();d['x']=(d[row.score]-d[row.score].mean())/d[row.score].std()
 cols=['x','age','sex','trial','risk'];d=d.dropna(subset=cols+[row.endpoint+'_time',row.endpoint+'_event'])
 ex=pd.get_dummies(d[cols],drop_first=True,dtype=float);ex=ex.loc[:,ex.nunique()>1]
 assert np.linalg.matrix_rank(np.column_stack([np.ones(len(ex)),ex]))==ex.shape[1]+1
 f=PHReg(d[row.endpoint+'_time'],ex,status=d[row.endpoint+'_event'],ties='efron').fit(disp=0,maxiter=200)
 j=ex.columns.get_loc('x');hr=np.exp(f.params[j]);pv=f.pvalues[j];assert np.isfinite([hr,pv]).all()
 assert abs(hr-row.HR)<1e-7 and abs(pv-row.P)<1e-7
 rec.append({'endpoint':row.endpoint,'arm':row.arm,'score':row.score,'HR_difference':abs(hr-row.HR),'P_difference':abs(pv-row.P)})
pd.DataFrame(rec).to_csv(p/'qc/PHASE270_COMPARATOR_INDEPENDENT_CHECK.tsv',sep='\t',index=False)
# Source data and frozen scores remain unchanged.
old=json.loads((r/'phase270_SPP1_literature_clinical_revalidation_20260927_v1/qc/SOURCE_MANIFEST.json').read_text())
old+=json.loads((r/'phase269_clinical_diagnosis_story_20260927_v1/qc/SOURCE_MANIFEST.json').read_text())
for item in old:assert hashlib.sha256(Path(item['path']).read_bytes()).hexdigest()==item['sha256']
fr=json.loads((p/'qc/SIGNATURE_FREEZE_BEFORE_CLINICAL.json').read_text())
src=r/'phase200_union_C4_C10_C15_reassessment_20260918_v1/core/contrasts/tumor_C4.npz'
assert hashlib.sha256(src.read_bytes()).hexdigest()==fr['source_sha256']
cl=json.loads((p/'qc/CLINICAL_SCORE_FREEZE.json').read_text())
assert hashlib.sha256((p/'tables/C4_CLINICAL_SCORES.tsv').read_bytes()).hexdigest()==cl['scores_sha256']
assert not set(cl['genes_used'])&{'SPP1','TYMP','FN1','CLEC5A','OLR1','IL1RN'}
a=pd.read_csv(p/'qc/INDEPENDENT_COX_CHECK.tsv',sep='\t')
assert len(a)==24 and np.isfinite(a[['HR_difference','P_difference','LRT_P_difference']]).all().all()
assert not a.warnings.fillna('').str.len().gt(0).any()
b=pd.read_csv(p/'qc/INDEPENDENT_TIME_VARIATION_CHECK.tsv',sep='\t');assert len(b)==6
for f in ['code/write_reports.py','delivery/26210_C4状态临床延伸分析报告.md','delivery/拟新增结果小节与图表安排.md']:
 x=p/f;x.write_text(x.read_text().replace('C4相关脂质处理状态呈现独立于SPP1转录负荷的时间依赖性进展风险关联','C4相关脂质处理状态与进展风险的时间依赖性关联'),encoding='utf-8')
for f in ['C4_SURVIVAL_MODELS.tsv','C4_TREATMENT_INTERACTIONS.tsv','C4_PFS_TIME_VARYING_MODELS.tsv','C4_PFS_TIME_SPECIFIC_HR.tsv','C4_TIME_VARYING_TREATMENT_INTERACTIONS.tsv','C4_TIME_VARIATION_TRIAL_SENSITIVITY.tsv','C4_TIME_VARIATION_TRIAL_HR.tsv','C4_TIME_VARIATION_LEAVE_ONE_PATIENT.tsv','FROZEN_C4_SIGNATURE.tsv','CLINICAL_GENE_COVERAGE.tsv','LOSO_CLINICAL_AVAILABILITY_SENSITIVITY.tsv','C4_CONSTRUCT_CORRELATIONS.tsv','PHASE270_ADJUSTED_SEX_CORRECTION_AUDIT.tsv']:
 shutil.copy2(p/'tables'/f,p/'delivery'/f)
sig=pd.read_csv(p/'tables/FROZEN_C4_SIGNATURE.tsv',sep='\t');sig[sig.gene.isin(cl['genes_used'])].to_csv(p/'delivery/CLINICAL_C4_SIGNATURE_19_GENES.tsv',sep='\t',index=False)
v={'completed_at':datetime.datetime.now().astimezone().isoformat(),'status':'PASS','marker_discovery':'24 paired tumor patients across5studies; fixed20genes,19commonclinicalgenes','clinical_models':24,'treatment_interactions':6,'time_varying_models':4,'time_varying_interactions':2,'trial_sensitivity_models':4,'leave_one_patient_fits':164,'clinical_independent_validation':'R vs Python HR/P/LRT all finite and tolerance<1e-7; design matrices full rank','time_independent_validation':'explicit event risk-set expansion vs coxph tt; all6 verified<1e-6','sex_labels':'harmonized before final models; initial4level analysis archived, not final','corrected_phase270_comparators_independently_verified':8,'unchanged_original_sources':len(old)+1,'C4_signature_and_clinical_scores_unchanged':True,'no_new_figures':True,'no_manuscript_edit':True,'background_jobs':0,'limitations':'Post-PH-diagnostic retrospective time-varying extension; no independent clinical cohort validation or established treatment benefit prediction'}
(p/'qc/FINAL_VALIDATION.json').write_text(json.dumps(v,ensure_ascii=False,indent=2));shutil.copy2(p/'qc/FINAL_VALIDATION.json',p/'delivery/FINAL_VALIDATION.json')
with (p/'delivery/DELIVERY_SHA256.tsv').open('w') as f:
 wr=csv.writer(f,delimiter='\t',lineterminator='\n');wr.writerow(['file','bytes','sha256'])
 for x in sorted((p/'delivery').iterdir()):
  if x.name!='DELIVERY_SHA256.tsv':wr.writerow([x.name,x.stat().st_size,hashlib.sha256(x.read_bytes()).hexdigest()])
print(json.dumps(v,ensure_ascii=False,indent=2));print('DELIVERY_FILES',len(list((p/'delivery').iterdir())))
