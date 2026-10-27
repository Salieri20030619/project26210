from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,shutil,datetime,csv
from openpyxl import load_workbook
h=Path('/mnt/BioDisk/wangmingyuan');r=h/'project26210_external_joint_v1';p=r/'phase270_SPP1_literature_clinical_revalidation_20260927_v1'
b=h/'project26210_external_joint_v1_phase64_large_immune_validation_v1'
w=load_workbook(b/'downloads_v1/CheckMate_Braun2020/Table_S1_clinical.xlsx',read_only=True,data_only=True)
it=w['Clinical_and_Immune_Data'].iter_rows(values_only=True);next(it);header=next(it);raw=pd.DataFrame(list(it),columns=header).replace({'NA':np.nan,'':np.nan});w.close()
a=pd.read_csv(p/'tables/CHECKMATE_ANALYSIS_INPUT.tsv',sep='\t')
raw=raw[raw.RNA_ID.isin(a.RNA_ID)].set_index('RNA_ID').loc[a.RNA_ID]
for col in ['PFS','PFS_CNSR','OS','OS_CNSR','Age']:
 np.testing.assert_allclose(pd.to_numeric(raw[col]),a[col],rtol=1e-10,equal_nan=True)
for col in ['SUBJID','Cohort','Arm','Sex','MSKCC']:assert raw[col].fillna('').tolist()==a[col].fillna('').tolist()
(p/'qc/CHECKMATE_RAW_CLINICAL_CHECK.json').write_text(json.dumps({'patients':len(raw),'endpoint_covariate_arm_match':True},indent=2))
old=json.loads((r/'phase269_clinical_diagnosis_story_20260927_v1/qc/SOURCE_MANIFEST.json').read_text())
for x in old:assert hashlib.sha256(Path(x['path']).read_bytes()).hexdigest()==x['sha256']
src=[b/'downloads_v1/CheckMate_Braun2020/Table_S1_clinical.xlsx',b/'downloads_v1/CheckMate_Braun2020/Table_S4_RNA.xlsx',b/'downloads_v1/JAVELIN_Renal101/Supplementary_Tables.xlsx',b/'inputs_v1/CheckMate_Braun2020/processed_v1/checkmate_rna_expression_v1.tsv.gz',b/'inputs_v1/JAVELIN_Renal101/processed_v1/javelin_clinical_v1.tsv']
manifest=[]
for f in src:manifest.append({'path':str(f),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
(p/'qc/SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
s=pd.read_csv(p/'tables/ALL_SURVIVAL_MODELS.tsv',sep='\t')
rep=[]
for score,ds,arm,ep,pap in [('Jiang5','CheckMate','NIVOLUMAB','OS',.015),('Jiang5','CheckMate','EVEROLIMUS','OS',.227),('Jiang5','CheckMate','NIVOLUMAB','PFS',.579),('Jiang5','CheckMate','EVEROLIMUS','PFS',.958)]:
 x=s[(s.score==score)&(s.dataset==ds)&(s.arm==arm)&(s.endpoint==ep)&(s.model=='median_cohort')].iloc[0]
 row=x.to_dict();row.update(paper='Jiang2025',paper_P=pap,paper_HR=np.nan,status='logrank_P_matches_to_3_decimals')
 assert round(x.logrank_P,3)==pap
 rep.append(row)
for mod in ['median_arm','median_cohort']:
 x=s[(s.score=='SPP1')&(s.dataset=='JAVELIN')&(s.arm=='Avelumab+Axitinib')&(s.model==mod)].iloc[0]
 row=x.to_dict();row.update(paper='Xu2024',paper_P=.004,paper_HR=1.61,status='direction_same_exact_effect_not_reproduced');rep.append(row)
pd.DataFrame(rep).to_csv(p/'tables/LITERATURE_REPLICATION.tsv',sep='\t',index=False)
f=p/'delivery/26210_SPP1临床文献复核与差距报告.md';txt=f.read_text()
z=s[(s.score=='Jiang5')&(s.dataset=='CheckMate')&(s.arm=='NIVOLUMAB')&(s.endpoint=='OS')&(s.model=='continuous_adjusted')].iloc[0]
txt=txt.replace('164/见完整表',f'164/{int(z.events)}')
z=s[(s.score=='Jiang5')&(s.dataset=='JAVELIN')&(s.arm=='Avelumab+Axitinib')&(s.model=='continuous_adjusted')].iloc[0]
txt=txt.replace('0.7528（以完整表为准）',f'{z.FDR:.4f}')
f.write_text(txt)
names=['ALL_SURVIVAL_MODELS.tsv','TREATMENT_INTERACTIONS.tsv','INTERACTION_DIAGNOSTICS.tsv','CONSTRUCT_CORRELATIONS.tsv','LITERATURE_REPLICATION.tsv']
for n in names:shutil.copy2(p/'tables'/n,p/'delivery'/n)
chk=json.loads((p/'qc/INDEPENDENT_CHECK_COMPLETE.json').read_text())
val={'status':'PASS_with_documented_replication_gap','completed_at':datetime.datetime.now().astimezone().isoformat(),'observed_exit_status':{'extraction':0,'models':0,'independent_check':0,'interaction_audit':0,'optimizer_check':0},'patients':{'CheckMate':311,'JAVELIN':726},'models':48,'interactions':6,'independent_check':chk,'original_phase269_inputs_unchanged':len(old),'CheckMate_paper_logrank_tests_reproduced':4,'JAVELIN_paper_exact_HR_not_reproduced':True,'source_clinical_checks':'both original workbooks match all compared endpoints, arm and covariates','manuscript_modified':False,'background_jobs':0,'optimizer_note':'one independent optimizer warning retained; max absolute score 1.75e-12 and R/Python target parameters agree','scope':'retrospective literature-guided audit, not new independent clinical validation'}
(p/'qc/FINAL_VALIDATION.json').write_text(json.dumps(val,ensure_ascii=False,indent=2))
shutil.copy2(p/'qc/FINAL_VALIDATION.json',p/'delivery/FINAL_VALIDATION.json')
print(json.dumps(val,ensure_ascii=False,indent=2));print('ADJUSTED_OS',int(s.loc[s.model_id==28,'events'].iloc[0]))
