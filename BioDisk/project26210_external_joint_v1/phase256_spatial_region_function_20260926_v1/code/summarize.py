from pathlib import Path
import sys,json,re
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'code'))
from extract import gate
gate('summarize_results')
assert json.loads((P/'qc/FINAL_VALIDATION.json').read_text())['status']=='PASS'
d=pd.read_csv(P/'tables/ALL_FUNCTIONAL_EVIDENCE.tsv.gz',sep='\t')
focus=json.loads((P/'resources/FOCUS_QUESTIONS_BEFORE_RESULT_REVIEW.json').read_text())
def name(term):return term.split('__',1)[-1]
def selected(row):
 if row.library=='Hallmark':return True
 names=focus.get(row.library,[])
 s=name(row.term).lower()
 return any(s==n.lower() or s.startswith(n.lower()+' ') for n in names)
d['display_name']=d.term.map(name)
f=d[d.apply(selected,axis=1)].copy()
f.to_csv(P/'tables/FOCUS_QUESTIONS_ALL_RESULTS.tsv',sep='\t',index=False)
# Explicit reference GO terms including any missing/non-estimable rows.
cov=pd.read_csv(P/'tables/GENESET_COVERAGE.tsv.gz',sep='\t')
cov['display_name']=cov.term.map(name)
ref=cov[(cov.library=='GO_BP')&cov.display_name.str.lower().isin([s.lower() for s in focus['GO_BP']])]
ref.merge(d,on=['dataset','model','version','library','term','display_name'],how='left',suffixes=('_coverage','')).to_csv(P/'tables/REFERENCE_FIG2G_ALL_TEN_TERMS.tsv',sep='\t',index=False)
primary=d[(d.model=='absolute')&(d.version=='exclude_definition_genes')].copy()
primary.to_csv(P/'tables/PRIMARY_FULL_LIBRARY_RESULTS.tsv.gz',sep='\t',index=False)
summ=[]
for key,g in d.groupby(['dataset','model','version','library']):
 summ.append(dict(zip(['dataset','model','version','library'],key),terms=len(g),patient_n=int(g.patient_n.iloc[0]),patient_FDR_positive=int(((g.FDR_patient_sign<.05)&(g.median_Z_difference>0)).sum()),patient_FDR_negative=int(((g.FDR_patient_sign<.05)&(g.median_Z_difference<0)).sum()),GSEA_FDR_positive=int(((g.FDR_GSEA<.05)&(g.NES>0)).sum()),GSEA_FDR_negative=int(((g.FDR_GSEA<.05)&(g.NES<0)).sum()),GSEA_NA=int(g.NES.isna().sum())))
pd.DataFrame(summ).to_csv(P/'tables/RESULT_COUNTS.tsv',sep='\t',index=False)
base=primary[primary.dataset=='GSE175540'];other=primary[primary.dataset=='TJRCC']
cols=['library','term','display_name','patient_n','median_Z_difference','positive','negative','P_patient_sign','FDR_patient_sign','NES','P_GSEA','FDR_GSEA']
cross=base[cols].merge(other[cols],on=['library','term','display_name'],suffixes=('_GSE175540','_TJRCC'))
cross['same_patient_direction']=np.sign(cross.median_Z_difference_GSE175540)==np.sign(cross.median_Z_difference_TJRCC)
cross['both_patient_FDR']=cross.FDR_patient_sign_GSE175540.lt(.05)&cross.FDR_patient_sign_TJRCC.lt(.05)
cross['both_GSEA_FDR']=cross.FDR_GSEA_GSE175540.lt(.05)&cross.FDR_GSEA_TJRCC.lt(.05)
cross['GSEA_patient_sign_agree']=(np.sign(cross.NES_GSE175540)==np.sign(cross.median_Z_difference_GSE175540))&(np.sign(cross.NES_TJRCC)==np.sign(cross.median_Z_difference_TJRCC))
for model in ['adjusted','macro_share']:
 for ds in ['GSE175540','TJRCC']:
  sub=d[(d.dataset==ds)&(d.model==model)&(d.version=='exclude_definition_genes')][['library','term','median_Z_difference','positive','negative','FDR_patient_sign','NES','FDR_GSEA']]
  cross=cross.merge(sub.rename(columns={c:c+'_'+model+'_'+ds for c in sub if c not in ['library','term']}),on=['library','term'],how='left')
 cross[model+'_direction_retained']=(np.sign(cross['median_Z_difference_'+model+'_GSE175540'])==np.sign(cross.median_Z_difference_GSE175540))&(np.sign(cross['median_Z_difference_'+model+'_TJRCC'])==np.sign(cross.median_Z_difference_TJRCC))
cross['max_primary_patient_FDR']=cross[['FDR_patient_sign_GSE175540','FDR_patient_sign_TJRCC']].max(1)
cross['max_primary_GSEA_FDR']=cross[['FDR_GSEA_GSE175540','FDR_GSEA_TJRCC']].max(1)
cross=cross.sort_values(['same_patient_direction','both_patient_FDR','both_GSEA_FDR','max_primary_patient_FDR','max_primary_GSEA_FDR'],ascending=[False,False,False,True,True])
cross.to_csv(P/'tables/CROSS_COHORT_AND_SENSITIVITY.tsv',sep='\t',index=False)
print('RESULT_COUNTS_PRIMARY')
print(pd.DataFrame(summ).query("model=='absolute' and version=='exclude_definition_genes'").to_string(index=False))
print('REFERENCE_GO_PRIMARY')
q=f[(f.library=='GO_BP')&(f.model=='absolute')&(f.version=='exclude_definition_genes')]
print(q[['dataset','display_name','patient_n','median_Z_difference','positive','negative','P_patient_sign','FDR_patient_sign','NES','P_GSEA','FDR_GSEA']].to_string(index=False))
print('HALLMARK_PRIMARY')
q=f[(f.library=='Hallmark')&(f.model=='absolute')&(f.version=='exclude_definition_genes')&~f.dataset.str.contains('GLMF')]
print(q[['dataset','display_name','median_Z_difference','positive','negative','FDR_patient_sign','NES','FDR_GSEA']].to_string(index=False))
print('CROSS_COHORT_FOCUS')
cf=cross[cross.apply(selected,axis=1)]
print(cf[['library','display_name','median_Z_difference_GSE175540','median_Z_difference_TJRCC','FDR_patient_sign_GSE175540','FDR_patient_sign_TJRCC','NES_GSE175540','NES_TJRCC','FDR_GSEA_GSE175540','FDR_GSEA_TJRCC','adjusted_direction_retained','macro_share_direction_retained']].head(35).to_string(index=False))
