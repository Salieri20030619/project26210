from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy.stats import rankdata
w=Path(__file__).resolve().parents[1];q=w/'qc';root=w.parent
b=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1')
m=pd.read_csv(b/'07_pyscenic_v1/input/cell_metadata_v1.tsv.gz',sep='\t',index_col=0)
a=pd.read_csv(b/'09_metabolism_trajectory/metabolism_v1/source_data/four_pathway_AUCell_cell_scores_v1.tsv.gz',sep='\t',index_col=0).loc[m.index]
kegg=pd.read_csv(b/'09_metabolism_trajectory/metabolism_v1/source_data/KEGG_AUCell_all_pathways_v1.tsv.gz',sep='\t',index_col=0).loc[m.index]
e=pd.read_csv(root/'phase129_sections2_6_2_7_redraw_remote_v1/data/ENZYME_LOG1P_CP10K.tsv.gz',sep='\t',index_col=0).loc[m.index]
eraw=pd.read_csv(root/'phase129_sections2_6_2_7_redraw_remote_v1/data/ENZYME_RAW_COUNTS.tsv.gz',sep='\t',index_col=0).loc[m.index]
assert len(m)==11950 and m.index.is_unique and m.curated_patient_id.nunique()==26
m['TYMP_group']=a.TYMP_group
assert np.array_equal(m.TYMP_group=='TYMP_positive',eraw.TYMP>0)
factors=['F'+str(i) for i in range(1,8)];groups=['TYMP_negative','TYMP_positive']
count=m.groupby(['curated_patient_id','dominant_factor_rank7_v1','TYMP_group']).size().unstack(fill_value=0)
# Paired descriptive display within each cNMF program; no cell-level P values.
keep=count[(count>=10).all(axis=1)];keep.to_csv(q/'PAIRED_PATIENT_PROGRAM_INCLUSION.tsv',sep='\t')
rows=[]
for (patient,factor),cnt in keep.iterrows():
 for group in groups:
  ix=m.index[(m.curated_patient_id==patient)&(m.dominant_factor_rank7_v1==factor)&(m.TYMP_group==group)]
  for path in kegg.columns:
   vals=kegg.loc[ix,path]
   rows.append({'patient':patient,'factor':factor,'group':group,'pathway':path,'cells':len(ix),'mean_AUC':vals.mean(),'fraction_AUC_positive':(vals>0).mean()})
d=pd.DataFrame(rows);d.to_csv(w/'data/PATHWAY_PATIENT_GROUP_DOT_VALUES.tsv.gz',sep='\t',index=False)
summary=d.groupby(['factor','group','pathway']).agg(mean_AUC=('mean_AUC','mean'),fraction=('fraction_AUC_positive','mean'),patients=('patient','nunique')).reset_index()
summary.to_csv(w/'data/PATHWAY_DOTPLOT_SUMMARY.tsv',sep='\t',index=False)
npats=keep.groupby(level=1).size().reindex(factors,fill_value=0)
assert npats.tolist()==[6,9,0,2,3,2,2]
# Independent aggregation: each matrix column represents one patient-group mean.
err=0.
for fac in factors:
 patients=keep.index.get_level_values(0)[keep.index.get_level_values(1)==fac]
 for group in groups:
  if not len(patients):continue
  means=[];freq=[]
  for pid in patients:
   ix=m.index[(m.curated_patient_id==pid)&(m.dominant_factor_rank7_v1==fac)&(m.TYMP_group==group)]
   means.append(kegg.loc[ix].mean().values);freq.append((kegg.loc[ix]>0).mean().values)
  obs=summary[(summary.factor==fac)&(summary.group==group)].set_index('pathway').loc[kegg.columns]
  err=max(err,float(np.max(abs(np.mean(means,axis=0)-obs.mean_AUC))),float(np.max(abs(np.mean(freq,axis=0)-obs.fraction))))
assert err<1e-12
# Full circle correlation matrix, descriptive patient-level effects. No new tests.
cols=['F1','AUC_Pyrimidine','AUC_CysMet','AUC_ArgPro','AUC_Glycerophospholipid'];cors=[]
for pid,sub in a.groupby('curated_patient_id'):
 if len(sub)<30:continue
 for i,x in enumerate(cols):
  for j,y in enumerate(cols):
   xx=sub[x].values;yy=sub[y].values
   if np.ptp(xx)==0 or np.ptp(yy)==0:continue
   rr=float(np.corrcoef(rankdata(xx),rankdata(yy))[0,1])
   cors.append({'patient':pid,'x':x,'y':y,'rho':rr})
c=pd.DataFrame(cors);c.to_csv(w/'data/CORRPLOT_PATIENT_EFFECTS.tsv',sep='\t',index=False)
cs=c.groupby(['x','y']).agg(rho=('rho','median'),patients=('patient','nunique')).reset_index();cs.to_csv(w/'data/CORRPLOT_SUMMARY.tsv',sep='\t',index=False)
original=pd.read_csv(b/'09_metabolism_trajectory/metabolism_v1/tables/four_pathway_AUCell_correlation_summary_v1.tsv',sep='\t')
for r in original[original.predictor=='F1'].itertuples():
 check=cs[(cs.x=='F1')&(cs.y==r.pathway)].iloc[0]
 assert abs(check.rho-r.median_patient_rho)<1e-12 and check.patients==r.n_patients
validation={'status':'PASS','cells':11950,'patients':26,'paired_patient_counts_per_program':npats.to_dict(),'dotplot_aggregation_error':err,'dotplot_size':'patient-equal fraction of cells with AUCell >0; not P values','no_F3_eligible_pairs':True,'corrplot':'median within-patient rho; off-diagonal pathway-pathway entries descriptive only','original_F1_four_pathway_effects_unchanged':True,'no_new_inferential_tests':True}
(q/'SECTION26_VALIDATION.json').write_text(json.dumps(validation,indent=2));print(json.dumps(validation),flush=True)
