from common import *
def main():
 src=read(P/'tables/SOURCE_GENE_PATIENT_DESCRIPTIVE.tsv')
 src=src[src.tissue.eq('Tumor')&src.eligible_descriptive&src.source_group.isin(['Conservative_malignant','Frozen_C4','Frozen_other_mature_macrophage','Refined_stromal_pool'])]
 out=src.groupby(['source_group','gene']).agg(patients=('patient','nunique'),median_detection=('detection_fraction','median'),median_logCP10K=('logCP10K','median')).reset_index()
 tab(out,'SOURCE_COMPARTMENT_DESCRIPTIVE_SUMMARY.tsv')
 g=read(R/'phase157_ordered_bridge_reassessment_20260915_v1/02_inputs/GROUP_INDEX.tsv');g=g[g.source_tissue.eq('Tumor')]
 x=g.pivot_table(index=['study','patient','sample_id'],columns='role',values='cells',aggfunc='sum',fill_value=0).reset_index();x['macro_cells']=x.get('C4',0)+x.get('OtherMacrophage',0)
 s=read(P/'tables/STROMAL_STATE_PSEUDOBULK.tsv');s=s[s.tissue.eq('Tumor')&s.population.eq('stromal_pool')][['study','patient','library','cells']].rename(columns={'library':'sample_id','cells':'stromal_cells'})
 z=x.merge(s,on=['study','patient','sample_id'],how='left');z['eligible_all_current']=(z.Malignant>=50)&(z.macro_cells>=100)&(z.stromal_cells>=50)
 old=read(P/'tables/MATCHED_SOURCE_RECEIVER_INPUTS.tsv');z['in_legacy16']=z.patient.isin(old.patient)
 assert set(z.loc[z.eligible_all_current,'patient'])==set(old.loc[old.eligible,'patient'])
 tab(z,'ALL_LIBRARY_SOURCE_CAPACITY_AUDIT.tsv')
 rec=[]
 for g in ['CD44','ITGAV','ITGB1','TGFBR1','TGFBR2']:
  d=old[old.eligible];rec.append(dict(gene=g,patients=len(d),patients_positive=int(d[g+'_positive_cells'].gt(0).sum()),minimum_positive_cells=int(d[g+'_positive_cells'].min()),interpretation='RNA_detection_not_protein_activity'))
 tab(rec,'RECEPTOR_DETECTION_QUALIFICATION.tsv')
 print('REVIEW_TABLES_COMPLETE',flush=True)
if __name__=='__main__':main()
