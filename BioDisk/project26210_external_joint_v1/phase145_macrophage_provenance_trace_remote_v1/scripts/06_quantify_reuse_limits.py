from pathlib import Path
import pandas as pd,numpy as np,json
O=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase145_macrophage_provenance_trace_remote_v1');a=pd.read_csv(O/'legacy_object/metadata.tsv.gz',sep='\t',low_memory=False);j=pd.read_csv(O/'identity/cell_label_crosswalk.tsv.gz',sep='\t')
res=[]
for sub,g in a.groupby('macro_subtype'):
 counts=g.study.value_counts();jj=j[j.macro_subtype==sub];m=jj[jj.cell_type_broad_v1=='Macrophage']
 res.append(dict(subtype=sub,n_old=len(g),dominant_study=counts.index[0],dominant_study_n=int(counts.iloc[0]),dominant_study_fraction=float(counts.iloc[0]/len(g)),old_libraries=g.sample_id.nunique(),n_current_retained=len(jj),n_current_macrophage=len(m),n_current_other_lineage=int((jj.cell_type_broad_v1!='Macrophage').sum()),n_missing_from_current=len(g)-len(jj),n_current_spp1_state=int((m.state_v2=='SPP1-state Macrophage').sum()),current_patients=jj.curated_patient_id.nunique()))
r=pd.DataFrame(res);r.to_csv(O/'identity/LEGACY_REUSE_SUMMARY.tsv',sep='\t',index=False);print(r.to_string(index=False))
qcpath=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/01_raw_counts_qc/qc_v1/tables/cell_qc_metrics_v1.tsv.gz')
qc=pd.read_csv(qcpath,sep='\t',low_memory=False);print('QC_SCHEMA',list(qc.columns));print('QC_ID',qc.iloc[0,:3].to_dict())
idcol=next((c for c in ['exact_cell_id','cell_id'] if c in qc),qc.columns[0]);missing=set(a.cell_id)-set(j.cell_id);m=qc[qc[idcol].isin(missing)].copy();m.to_csv(O/'identity/OLD_CELLS_ABSENT_CURRENT_QC.tsv.gz',sep='\t',index=False,compression='gzip')
print('MISSING_QC_MATCH',len(m),'of',len(missing))
for c in m:
 if any(t in c.lower() for t in ['reason','keep','doublet','quality','filter','exclude']):print('QC',c,m[c].value_counts(dropna=False).to_dict())
print('MATCHED_GENE_COUNTS',j.groupby('macro_subtype').size().to_dict())
