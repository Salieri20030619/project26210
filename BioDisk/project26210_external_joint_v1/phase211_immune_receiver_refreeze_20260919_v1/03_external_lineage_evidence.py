from pathlib import Path
import pandas as pd,numpy as np
O=Path.home()/'project26210_external_joint_v1/phase211_immune_receiver_refreeze_20260919_v1';rows=[]
for ds in ['GSE264586','GSE294109']:
 d=pd.read_csv(O/f'external/{ds}_CELL_MAPPING.tsv.gz',sep='\t',index_col='cell_id');c=pd.read_csv(O/f'external/{ds}_TNK_PANEL_RAW_COUNTS.tsv.gz',sep='\t',index_col='cell_id');d=d.loc[c.index]
 for cl,g in d.groupby('rds_phase7_cluster_res0_8'):
  x=c.loc[g.index];T=(x[['CD3D','CD3E','TRAC','CD3G','TRBC1','TRBC2']]>0).sum(axis=1)>=2;gdt=(x[['TRDC','TRGC1','TRGC2']]>0).sum(axis=1)>=1
  row={'dataset':ds,'cluster':cl,'old_label':g.phase7_fine_annotation.iloc[0],'cells':len(g),'patients':g.patient_key.nunique(),'T_backbone_2plus':T.mean(),'T_and_CD8B':(T&(x.CD8B>0)).mean(),'T_and_CD4':(T&(x.CD4>0)).mean(),'T_and_FOXP3':(T&(x.FOXP3>0)).mean(),'T_and_gdT':(T&gdt).mean(),'NK_specific_no_T':(((x[['KLRD1','NCR1','KLRF1','NCAM1','FCGR3A']]>0).sum(axis=1)>=2)&~T).mean()}
  for gene in ['CD3D','CD3E','TRAC','TRBC1','TRBC2','CD4','IL7R','CD8A','CD8B','FOXP3','IKZF2','IL2RA','CTLA4','TRDC','TRGC1','TRGC2','NCR1','KLRD1','KLRF1','FCGR3A']:row[gene+'_pct']=(x[gene]>0).mean()
  rows.append(row)
a=pd.DataFrame(rows);a.to_csv(O/'external/EXTERNAL_RAW_CLUSTER_LINEAGE_EVIDENCE.tsv',sep='\t',index=False);print(a[['dataset','cluster','old_label','cells','T_backbone_2plus','T_and_CD8B','T_and_CD4','T_and_FOXP3','T_and_gdT','NK_specific_no_T','NCR1_pct','KLRD1_pct']].round(3).to_string(index=False))
