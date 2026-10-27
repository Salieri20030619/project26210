from pathlib import Path
import numpy as np,pandas as pd,json,shutil
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase199_fig3_marker_GSEA_20260918_v5'
old=d/'wilcoxon_rank_sensitivity';old.mkdir()
for folder in ['tables','panels','qc']:
 (old/folder).mkdir()
 for p in (d/folder).iterdir():
  if p.is_file():shutil.copy2(p,old/folder/p.name)
a=pd.read_csv(d/'tables/RANKED_PATIENT_GENES.tsv',sep='\t')
z=np.load(r/'part2_rebuild_v1/04_markers_annotation/run_v4_macrophage_only/paired_patient_effects/cluster_4.npz',allow_pickle=True)
idx={str(g):i for i,g in enumerate(z['genes'])}
v=z['effects'][:,[idx[g] for g in a.gene]]
assert v.shape==(26,len(a)) and np.isfinite(v).all()
mean=v.mean(0);sd=v.std(0,ddof=1)
assert (sd>0).all()
a['ranking_score']=mean/(sd/np.sqrt(26))
a['paired_mean_delta']=mean
a['paired_sd_delta']=sd
a=a.sort_values(['ranking_score','gene'],ascending=[False,True],kind='stable')
a.to_csv(d/'tables/RANKED_PATIENT_GENES.tsv',sep='\t',index=False)
(d/'qc/CONTINUOUS_RANK_AUDIT.json').write_text(json.dumps({'patients':26,'genes':len(a),'ties':int(a.ranking_score.duplicated().sum()),'ranking':'paired effect mean / (sample SD / sqrt(26))','same_eligible_genes':True,'new_gene_P_tests':False},indent=2))
print('CONTINUOUS_RANK_READY',len(a),'ties',a.ranking_score.duplicated().sum())
