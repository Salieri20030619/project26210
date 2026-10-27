from pathlib import Path
import pandas as pd,numpy as np,json
O=Path.home()/'project26210_external_joint_v1/phase211_immune_receiver_refreeze_20260919_v1';d=pd.read_csv(O/'RECEIVER_CELL_METADATA_WITH_CLUSTERS.tsv.gz',sep='\t',index_col=0);a=np.load(O/'RAW_PANEL_CACHE.npz');genes=list(a['genes']);x=a['counts'];rows=[]
for cl in [3,16]:
 for seed,g in d[d.receiver_cluster==cl].groupby('receiver_cluster_seed42'):
  mask=((d.receiver_cluster==cl)&(d.receiver_cluster_seed42==seed)).to_numpy();r={'unit':f'{cl}.{seed}','cells':len(g),'patients':g.curated_patient_id.nunique(),'old_candidates':json.dumps(g.legacy_candidate.value_counts().to_dict())}
  for gene in ['CD3D','CD3E','TRAC','CD4','IL7R','CD8A','CD8B','FOXP3','KLRD1','NCR1','TRDC','TRGC1','TRGC2']:r[gene]=(x[mask,genes.index(gene)]>0).mean()
  rows.append(r)
f=pd.DataFrame(rows);f.to_csv(O/'MIXED_PARENT_SUBUNIT_EVIDENCE.tsv',sep='\t',index=False);print(f[f.cells>=50].round(3).to_string(index=False))
