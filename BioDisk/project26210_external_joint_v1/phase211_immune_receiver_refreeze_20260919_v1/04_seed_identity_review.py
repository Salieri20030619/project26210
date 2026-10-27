from pathlib import Path
import pandas as pd,numpy as np
O=Path.home()/'project26210_external_joint_v1/phase211_immune_receiver_refreeze_20260919_v1';d=pd.read_csv(O/'RECEIVER_CELL_METADATA_WITH_CLUSTERS.tsv.gz',sep='\t',index_col=0);a=np.load(O/'RAW_PANEL_CACHE.npz');genes=list(a['genes']);x=a['counts'];rows=[]
for cl,g in d.groupby('receiver_cluster_seed42'):
 mask=(d.receiver_cluster_seed42==cl).to_numpy();r={'cluster':cl,'cells':len(g),'primary_clusters':str(g.receiver_cluster.value_counts().head(5).to_dict())}
 for gene in ['CD3D','CD3E','TRAC','CD4','IL7R','CD8A','CD8B','FOXP3','IKZF2','IL2RA','KLRD1','NCR1','TRDC','TRGC1','TRGC2','KIT','GATA3','RORA','LST1','LYZ']:
  r[gene]=(x[mask,genes.index(gene)]>0).mean()
 rows.append(r)
f=pd.DataFrame(rows);f.to_csv(O/'SEED42_IDENTITY_MARKER_REVIEW.tsv',sep='\t',index=False);print(f.round(2).to_string(index=False))
