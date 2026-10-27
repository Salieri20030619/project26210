from pathlib import Path
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
d=r/'phase198_fig2_contrast_corr_revision_20260918_v5'
p=d/'scripts/draw_top20.py';s=p.read_text()
s=s.replace("fig=plt.figure(figsize=(8.7,11.0))","""row_stats=ds[ds.pathway.isin(paths)].groupby('pathway').mean_AUC.agg(['mean','std'])
ds['row_z']=(ds.mean_AUC-ds.pathway.map(row_stats['mean']))/ds.pathway.map(row_stats['std'])
ds[ds.pathway.isin(paths)].to_csv(d/'tables/TOP20_DISPLAY_VALUES.tsv',sep='\\t',index=False)
zlim=float(abs(ds.loc[ds.pathway.isin(paths),'row_z']).max())
fig=plt.figure(figsize=(8.7,11.0))""")
s=s.replace("c=[z.mean_AUC],cmap='turbo',norm=PowerNorm(gamma=.45,vmin=0,vmax=vmax)","c=[z.row_z],cmap='RdBu_r',vmin=-zlim,vmax=zlim")
s=s.replace("cb.set_ticks([0,.01,.03,.06,.10,.20,.30,.40])","cb.set_ticks([-2,-1,0,1,2])")
s=s.replace("Patient-equal mean AUC (power scale)","Within-pathway z-score of mean AUC")
s=s.replace("Shared power scale: gamma = 0.45 | Original AUC values","Colors compare groups within each pathway | Original AUC values retained")
p.write_text(s)
