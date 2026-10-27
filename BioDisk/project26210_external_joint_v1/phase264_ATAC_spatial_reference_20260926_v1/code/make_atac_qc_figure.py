from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42})
fig,axes=plt.subplots(2,2,figsize=(8.27,9.7));fig.subplots_adjust(left=.12,right=.94,bottom=.23,top=.9,hspace=.38,wspace=.4)
fig.suptitle('ATAC chr1 pilot: projection coverage and score components',fontsize=12,y=.96)
for ax,name in zip(axes[0],['GLMF1','GLMF2']):
 d=pd.read_csv(P/f'plot_data/{name}_spatial_ATAC_projection.tsv',sep='\t')
 sc=ax.scatter(d.x_hires_px,d.y_hires_px,c=d.represented_fraction,s=3,cmap='viridis',vmin=0,vmax=1,rasterized=True)
 ax.invert_yaxis();ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([])
 ax.set_title(f'{name}: {d.projection_eligible.sum()}/{len(d)} eligible')
 fig.colorbar(sc,ax=ax,shrink=.8,label='Represented abundance fraction')
a=pd.read_csv(P/'tables/ATAC_CHR1_TF_ACTIVITY.tsv',sep='\t',index_col=0)
s=pd.read_csv(P/'differential/RCC100_chr1_statistics.txt',sep='\t').set_index('Motif')
ps=pd.DataFrame({tf:s.loc[[i for i in s.index if i.endswith('.'+tf)][0],['Protection_Score_'+l for l in a.index]].to_numpy() for tf in a.columns},index=a.index).astype(float)
for ax,m,title,cm in [(axes[1,0],a,'HINT activity = PS + TC','viridis'),(axes[1,1],ps,'Protection score (PS)','RdBu_r')]:
 limit=float(np.abs(m.to_numpy()).max()); kw=dict(vmin=-limit,vmax=limit) if cm=='RdBu_r' else {}
 im=ax.imshow(m,cmap=cm,aspect='auto',**kw)
 ax.set_xticks(range(4),m.columns,rotation=35,ha='right');ax.set_yticks(range(4),m.index,fontsize=8);ax.set_title(title,fontsize=10)
 for i in range(4):
  for j in range(4):ax.text(j,i,f'{m.iloc[i,j]:.3f}',ha='center',va='center',fontsize=7,color='black',bbox=dict(facecolor='white',alpha=.7,edgecolor='none',pad=.6))
 fig.colorbar(im,ax=ax,shrink=.8)
fig.text(.1,.055,'Projection gate: represented fraction >= 0.80; conditional weights among four lineages.\nHNF4A and STAT3 have negative protection scores in all four lineages.\nPS + TC values do not establish protected binding for every displayed TF.\nOne ATAC donor, chromosome 1 and four fixed motifs; no inferential P values.',fontsize=8,linespacing=1.6)
fig.savefig(P/'delivery/26210_ATAC_chr1_coverage_and_scores_v1.pdf',dpi=300)
fig.savefig(P/'delivery/26210_ATAC_chr1_coverage_and_scores_v1.png',dpi=200)
print('QC_FIGURE_READY')
