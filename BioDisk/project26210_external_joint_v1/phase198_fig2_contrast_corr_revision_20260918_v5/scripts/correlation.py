from pathlib import Path
import json
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase198_fig2_contrast_corr_revision_20260918_v5'
o=r/'phase131_sections2_5_2_6_reference_correction_remote_v1'
cs=pd.read_csv(o/'data/CORRPLOT_SUMMARY.tsv',sep='\t')
eff=pd.read_csv(o/'data/CORRPLOT_PATIENT_EFFECTS.tsv',sep='\t')
cols=['F1','AUC_Pyrimidine','AUC_CysMet','AUC_ArgPro','AUC_Glycerophospholipid']
names=['F1','Pyrimidine','Cys / Met','Arg / Pro','GPL']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none'})
mat=cs.pivot(index='x',columns='y',values='rho').loc[cols,cols].values
fig=plt.figure(figsize=(7.2,6.2))
ax=fig.add_axes([.19,.29,.65,.66])
rec=[]
for i in range(5):
 for j in range(5):
  if i>j:
   rho=mat[i,j]
   ax.scatter(j,i,s=abs(rho)*2100,c=[rho],cmap='RdBu_r',vmin=-.5,vmax=.5,edgecolors='none',zorder=3)
   ax.text(j,i,f'{rho:.2f}',ha='center',va='center',fontsize=10,color='white' if abs(rho)>.34 else '#30353A',zorder=4)
  elif i<j:
   v=eff[(eff.x==cols[i])&(eff.y==cols[j])].rho.dropna()
   assert len(v)==20
   assert abs(float(v.median())-mat[i,j])<1e-12
   k=int((v>0).sum());fraction=k/len(v)
   ax.add_patch(Rectangle((j-.34,i+.08),.68,.12,fc='#E5E7E8',ec='none'))
   ax.add_patch(Rectangle((j-.34,i+.08),.68*fraction,.12,fc='#408E91',ec='none'))
   ax.text(j,i-.11,f'{k}/{len(v)}',ha='center',va='center',fontsize=9,color='#285B5D')
   rec.append({'x':cols[i],'y':cols[j],'median_rho':float(mat[i,j]),'positive':k,'n':len(v)})
  else:
   ax.text(j,i,'-',ha='center',va='center',color='#AAAAAA')
ax.set(xlim=(-.5,4.5),ylim=(4.5,-.5),aspect='equal')
ax.set_xticks(range(5),names,fontsize=9,rotation=35,ha='right')
ax.set_yticks(range(5),names,fontsize=9)
ax.set_xticks(np.arange(-.5,5,1),minor=True);ax.set_yticks(np.arange(-.5,5,1),minor=True)
ax.grid(which='minor',color='#E6E8EA',lw=.6);ax.tick_params(which='both',length=0)
for sp in ax.spines.values():sp.set_visible(False)
ca=fig.add_axes([.87,.48,.022,.36])
cb=fig.colorbar(ScalarMappable(norm=Normalize(-.5,.5),cmap='RdBu_r'),cax=ca,ticks=[-.5,-.25,0,.25,.5])
cb.set_label('Median patient Spearman rho',fontsize=8);cb.ax.tick_params(labelsize=8)
fig.text(.20,.982,'Lower: median rho',fontsize=9,weight='bold')
fig.text(.57,.982,'Upper: positive patients',fontsize=9,color='#285B5D',weight='bold')
summ=pd.read_csv('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/09_metabolism_trajectory/metabolism_v1/tables/four_pathway_AUCell_correlation_summary_v1.tsv',sep='\t')
for k,(col,label) in enumerate(zip(cols[1:],names[1:])):
 z=summ[(summ.pathway==col)&(summ.predictor=='F1')].iloc[0]
 row=next(t for t in rec if t['x']=='F1' and t['y']==col)
 assert row['positive']==round(z.fraction_positive*z.n_patients)
 fig.text(.06+(k%2)*.49,.15-(k//2)*.068,f'F1 - {label}: P={z.signed_rank_PValue:.3g}; FDR={z.FDR:.3g}',fontsize=8)
fig.text(.06,.016,'20 patients | Circle area = |rho| | Bar = fraction with rho > 0\nGPL: glycerophospholipid | Original four-pathway FDR family retained',fontsize=7,linespacing=1.5)
fig.savefig(d/'panels/correlation_consistency.pdf',bbox_inches='tight',pad_inches=.035)
fig.savefig(d/'figures/Standalone_F1_pathway_correlation_consistency_v5.pdf',bbox_inches='tight',pad_inches=.035)
plt.close(fig)
pd.DataFrame(rec).to_csv(d/'tables/CORRELATION_AND_CONSISTENCY.tsv',sep='\t',index=False)
print('CORRELATION_READY')
