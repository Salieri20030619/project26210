from pathlib import Path
import json,hashlib,shutil
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
from matplotlib.colors import to_rgb
r=Path.cwd();d=r/'phase204_concentric_comparison_20260918_v8'
src=r/'phase202_candidate_labels_comparison_20260918_v6/tables/CANDIDATE_VS_WHOLE_SUMMARY.tsv'
s=pd.read_csv(src,sep='\t');shutil.copy2(src,d/'tables'/src.name)
plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'svg.fonttype':'none','font.size':8})
fig=plt.figure(figsize=(8.2,3.1));ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,8.2),ylim=(0,3.1),aspect='equal');ax.axis('off')
cols=['#993F80','#146A98','#397348','#B9852D','#7D8993'];names=['C4','C10','C15','Union','All macrophages'];xs=np.linspace(.86,7.34,5)
ax.text(.26,2.91,'SPP1 expression and source contribution',fontsize=11,fontweight='bold',va='center')
for x,z,col,name in zip(xs,s.itertuples(),cols,names):
 ax.text(x,2.60,name,ha='center',va='center',fontsize=9,fontweight='bold',color=col)
 ax.text(x,2.36,f'n = {int(z.n_expression)}'+('  |  C4+C10+C15' if name=='Union' else ''),ha='center',va='center',fontsize=6.3,color='#626F7B')
 light=tuple(.52*np.array(to_rgb(col))+.48)
 for radius,width,value,color in [(.58,.105,z.contribution,col),(.435,.095,z.detection,light)]:
  ax.add_patch(Wedge((x,1.61),radius,0,360,width=width,facecolor='#E8EDF0',edgecolor='none'))
  if value>0:ax.add_patch(Wedge((x,1.61),radius,90-3.6*value,90,width=width,facecolor=color,edgecolor='none'))
 ax.text(x,1.65,f'{z.expression:.2f}',ha='center',va='center',fontsize=13,fontweight='bold',color='#283B49')
 ax.text(x,1.40,'SPP1',ha='center',va='center',fontsize=6.2,color='#7A8690')
 ax.text(x,.80,f'{z.detection:.1f}%  /  {z.contribution:.1f}%',ha='center',va='center',fontsize=8,color=col)
ax.plot([.28,7.92],[.58,.58],color='#E5E9ED',lw=.7)
ax.text(4.1,.41,'Center: expression, log2(CPM + 1)   |   Inner ring: SPP1-positive cells   |   Outer ring: SPP1 RNA share',ha='center',va='center',fontsize=6.5,color='#53616D')
ax.text(4.1,.18,'Percentages below: inner / outer.  Patient medians; n above applies to expression and positivity.  RNA share: matched n = 16.',ha='center',va='center',fontsize=6.1,color='#6D7882')
for ext in ['pdf','svg','png']:fig.savefig(d/'figures'/f'SPP1_concentric_comparison_v8.{ext}',dpi=600 if ext=='png' else 150)
(d/'qc/VALIDATION.json').write_text(json.dumps({'source':str(src),'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'values_unchanged':True,'outer_ring':'SPP1 RNA share','inner_ring':'SPP1-positive fraction','center':'SPP1 log2(CPM+1)','ring_angles':{'inner':(3.6*s.detection).tolist(),'outer':(3.6*s.contribution).tolist()},'expression_n':s.n_expression.tolist(),'rna_n':s.n_contribution.tolist(),'size_inches':[8.2,3.1]},indent=2))
print('DRAW_READY')
