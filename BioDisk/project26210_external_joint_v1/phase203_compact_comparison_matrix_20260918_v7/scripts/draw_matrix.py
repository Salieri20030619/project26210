from pathlib import Path
import json,hashlib,shutil
import pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Wedge,FancyBboxPatch
from matplotlib.colors import LinearSegmentedColormap,Normalize
r=Path.cwd();d=r/'phase203_compact_comparison_matrix_20260918_v7'
src=r/'phase202_candidate_labels_comparison_20260918_v6/tables/CANDIDATE_VS_WHOLE_SUMMARY.tsv'
s=pd.read_csv(src,sep='\t');shutil.copy2(src,d/'tables'/src.name)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none'})
fig=plt.figure(figsize=(7.8,4.05));ax=fig.add_axes([0,0,1,1]);ax.set_xlim(0,7.8);ax.set_ylim(0,4.05);ax.set_aspect('equal');ax.axis('off')
cols=['#993F80','#146A98','#397348','#B9852D','#7D8993'];names=['C4','C10','C15','Union','All macrophages'];xs=np.linspace(2.15,7.15,5)
cmap=LinearSegmentedColormap.from_list('expression',['#E5EEF3','#A5C4D3','#527D9A','#D7A66C','#AD6434']);norm=Normalize(8,12)
ax.text(.22,3.84,'SPP1 expression and source contribution',fontsize=11,fontweight='bold',va='center')
for x,name,col,n in zip(xs,names,cols,s.n_expression):
 ax.plot([x-.30,x+.30],[3.57,3.57],color=col,lw=2.4,solid_capstyle='round')
 ax.text(x,3.39,name,ha='center',va='center',fontsize=8.2,fontweight='bold',color=col)
 ax.text(x,3.17,('C4 + C10 + C15' if name=='Union' else ''),ha='center',va='center',fontsize=6.4,color='#606975')
 ax.text(x,2.96,f'n = {int(n)}',ha='center',va='center',fontsize=7,color='#536170')
for y in [2.52,1.67,.79]:
 ax.add_patch(FancyBboxPatch((.18,y-.35),7.40,.70,boxstyle='round,pad=0.01,rounding_size=0.07',facecolor='#F6F8FA',edgecolor='none',zorder=0))
ax.text(.34,2.60,'Expression',fontweight='bold',fontsize=8.3);ax.text(.34,2.37,'log2(CPM + 1)',fontsize=6.8,color='#66717C')
ax.text(.34,1.78,'SPP1-positive',fontweight='bold',fontsize=8.3);ax.text(.34,1.55,'Cells within group',fontsize=6.8,color='#66717C')
ax.text(.34,.94,'RNA contribution',fontweight='bold',fontsize=8.3);ax.text(.34,.71,'Matched 16 patients',fontsize=6.8,color='#66717C');ax.text(.34,.51,'Whole-macrophage denominator',fontsize=5.6,color='#66717C')
for x,z,col in zip(xs,s.itertuples(),cols):
 ax.add_patch(FancyBboxPatch((x-.35,2.28),.70,.48,boxstyle='round,pad=0,rounding_size=.07',facecolor=cmap(norm(z.expression)),edgecolor='white',lw=.8))
 ax.text(x,2.52,f'{z.expression:.2f}',ha='center',va='center',fontsize=9,fontweight='bold',color='white' if z.expression>9.4 else '#244455')
 rad=.29*np.sqrt(z.detection/100);ax.add_patch(Circle((x,1.73),rad,facecolor=col,edgecolor='white',lw=.6,alpha=.9))
 ax.text(x,1.30,f'{z.detection:.1f}%',ha='center',va='center',fontsize=7.5,color='#263745')
 ax.add_patch(Wedge((x,.82),.275,0,360,width=.065,facecolor='#DEE4E8',edgecolor='none'))
 if z.contribution>0:ax.add_patch(Wedge((x,.82),.275,90-3.6*z.contribution,90,width=.065,facecolor=col,edgecolor='none'))
 ax.text(x,.82,f'{z.contribution:.1f}%',ha='center',va='center',fontsize=6.4 if z.contribution>=100 else 7.6,fontweight='bold',color=col)
cax=fig.add_axes([.045,.030,.145,.017]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation='horizontal');cb.set_ticks([8,10,12]);cb.ax.tick_params(labelsize=5.6,length=1,pad=1);cb.outline.set_visible(False)
ax.text(1.63,.17,'Color: expression  |  Bubble area: positive cells  |  Ring arc: RNA share',fontsize=6.3,va='center',color='#596773')
for ext in ['pdf','svg','png']:fig.savefig(d/'figures'/f'SPP1_candidate_comparison_matrix_v7.{ext}',dpi=600 if ext=='png' else 150)
fig.savefig(d/'qc/matrix.jpg',dpi=200)
checks={'source':str(src),'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'data_values_unchanged':True,'expression_scale':[8,12],'bubble_area_proportional_to_detection':True,'ring_angle_degrees':(s.contribution*3.6).tolist(),'expression_n':s.n_expression.tolist(),'contribution_n':s.n_contribution.tolist(),'size_inches':[7.8,4.05]}
(d/'qc/VALIDATION.json').write_text(json.dumps(checks,indent=2))
print('MATRIX_READY')
