from pathlib import Path
import numpy as np,pandas as pd,json,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
w=Path(__file__).resolve().parents[1]
genes=json.loads((w/'qc/VALIDATION.json').read_text())['gene_panel']
clmedian=pd.read_csv(w/'qc/K1_GROUP_MEDIANS.tsv',sep='\t',index_col=0)
cnmedian=pd.read_csv(w/'qc/K2_GROUP_MEDIANS.tsv',sep='\t',index_col=0)
maxrad=float(np.ceil(max(clmedian.to_numpy().max(),cnmedian.to_numpy().max())*2)/2);angles=np.linspace(0,2*np.pi,len(genes),endpoint=False);closed=np.r_[angles,angles[0]];colors=['#C95D6A','#4479A9']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'axes.linewidth':.7})
def draw(df,title,n,name,subtitle):
 fig=plt.figure(figsize=(7.1,6.7));ax=fig.add_axes([.14,.22,.66,.60],projection='polar');ax.set_theta_offset(np.pi/2);ax.set_theta_direction(-1);ax.set_xticks(angles,genes,fontsize=10,fontstyle='italic');ax.tick_params(axis='x',pad=9);ax.set_ylim(0,maxrad);ticks=np.linspace(0,maxrad,5)[1:];ax.set_yticks(ticks,[f'{t:g}' for t in ticks],fontsize=8,color='#777777');ax.set_rlabel_position(18);ax.grid(color='#C8CCD1',lw=.6);ax.spines['polar'].set_color('#9AA0A6');ax.spines['polar'].set_linewidth(.7)
 for (group,row),color in zip(df.iterrows(),colors):
  vals=row[genes].to_numpy();ax.plot(closed,np.r_[vals,vals[0]],color=color,lw=1.9,marker='o',markersize=4.2,label=group);ax.fill(closed,np.r_[vals,vals[0]],color=color,alpha=.10)
 fig.text(.025,.96,'K',fontsize=19,fontweight='bold',va='top');fig.text(.14,.956,title,fontsize=12,fontweight='bold',va='top');fig.text(.14,.912,f'n = {n} paired patients'+subtitle,fontsize=9,color='#555555');fig.legend(*ax.get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.49,.086),ncol=2,frameon=False,fontsize=10);fig.text(.14,.035,'Radius: median of patient-level mean log1p RNA expression\nSame ten genes and radial scale in both versions.',fontsize=8.5,color='#555555')
 for ext in ['png','pdf']:fig.savefig(w/'figures'/(name+'.'+ext),dpi=350,bbox_inches='tight',pad_inches=.12)
 plt.close(fig)
draw(clmedian,'F1-high clone vs other clones',14,'K1_F1_high_vs_other_clones_radar_v2','')
draw(cnmedian,'Malignant vs non-malignant epithelial cells',3,'K2_malignant_vs_nonmalignant_radar_v2','; exploratory')
