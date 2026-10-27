from pathlib import Path
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=Path(__file__).resolve().parents[2];w=Path(__file__).resolve().parents[1];d=Path('C:/Users/31937/Downloads/26210/02_2.1_G_I重绘_v2');d.mkdir(parents=True,exist_ok=True)
a=pd.read_csv(r/'phase32_part1_m03_cnv_clone_trees_local_v1/input_snapshot_local_v1/formal_cell_cnv_clone_assignments_v1.tsv.gz',sep='\t');b=pd.read_csv(r/'phase43_part1_client_figures_local_v1/source_remote_v1/tumor_epithelial_malignancy_calls_v1.tsv.gz',sep='\t');a=a.merge(b[['exact_cell_id','infercnv_cnv_gene_burden','infercnv_high_cnv_threshold']],on='exact_cell_id',validate='one_to_one');assert len(a)==11950
a['value']=np.log2((a.infercnv_cnv_gene_burden+1)/(a.infercnv_high_cnv_threshold+1));pts=['GSE207493__RCC104','GSE207493__RCC86','GSE156632__RCC4','GSE156632__RCC5'];colors=['#8BBBD0','#E3B546','#73A56B','#A76B9F','#BF8F7C','#A5B6BD','#8C9D63'];rng=np.random.default_rng(26210)
plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.titlesize':11,'axes.linewidth':.8,'pdf.fonttype':42})
for kind in ['violin_points','box_points']:
 fig,axes=plt.subplots(2,2,figsize=(10,8))
 for ax,pt in zip(axes.ravel(),pts):
  pc=a[(a.infercnv_unit==pt)&(a.clone_cells>=20)];clones=pc.groupby('formal_inferred_cnv_clone_v1').value.median().sort_values(kind='stable').index.tolist();counts=[]
  for j,cl in enumerate(clones):
   cc=pc[pc.formal_inferred_cnv_clone_v1==cl];v=cc.value.to_numpy();co=colors[int(cl.rsplit('C',1)[-1])-1];counts.append(len(v))
   if kind=='violin_points' and np.ptp(v)>1e-10:
    viol=ax.violinplot([v],positions=[j],widths=.8,showextrema=False,bw_method=.3);bb=viol['bodies'][0];bb.set_facecolor(co);bb.set_edgecolor(co);bb.set_alpha(.40);bb.set_linewidth(1)
   bp=ax.boxplot([v],positions=[j],widths=.18 if kind=='violin_points' else .50,showfliers=False,patch_artist=True,whis=1.5,medianprops={'color':'#444444','linewidth':1.3},boxprops={'facecolor':co,'edgecolor':co,'alpha':.55},whiskerprops={'color':co,'linewidth':1},capprops={'color':co,'linewidth':1})
   # Points represent source HMM groups, not duplicated cells sharing one burden.
   dots=cc.groupby('infercnv_group',sort=True).value.first().to_numpy();offset=np.linspace(-.16,.16,len(dots)) if len(dots)>1 else [0]
   ax.scatter(j+np.asarray(offset),dots,s=20,c=co,edgecolors=co,linewidths=.7,alpha=.8,zorder=4)
  ax.set_xticks(range(len(clones)),[cl.rsplit('_',1)[-1] for cl in clones]);ax.set_xlabel('CNV subclone');ax.set_ylabel('Relative CNV burden');ax.set_title(pt.replace('__','  '),fontweight='normal');ax.spines[['top','right']].set_visible(False);ax.tick_params(length=3);ax.grid(False);ax.margins(x=.12,y=.15)
 fig.subplots_adjust(left=.085,right=.98,top=.92,bottom=.10,hspace=.44,wspace=.3);fig.text(.015,.97,'G',fontsize=19,fontweight='bold');fig.text(.085,.025,'Points: original HMM groups. Clones ordered by median CNV burden within each patient.',fontsize=9,color='#555555')
 for ext in ['png','pdf']:fig.savefig(d/('G_CNV_'+kind+'_v2.'+ext),dpi=300,bbox_inches='tight',pad_inches=.12)
 plt.close(fig)
print('G variants complete')
