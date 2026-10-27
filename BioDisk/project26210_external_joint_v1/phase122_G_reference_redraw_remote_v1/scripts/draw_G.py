from pathlib import Path
import json,hashlib,sys
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
w=Path(__file__).resolve().parents[1];out=w/'figures';out.mkdir(exist_ok=True)
src=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/04_subclones_tymp/cnv_subclones_v1/formal_cell_cnv_clone_assignments_v1.tsv.gz')
a=pd.read_csv(src,sep='\t');assert len(a)==11950
pts=['GSE207493__RCC104','GSE207493__RCC86','GSE156632__RCC4','GSE156632__RCC5']
frames=[]
for pt in pts:
 c=a[(a.infercnv_unit==pt)&(a.clone_cells>=20)].copy()
 scores=pd.read_csv(w/'data'/f'{pt}.tsv',sep='\t')
 mapping=pd.read_csv(Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/infercnv_v2/patients')/pt/'cell_id_mapping_v2.tsv',sep='\t')
 scores=scores.rename(columns={'exact_cell_id':'tool_id'}).merge(mapping[['tool_id','exact_cell_id']],on='tool_id',validate='one_to_one')
 c=c.merge(scores[['exact_cell_id','cnv_score']],on='exact_cell_id',validate='one_to_one',how='left')
 assert c.cnv_score.notna().all(),(pt,c.exact_cell_id[c.cnv_score.isna()].head().tolist(),scores.exact_cell_id.head().tolist())
 c['clone']=c.formal_inferred_cnv_clone_v1.str.rsplit('_',n=1).str[-1]
 frames.append(c)
allc=pd.concat(frames,ignore_index=True);allc.to_csv(w/'data'/'plot_cells_server_only.tsv.gz',sep='\t',index=False)
colors={'C01':'#77A9C7','C02':'#E8AD65','C03':'#82AD7E','C04':'#CB797C','C05':'#A28DBB','C06':'#A99A90','C07':'#CF98BD'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.linewidth':.7,'pdf.fonttype':42,'svg.fonttype':'none','axes.labelcolor':'#333333','xtick.color':'#333333','ytick.color':'#333333'})
upper=float(allc.cnv_score.max())*1.055
summ=[]
def panel(ax,c,pt):
 clones=sorted(c.clone.unique());rng=np.random.default_rng(26210)
 for j,cl in enumerate(clones):
  v=c.loc[c.clone==cl,'cnv_score'].to_numpy();co=colors[cl]
  if np.ptp(v)>1e-12:
   vp=ax.violinplot(v,positions=[j],widths=.76,points=180,bw_method=.28,showextrema=False)
   body=vp['bodies'][0];body.set_facecolor(co);body.set_alpha(.90);body.set_edgecolor('#303030');body.set_linewidth(.7)
  # Show actual observations on one side of each violin; no invented y jitter.
  ids=np.linspace(0,len(v)-1,min(60,len(v)),dtype=int);dots=np.sort(v)[ids]
  ax.scatter(j+rng.uniform(-.14,.14,len(dots)),dots,s=2.7,c='#222222',alpha=.25,linewidths=0,rasterized=True,zorder=3)
  q=np.quantile(v,[.25,.5,.75]);ax.plot([j,j],[q[0],q[2]],c='#222222',lw=1.6,zorder=4);ax.scatter(j,q[1],s=8,c='white',edgecolors='#222222',linewidths=.5,zorder=5)
 ax.set_xlim(-.65,len(clones)-.35);ax.set_ylim(0,upper);ax.set_xticks(range(len(clones)),clones,fontsize=7.5)
 ax.tick_params(axis='both',length=2.5,width=.7,pad=3)
 for spine in ax.spines.values():spine.set_color('#666666')
 ax.text(.965,.95,pt.split('__')[1],transform=ax.transAxes,ha='right',va='top',fontsize=10,fontweight='bold',color='#8F4037')
 ax.text(.035,.95,pt.split('__')[0],transform=ax.transAxes,ha='left',va='top',fontsize=7,color='#666666')
 ax.set_xlabel('CNV subclone',labelpad=5);ax.set_ylabel('CNV score',labelpad=5)
def save(fig,name):
 for ext in ['png','pdf']:fig.savefig(out/(name+'.'+ext),dpi=350,facecolor='white',bbox_inches='tight',pad_inches=.07)
 plt.close(fig)
fig,axs=plt.subplots(2,2,figsize=(7.6,5.6),sharey=True)
for ax,c,pt in zip(axs.ravel(),frames,pts):panel(ax,c,pt)
fig.subplots_adjust(left=.09,right=.99,bottom=.095,top=.925,wspace=.21,hspace=.36)
fig.text(.015,.97,'G',fontsize=18,fontweight='bold',va='top')
save(fig,'G_CNV_reference_complete_v3')
for c,pt in zip(frames,pts):
 fig,ax=plt.subplots(figsize=(3.7,2.65));panel(ax,c,pt);fig.subplots_adjust(left=.17,right=.98,bottom=.19,top=.97);save(fig,'G_'+pt+'_v3')
 for cl,g in c.groupby('clone'):
  summ.append({'patient':pt,'clone':cl,'cells':len(g),'unique_continuous_scores':g.cnv_score.nunique(),'min':g.cnv_score.min(),'median':g.cnv_score.median(),'max':g.cnv_score.max()})
pd.DataFrame(summ).to_csv(w/'qc'/'G_aggregate_summary.tsv',sep='\t',index=False)
val={'status':'PASS_DATA_CHECKS','cells':len(allc),'patients':4,'all_source_assignments':len(a),'clone_labels_unchanged':True,'matched_every_selected_cell':True,'finite_scores':bool(np.isfinite(allc.cnv_score).all()),'score':'mean(abs(log2(preliminary_infercnv_relative_signal))) across retained genes per cell','old_score':'log2((HMM altered-gene burden+1)/(patient threshold+1)); not the same metric','points':'up to 60 evenly spaced order statistics per clone; actual y values, horizontal display jitter only','violin':'all retained cells, KDE bw .28 clipped at observed min/max','inference':'none; patients are biological units; cells are descriptive observations','versions':{'python':sys.version,'matplotlib':matplotlib.__version__,'numpy':np.__version__,'pandas':pd.__version__}}
(w/'qc'/'G_VALIDATION.json').write_text(json.dumps(val,indent=2))
(out/'G_README.txt').write_text('G 图重绘 v3\n\n按参考图重排为紧凑四分面；统一纵轴、细黑轮廓、柔和实色、四分位线与中位数白点。保留原有四位患者和正式亚克隆标签，克隆按编号排列。\n\n小提琴使用服务器原始 preliminary inferCNV 连续相对信号：每个细胞的 CNV score = 各保留基因 |log2(relative signal)| 的均值。与上一版基于 HMM 分组的 relative burden 不是同一个指标，数值不能直接比较；原版及阈值仍保留，未重新判定恶性或亚克隆。\n\n小提琴涵盖全部合格细胞，叠加点最多每克隆60个实际观测分位位置，仅横向错开。细胞是描述性观测，不是独立生物学重复；本图不新增 P 值/FDR，也不代表恶性程度或因果性。\n\n主图：G_CNV_reference_complete_v3。另附四张患者单图便于排版。\n',encoding='utf8')
print(json.dumps(val))
