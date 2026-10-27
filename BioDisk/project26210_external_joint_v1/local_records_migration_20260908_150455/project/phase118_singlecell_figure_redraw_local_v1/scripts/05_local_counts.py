from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import pandas as pd,hashlib,json
r=Path(__file__).resolve().parents[1];src=r.parent/'phase43_part1_client_figures_local_v1/source_remote_v1/celltype_counts_v1.tsv'
d=pd.read_csv(src,sep='\t').sort_values('n_cells',ascending=False);assert d.n_cells.sum()==258120
palette={'T cell':'#81729D','NK cell':'#A5A76B','B cell':'#8D719D','Plasma cell':'#BC8CA6','Macrophage':'#D76A73','Monocyte':'#E1A078','Dendritic cell':'#7AAFA0','Mast cell':'#AF7694','Epithelial cell':'#D4AA62','Endothelial':'#6AAE99','Pericyte/SMC':'#7DA2BF','Fibroblast':'#A3BC78','Erythrocyte':'#C57B65','Epithelial/Myeloid mixed':'#A7A7A7'}
plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.linewidth':.7,'pdf.fonttype':42})
cs=d.n_cells.to_numpy();ns=d.cell_type_broad_v1.tolist();co=[palette[n] for n in ns]
for h in [True,False]:
 fig,ax=plt.subplots(figsize=(7.4,6) if h else (10,6));ix=range(len(d))
 if h:
  ax.barh(ix,cs,color=co,height=.65);ax.set_yticks(ix,ns);ax.invert_yaxis();ax.set_xlim(0,80000);ax.set_xlabel('Number of cells');ax.xaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{v/1000:.0f}k'))
  for i,v in enumerate(cs):ax.text(v+900,i,f'{v:,}  ({v/cs.sum()*100:.2f}%)',va='center',fontsize=8)
 else:
  ax.bar(ix,cs,color=co,width=.66);ax.set_xticks(ix,ns,rotation=48,ha='right');ax.set_ylim(0,74000);ax.set_ylabel('Number of cells');ax.yaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{v/1000:.0f}k'))
  for i,v in enumerate(cs):ax.text(i,v+850,f'{v:,}',ha='center',va='bottom',fontsize=8,rotation=90)
 ax.spines[['top','right']].set_visible(False);ax.tick_params(length=3);ax.set_title('Cell composition',loc='left',fontweight='bold',pad=18);ax.text(1,1.02,'258,120 cells',transform=ax.transAxes,ha='right',fontsize=9,color='#666666');fig.tight_layout()
 stem='01_counts_horizontal_v1' if h else '02_counts_vertical_v1'
 for ext in ['png','pdf']:fig.savefig(r/'figures'/f'{stem}.{ext}',dpi=300,bbox_inches='tight',pad_inches=.12)
 plt.close(fig)
meta={'status':'LOCAL_COUNTS_VALIDATED','n_cells':int(cs.sum()),'n_groups':len(d),'source':str(src),'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'versions':{'matplotlib':matplotlib.__version__,'pandas':pd.__version__},'remote_transfer':'BLOCKED_PENDING_EXPLICIT_DESTINATION_APPROVAL'}
(r/'qc'/'local_counts_validation.json').write_text(json.dumps(meta,indent=2),encoding='utf8')
(r/'FIRST_BATCH_REVIEW.md').write_text('''# 单细胞第一组重绘\n\n当前已完成细胞数量横版和竖版柱形图，使用原有正式计数表，细胞总数258120，包含13类明确注释和1类上皮髓系混合注释。\n\n图注：各柱表示原有单细胞注释对象中对应细胞群的数量。横版同时标注占总细胞数的百分比。各群采用固定颜色，混合注释群标灰。细胞数量仅为描述性统计。\n\n结果：T细胞63071个（24.43%），上皮细胞62652个（24.27%），巨噬细胞42546个（16.48%），为数量最多的三类。混合注释918个（0.36%）单独保留。\n\n待完成：按参考图重绘细胞类型UMAP及13类代表Marker UMAP。自动审批阻止从服务器向此目录传输逐细胞数据，也阻止仅传输图件和无标识汇总结果。需要用户明确授权图件和汇总结果接收目录后继续；逐细胞表达、细胞ID和患者样本元数据保持在服务器。\n''',encoding='utf8')
print(json.dumps(meta))
