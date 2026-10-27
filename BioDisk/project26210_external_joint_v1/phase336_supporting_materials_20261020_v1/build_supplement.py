from pathlib import Path
import json,hashlib,shutil,zipfile,html
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
import pymupdf as fitz
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Paragraph,Table,TableStyle
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors

w=Path(__file__).resolve().parent;i=w/'inputs';q=w/'qc';figs=w/'figures';out=w/'delivery'
for p in [q,figs,out]:p.mkdir(exist_ok=True)
manifest=json.loads((q/'SOURCE_MANIFEST.json').read_text(encoding='utf-8'))
assert all(hashlib.sha256((i/z['name']).read_bytes()).hexdigest()==z['sha256'] for z in manifest)
def read(n):return pd.read_csv(i/n,sep='\t')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6})
validation={'source_files_verified':len(manifest),'new_models':False,'new_statistical_tests':False,'manuscript_baseline':'Phase335 v6','figures':{}}
def save(fig,n):
 fig.savefig(figs/f'{n}.pdf',dpi=350);plt.close(fig)
shutil.copy2(i/'S3_genes.pdf',figs/'S3.pdf')
rank=read('D02_rank_qc.tsv');ct=pd.read_csv(i/'D02_study_counts.tsv',sep='\t',index_col=0)
fig=plt.figure(figsize=(190/25.4,155/25.4));ax=fig.add_axes([.08,.63,.38,.29]);bx=fig.add_axes([.60,.63,.35,.29])
ax.plot(rank['rank'],rank.mean_component_silhouette,'o-',color='#467DA2',label='Mean silhouette');ax.plot(rank['rank'],rank.min_run_fraction,'s-',color='#B97B58',label='Minimum run coverage')
ax.set(xlabel='Number of programs',ylabel='Stability / coverage',xticks=rank['rank'],ylim=(0,1.05));ax.legend(fontsize=6.5,frameon=False,loc='lower left');ax.set_title('A  Stability across 20 starts',loc='left',weight='bold')
bx.errorbar(rank['rank'],rank.median_run_mse,yerr=rank.iqr_run_mse/2,fmt='o-',capsize=3,color='#4E8C7A');bx.set(xlabel='Number of programs',ylabel='Reconstruction MSE',xticks=rank['rank']);bx.set_title('B  Reconstruction',loc='left',weight='bold')
for a in [ax,bx]:a.axvspan(6.85,7.15,color='#D8BC72',alpha=.25,zorder=0)
ax=fig.add_axes([.08,.14,.83,.34]);start=np.zeros(len(ct));pal=['#7296BC','#E3AD66','#76A893','#9F86B4','#C47F85']
for col,c in zip(ct.columns,pal):
 v=ct[col].values/ct.sum(axis=1).values;ax.barh(np.arange(len(ct)),v,left=start,color=c,label=col,height=.65);start+=v
ax.set_yticks(range(7),ct.index);ax.invert_yaxis();ax.set(xlim=(0,1),xlabel='Fraction of cells assigned to each dominant program');ax.set_title('C  Study composition of all 11,950 malignant cells',loc='left',weight='bold')
for j,n in enumerate(ct.sum(axis=1)):ax.text(1.01,j,f'n={n:,}',va='center',size=6.8)
fig.legend(handles=[Patch(color=c,label=s) for c,s in zip(pal,ct.columns)],loc='lower center',bbox_to_anchor=(.5,.01),ncol=5,fontsize=6.7,frameon=False)
save(fig,'S4');validation['figures']['S4']={'cells':int(ct.values.sum()),'F1_cells':int(ct.loc['F1'].sum()),'F1_GSE207493':int(ct.loc['F1','GSE207493']),'ranks':rank['rank'].tolist()}

de=read('D03_clone_DE.tsv.gz');deltas=read('D03_patient_deltas.tsv.gz');ora=read('D03_pathway_ORA.tsv')
sel=de[de.gene.isin(deltas.gene)].sort_values(['edgeR_FDR','gene']).head(24).gene.tolist();hm=deltas[deltas.gene.isin(sel)].pivot(index='gene',columns='curated_patient_id',values='logCPM_delta').reindex(sel)
assert len(de)==15308 and hm.shape==(24,14)
fig=plt.figure(figsize=(190/25.4,186/25.4));ax=fig.add_axes([.13,.64,.77,.29]);nom=de.edgeR_PValue<.05;sig=de.edgeR_FDR<.05
cl=np.where(sig,'#A94F60',np.where(nom,'#D5A45A','#C5CDD3'));ax.scatter(de.edgeR_logFC,-np.log10(de.edgeR_PValue.clip(lower=1e-300)),s=4,c=cl,lw=0,rasterized=True,alpha=.7)
ax.axvline(0,color='#AAA',lw=.5);ax.axhline(-np.log10(.05),color='#777',ls='--',lw=.6);ax.set(xlabel='F1-high clone vs other: log2 fold change',ylabel='-log10(raw P)');ax.set_title('A  Patient-blocked differential expression (14 patients)',loc='left',weight='bold')
ax.legend(handles=[Patch(color=c,label=t) for c,t in [('#A94F60','FDR < 0.05'),('#D5A45A','P < 0.05; FDR >= 0.05'),('#C5CDD3','P >= 0.05')]],fontsize=6.5,frameon=False,loc='upper right')
ax=fig.add_axes([.17,.09,.70,.43]);im=ax.imshow(hm,cmap='RdBu_r',vmin=-3,vmax=3,aspect='auto');ax.set_yticks(range(24),sel,fontsize=6.5);ax.set_xticks(range(14),[f'P{k+1:02}' for k in range(14)],fontsize=6.5);ax.tick_params(length=0);ax.set_title('B  Paired logCPM differences: frozen 24-gene display',loc='left',weight='bold',pad=8);ax.set_xlabel('Independent patients; original identities in Data D03')
ca=fig.add_axes([.91,.16,.015,.27]);cb=fig.colorbar(im,cax=ca,ticks=[-3,0,3]);cb.ax.tick_params(labelsize=7);ca.set_title('Delta',fontsize=7)
save(fig,'S5');pd.DataFrame({'display_id':[f'P{k+1:02}' for k in range(14)],'original_patient':hm.columns}).to_csv(q/'D03_figure_patient_key.tsv',sep='\t',index=False,encoding='utf-8')
validation['figures']['S5']={'tested_genes':len(de),'patients':14,'displayed_genes':sel,'nominal_significant':int(nom.sum()),'FDR_significant':int(sig.sum()),'ORA_up_FDR':int(((ora.direction=='up')&(ora.FDR<.05)).sum()),'ORA_down_FDR':int(((ora.direction=='down')&(ora.FDR<.05)).sum())}

immune=read('D04_FORMAL_SUMMARIES.tsv');five=read('D04_FIVE_GENE_SUMMARIES.tsv');pairs=list(immune[['source','target']].drop_duplicates().itertuples(index=False,name=None));assert len(pairs)==10
def short(t):
 if t.startswith('GO:0002474'):return 'MHC-I'
 if t.startswith('GO:0002495'):return 'MHC-II'
 if t.startswith('GO:0006909'):return 'Phagocytosis'
 return {'Interferon Alpha Response':'IFN-alpha','Interferon Gamma Response':'IFN-gamma','Inflammatory Response':'Inflammation'}.get(t,t)
variants=['original_core','conditional','metabolic_immune_nonoverlap','immune_program_nonoverlap','joint_five_gene_deletion']
fig=plt.figure(figsize=(190/25.4,150/25.4))
for j,cond in enumerate(['L2_broad','L3_strict']):
 # Resolve actual strict condition from frozen table, without changing membership.
 if cond not in set(immune.condition):cond=next(x for x in immune.condition.unique() if x!='L2_broad')
 arr=np.full((10,5),np.nan)
 for row,(source,target) in enumerate(pairs):
  for col,var in enumerate(variants):
   table=five if col==4 else immune;sub=table[(table.condition==cond)&(table.source==source)&(table.target==target)&(table.variant==var)]
   if len(sub):assert len(sub)==1;arr[row,col]=sub.iloc[0].median_rho
 ax=fig.add_axes([.24+j*.38,.23,.34,.65]);cm=plt.get_cmap('RdBu_r').copy();cm.set_bad('#ECEFF2');im=ax.imshow(arr,cmap=cm,vmin=-.4,vmax=.4,aspect='auto')
 ax.set_yticks(range(10),[a+' / '+short(b) for a,b in pairs] if j==0 else ['']*10,fontsize=7);ax.set_xticks(range(5),['O','C','M','I','5G'],fontsize=8);ax.tick_params(length=0);ax.set_title(('A  Broad' if j==0 else 'B  Strict')+' | n=11',loc='left',weight='bold',pad=7)
 for row in range(10):
  for col in range(5):ax.text(col,row,'NA' if np.isnan(arr[row,col]) else f'{arr[row,col]:.3f}',ha='center',va='center',size=6.5,color='white' if abs(arr[row,col])>.25 else '#283845')
 ax.set_xticks(np.arange(-.5,5),minor=True);ax.set_yticks(np.arange(-.5,10),minor=True);ax.grid(which='minor',color='white',lw=.7);ax.tick_params(which='minor',length=0)
ca=fig.add_axes([.39,.15,.42,.025]);cb=fig.colorbar(im,cax=ca,orientation='horizontal',ticks=[-.4,0,.4]);cb.set_label('Median within-patient partial Spearman rho',size=7)
fig.text(.03,.035,'O: baseline   C: control the other gene   M: metabolic/immune overlap removed\nI: immune-program overlap removed   5G: joint deletion of five influential genes (descriptive)',size=7,linespacing=1.7)
save(fig,'S6');validation['figures']['S6']={'formal_rows':len(immune),'five_gene_rows':len(five),'all_relations':10,'joint_deletion_P_not_invented':bool(five.P.isna().all())}

sp=read('D05_PATIENT_EFFECTS.tsv');ss=read('D05_PROGRAM_CHECK_SUMMARY.tsv');ident=read('D05_IDENTIFIABILITY_RESULTS.tsv');wide=sp.pivot(index='patient',columns='variant',values='rho');assert len(wide)==23
fig=plt.figure(figsize=(190/25.4,154/25.4));ax=fig.add_axes([.12,.55,.33,.38]);
for row in wide.itertuples():ax.plot([0,1],[row.original,row.no_LGALS3],'-o',lw=.65,ms=2.8,color='#8C9CA6',alpha=.7)
med=wide[['original','no_LGALS3']].median();ax.plot([0,1],med,'-o',lw=2,color='#A94F60',ms=6,label='Patient median');ax.axhline(0,color='#888',ls='--',lw=.6);ax.set_xticks([0,1],['Original','Without LGALS3']);ax.set(ylabel='Partial Spearman rho',xlim=(-.25,1.25));ax.set_title('A  Same 23 patients',loc='left',weight='bold');ax.legend(frameon=False,fontsize=7)
ax=fig.add_axes([.65,.55,.28,.38]);pop=['ffpe','frozen'];x=np.arange(2)
for k,(var,label,col) in enumerate([('original','Original','#668DAA'),('no_LGALS3','Without LGALS3','#A94F60')]):
 z=ss[ss.variant==var].set_index('population').loc[pop];ax.bar(x+(k-.5)*.28,z.median_rho,width=.28,color=col,label=label)
 for j,row in enumerate(z.itertuples()):ax.text(x[j]+(k-.5)*.28,row.median_rho+.005,f'{row.positive}/{row.n}',ha='center',fontsize=7)
ax.set(xticks=x,xticklabels=['FFPE','Frozen'],ylabel='Median rho',ylim=(0,.27));ax.set_title('B  Preparation strata',loc='left',weight='bold');ax.legend(frameon=False,fontsize=6.5)
idc=ident[ident.variant=='complete'].copy();ax=fig.add_axes([.25,.12,.68,.29]);names=['C4','Other_mature_macrophage','CD8_T','NK','gdT_enriched','Malignant'];idc=idc.set_index('label').loc[names];bars=ax.barh(range(6),idc.median_margin,color=['#A94F60' if not v else '#518B80' for v in idc.gate_pass]);ax.set_yticks(range(6),['C4','Other macrophages','CD8 T','NK','gdT-enriched','Malignant'],fontsize=7);ax.invert_yaxis();ax.axvline(.05,ls='--',c='#333',lw=.8);ax.set_xlabel('Median classification margin');ax.set_title('C  Frozen reference identifiability',loc='left',weight='bold',pad=8)
for j,row in enumerate(idc.itertuples()):ax.text(row.median_margin+.002,j,f'{int(round(row.recall*row.n_patient_state_units))}/{row.n_patient_state_units}',va='center',size=7)
ax.set_xlim(0,.15);save(fig,'S7');validation['figures']['S7']={'patients':23,'original_median':float(med.original),'no_LGALS3_median':float(med.no_LGALS3),'C4_margin':float(idc.loc['C4','median_margin']),'complete_identifiability_rows':len(ident)}

# Chinese reading document: no new statistical estimates, source tables copied byte for byte.
pdfmetrics.registerFont(TTFont('Song',r'C:\Windows\Fonts\simsun.ttc',subfontIndex=0));pdfmetrics.registerFont(TTFont('Hei',r'C:\Windows\Fonts\simhei.ttf'))
style=ParagraphStyle('body',fontName='Song',fontSize=10,leading=15,wordWrap='CJK',spaceAfter=6)
small=ParagraphStyle('small',parent=style,fontSize=9,leading=13)
heading=ParagraphStyle('heading',fontName='Hei',fontSize=14,leading=21,spaceAfter=8)
mm=72/25.4;W,H=210*mm,297*mm;layout=w/'layout.pdf';cv=canvas.Canvas(str(layout),pagesize=(W,H));placements=[];page=0
def para(s,y,sty=style):
 p=Paragraph(html.escape(s).replace('\n','<br/>'),sty);ww,hh=p.wrap(W-80,H);assert y-hh>38,(page,s[:30],y,hh);p.drawOn(cv,40,y-hh);return y-hh-7
def begin(title):
 global page
 page+=1;cv.setFillColor(colors.HexColor('#65717D'));cv.setFont('Song',8);cv.drawString(40,H-27,'26210 ｜ 与完整稿v6配套的补充材料');cv.setFillColor(colors.black);return para(title,H-48,heading)
def end():
 cv.setFont('Song',8);cv.drawCentredString(W/2,23,f'S-{page}');cv.showPage()
y=begin('补充材料：程序定义、稳健性与分析边界')
y=para('本材料配套于已更新Fig2的完整稿v6。补充图从S3开始编号，原文中的S1、S2及原附表S1–S20保持。这里集中展示已完成但未作为主图展开的证据；完整稿中的主要阴性及关键解释边界继续保留。',y)
modules=[('S3','程序代表基因','固定35个基因及F1–F7相对载荷'),('S4','程序稳定性与来源','20次起始运行、rank 4–8、全部队列构成'),('S5','F1高主克隆差异','14患者配对、完整15308基因与402条富集检验'),('S6','免疫关联稳健性','11患者、控制另一基因、去共享基因及5基因删除'),('S7','空间定义敏感性','去LGALS3前后对照及参考可识别性')]
for num,title,desc in modules:y=para(f'{num}  {title}：{desc}。',y)
y=para('结果保留原则',y,heading)
y=para('以独立患者为统计单位，保留原效应量、原始P值、FDR及检验家族。不重新拟合模型、不重新选择患者，也不把细胞数或spot数作为独立样本量。P<0.05且FDR未通过的结果保持名义显著定位。五基因联合删除属于既有数据内选择后的描述性检查，源表没有推断P值，本材料不补造。',y)
y=para('补充数据D01–D06随附为完整统计表压缩包。新增图只改变展示方式，既有数值文件逐字节保留；图中简写、患者键和来源索引附在压缩包内。D06包含完整候选调控子靶基因表，可核查ILF2→TYMP在完整网络中的预测连接；它不属于主图2E展示的71条边。',y)
end()

pages=[
('S3','F1–F7的固定代表基因区分程序组成',116,
'每个程序沿既有载荷排序依次选取5个尚未被前面程序选中的代表基因，共35个；未为本次补充图重新挑选。横轴为F1–F7，纵轴为基因，颜色为同一基因的载荷除以其七程序最大值。两块热图依次展示F1–F4及F5–F7来源的代表基因，色条标示选择来源。',
'这张图描述的是程序组成，不是差异表达检验。TYMP不在原2500个cNMF载荷输入基因中，其RNA来自完整表达矩阵，作为独立读出与F1使用量关联；不把TYMP添加进载荷矩阵。对应数据D01。'),
('S4','rank 7沿用冻结选择，F1主导细胞集中于单一队列',155,
'A、B分别展示20次原始起始运行的稳定性指标与重建误差；阴影标示冻结rank 7，误差棒沿原图为运行MSE四分位距的一半，并非置信区间。C展示全部11,950个恶性细胞的主导程序及队列构成，比例分母为各程序细胞数。',
'F1主导细胞为4,079个，其中4,065个来自GSE207493（99.66%）。该来源集中性需与患者内连续F1相关分析区分，不能将主导群的来源组成写成跨队列均衡复现。其他rank也可能通过稳定性门槛；本图不宣称rank 7为唯一最优解。数据D02。'),
('S5','F1高主克隆具有患者配对的差异表达证据',186,
'A为14名患者、28个患者×组别pseudobulk样本的原edgeR结果，横轴log2倍数变化，纵轴原始P的负对数，三种颜色区分FDR显著、仅名义显著及未达P<0.05。B沿用既有24基因展示集合，颜色为每位患者F1高主克隆减其他主克隆的logCPM差。',
'完整表保留15,308个受检基因，原FDR<0.05共760个（上调233、下调527）。原方向性富集以各方向按P排序的前250基因为输入，上调17/201条FDR<0.05、下调0/201条。完整两方向结果与患者标识见D03；不把载荷、群内相关和克隆差异表达互换。'),
('S6','免疫程序的稳健性具有明确的关系依赖',150,
'两块热图分别为既有宽松与严格成员版本，均以11名独立患者为单位；全部10条预定关系均保留。单元格显示患者内偏Spearman相关的中位数，灰色NA代表该组合原先未定义，不能当作零。O为继承基线，C控制另一基因，M去代谢/免疫重叠，I去免疫程序重叠，5G为联合删除5个影响基因。',
'TYMP–IFNγ/MHC-I/吞噬、SPP1–IFNγ/MHC-II主要方向在核查后保持。IFNα在免疫程序去重叠后减弱；TYMP–炎症在宽松版本的相应rho降至0.003014，P=0.548829、FDR=0.564509。宽松TYMP–糖酵解从条件rho=0.102671降至5G描述值0.021528。原224、条件20、去重叠36检验家族分别保留，5G不补P值。全部76行正式统计、20行5G汇总及逐患者/逐基因结果见D04。'),
('S7','去除共享LGALS3后空间程序关联明显减弱',154,
'A以相同23名患者配对展示原程序与双程序同时去除LGALS3后的效应；一名患者的两种制备按原Fisher z等权合并。B分FFPE与frozen展示，各12名患者，柱顶为正向人数/有效人数，两制备共享1名患者。C展示原参考可识别性检查中的六类关键身份，数字为识别正确数/患者状态单位数，虚线为预定0.05分离余量，红色为未通过整体门槛。',
'FFPE原/去共享rho为0.21206/0.01875，frozen为0.20835/0.01793。去共享P=0.17627/0.33936，两项新核查家族FDR均0.33936；原程序FDR分别0.0021973/0.0043945，沿原家族保留。23患者汇总去共享P=0.104543，此单个预设汇总不另作FDR。C4虽26/26分类正确，分离余量0.02938未过门槛；不能据此声称混合spot的C4丰度已验证。完整16类×3基因版本及全部患者效应见D05。')]
for n,title,height,cap,interpretation in pages:
 y=begin(f'补充图{n}  {title}');y-=5
 placements.append({'figure':n,'page':page-1,'rect':[10*mm,H-y,(200)*mm,H-y+height*mm]});y-=height*mm+12
 y=para(cap,y,small);y=para(interpretation,y,small);end()
y=begin('补充数据索引与阅读口径')
for label,desc in [('D01','35×7基因相对载荷矩阵。'),('D02','rank 4–8的全部稳定性/重建指标，以及七程序×五队列的全部细胞计数。'),('D03','15,308基因的edgeR与配对Wilcoxon结果、2,100行患者差值、402条方向性富集，以及图中患者简写对应。'),('D04','全部76行正式效应、836行患者效应、单基因剔除与5基因联合剔除完整结果、基因重叠和保留说明。源表保留原检验家族，5G中的P空值为不做该推断检验。'),('D05','原/去LGALS3程序的逐片和逐患者效应、汇总统计、全部参考覆盖与可识别性结果。各制备和患者汇总分开，空FDR不当作不显著或零。'),('D06','GENIE3/GRNBoost2完整候选调控子靶基因表，包含权重和算法来源；预测边不等同于实验确认。')]:y=para(label+'：'+desc,y)
y=para('与正文的对应',y,heading)
y=para('S3–S5补充正文2.1–2.2及图2的程序/克隆证据；S6补充免疫程序章节的完整稳健性结果；S7补充空间程序定义和参考可识别性边界。主文已保留的数值与关键阴性不因进入补充材料而删除。仅将选定模块制成配套材料，未把历史失败版本、无关探索或原始表达矩阵加入交付包。',y)
y=para('文件采用UTF-8编码；TSV列名及原数值保持。源文件清单、逐文件SHA256和图件显示数据存于服务器Phase336目录，可按表中原字段回溯。当前材料只整理既有分析，不构成新增独立队列验证。',y)
end();cv.save()
doc=fitz.open(layout)
for pl in placements:
 src=fitz.open(figs/f"{pl['figure']}.pdf");doc[pl['page']].show_pdf_page(fitz.Rect(pl['rect']),src,0)
doc.set_toc([[1,'补充材料概览',1]]+[[1,f'补充图{z[0]} {z[1]}',k+2] for k,z in enumerate(pages)]+[[1,'补充数据索引',7]])
final=out/'26210_Supplementary_Materials_S3-S7_v1.pdf';doc.save(final,garbage=4,deflate=True)
for k,p in enumerate(doc):p.get_pixmap(matrix=fitz.Matrix(1.4,1.4)).save(q/f'page_{k+1:02}.png')
readme='26210 补充数据D01–D06\n与完整稿v6及补充材料S3–S7配套。\n所有D开头的原始统计文件按来源逐字节复制；没有重新计算P或FDR。\nD03_figure_patient_key.tsv为绘图患者简写映射。\nsource_manifest.json记录服务器来源、文件大小和SHA256；Q开头内部核查文件不进入此统计交付包。\n详细检验家族与空值含义见补充材料PDF。\n'
with zipfile.ZipFile(out/'26210_Supplementary_Data_D01-D06_v1.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in sorted(i.glob('D*')):z.write(p,p.name)
 z.write(q/'D03_figure_patient_key.tsv','D03_figure_patient_key.tsv');z.writestr('阅读说明.txt',readme.encode('utf-8'));z.writestr('source_manifest.json',json.dumps([x for x in manifest if x['name'].startswith('D')],ensure_ascii=False,indent=2).encode('utf-8'))
 z.writestr('D05_FDR_family_notes.txt',('原程序FDR沿原Phase70家族保留：FFPE 0.0021973；frozen 0.0043945。\n去LGALS3两制备为Phase219新核查家族，两项FDR均0.33935546875。\n23患者合并汇总是单个预设检查，不另作FDR。\nD05_PROGRAM_CHECK_SUMMARY.tsv中空FDR未填零，也未重新校正。\n原FDR来源：'+str(next(x['source'] for x in manifest if x['name']=='Q_spatial_notes.md'))+'，第三节。\n').encode('utf-8'))
validation.update({'pages':len(doc),'placements':placements,'pdf_sha256':hashlib.sha256(final.read_bytes()).hexdigest(),'visual_review':'PENDING'})
(q/'VALIDATION.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2),encoding='utf-8')
print('BUILT',len(doc),'pages',validation['figures'])
