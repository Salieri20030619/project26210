from pathlib import Path
import json,hashlib,shutil,subprocess,concurrent.futures
import pandas as pd,numpy as np,fitz
r=Path.cwd();d=r/'phase206_all_mainfigures_A4_20260918_v10';out=d/'delivery'
notes='''# 三张主图 A4 单页版 v10

三张图均为210 × 297 mm，每张占一页，横向有效范围约190 mm。打印或插入文档请用100%/实际大小。版式和字大小已按最终页面制作；PDF可提取文字最小6 pt，面板编号10 pt，短标题7.5 pt。提供矢量PDF/SVG及600 dpi PNG。

## Fig1：TYMP相关恶性上皮程序

A：原258,120个细胞、75个文库、40名患者的整合UMAP，保留原注释与坐标。颜色表示细胞类别；图中四大类及主要细胞图例为读图简化，红细胞为棕色，上皮/髓系混合为灰色。

B：GSE156632 RCC2的原染色体CNV示例。列按染色体1–22排列，主图仅稀疏显示1、5、10、15、20号标签以避免拥挤，不删除其余染色体。上方200个参考细胞，下方701个肿瘤上皮查询细胞；颜色为原log2相对信号。

C：同一原恶性上皮UMAP，分别显示七个cNMF主导程序、F1 usage及TYMP RNA。未重新拟合cNMF、UMAP或表达值；F1–F7颜色沿用原图。

D：31个患者×程序表达profile之间的Spearman相关，来自19名患者；每个profile至少20个细胞，原2,500变异基因及基因内标准化平均RNA保持。保留原矩阵、排序、树和C1–C7展示簇。P01–P31是profile短索引，不是31名独立患者，完整患者/研究/程序对应见索引表。细条分别对应原程序及研究分组。色条为Spearman相关。

E：原24个regulon、七个程序、两种算法的患者平衡AUCell矩阵，顺序、值及原行树保持。每个程序的相邻两列为G3（GENIE3）和GB2（GRNBoost2）；各算法内按行Z标准化。21名独立患者形成37个患者×程序profile，F1–F7各自有效患者数9、12、2、4、4、2、4。粗体保留ILF2、HOXD3及HNF4A。行树的四个展示簇沿原定义，不是新细胞亚型。

F：原3个TF和68个靶基因，71个节点及71条有向连接全部保留。菱形为TF，灰点为靶基因，金色为共享靶点；紫线为两算法共同支持，棕线仅GENIE3，蓝线仅GRNBoost2。网络为预测调控连接，不是新的因果推断。

G：F1 usage与四个蛋白metaVIPER NES的原密度网格。保留CASP9、CHMP1A、PRSS8及PXN；色深表示相对密度。未重新估计KDE，统计单位为患者；完整原效应量、n、P和FDR附于统计表。

## Fig2：代谢背景、巨噬候选与SPP1来源

A：原Top20通路及排序保留，列为TYMP-/TYMP+中的F1–F7。颜色是同一通路两组全部合格程序均值共同计算的行Z，不能解释为跨通路绝对活性或组间显著性。点面积为AUCell>0的细胞比例。Top20排序使用16名患者；F1–F7合格配对患者数6、9、0、2、3、2、2，无合格配对的F3以灰叉表示。通路长名作简写，原名和所有展示值保留在统计数据中。

B：原11,950个恶性上皮细胞及224段Monocle3主图。使用原UMAP及原全局伪时间，未重拟合。颜色从Early到Late。

C：F1及嘧啶、半胱氨酸/甲硫氨酸、精氨酸/脯氨酸及甘油磷脂（GPL）代谢之间的关系。下三角为20例患者内Spearman rho的中位数，圆面积为|rho|；上三角显示rho>0的患者数/20，条长表示比例。F1四通路原P/FDR依次为：嘧啶0.0145/0.0193、Cys/Met 0.00124/0.00248、Arg/Pro 0.000636/0.00248、GPL 0.49/0.49；沿用原四通路校正家族。完整非四舍五入统计见附表。

D：将原两幅巨噬定位展示合成一幅，保留12个冻结候选的密度底图，并叠加C4紫色虚线、C10蓝色点划线、C15绿色点线。边界为各群自身密度峰值的12%，表示群位置而非SPP1峰值或专属表达。细胞成员、坐标和密度算法未改变。

E：中心为SPP1 log2(CPM+1)患者中位数；浅色内环为SPP1阳性细胞比例中位数；深色外环为SPP1原始RNA贡献比例中位数，下方依次显示内/外百分比。各圆统一100%，从12点顺时针增加。完整发现队列的表达/阳性率有效患者数分别为C4 26、C10 2、C15 9、联合28、整体38；群及联合要求每例目标与其余各至少20个巨噬细胞，整体至少20个巨噬细胞。外环使用原16例同文库匹配患者，分母为每例全部成熟巨噬细胞的SPP1原始计数，缺失群记0。C10的0%是中位数，不表示每例均为0。组之间有包含关系且有效患者构成不同，不能相加或把中位数之差当作配对效应；仅作描述，不进行新的推断检验。

F：原255例组织样本，四队列分别97、52、50、56例。研究内按实测MTA排序，MTA先log2(x+1)，所有行在各队列内Z标准化。该图为组织层面，不是C4特异结果。原MTA–SPP1关系P=0.02926、FDR=0.09755，为名义显著。

## Fig3：C4状态、功能富集与髓系顺序

A：原38个marker×12个成熟巨噬状态，患者等权log2(CPM+1)汇总后在原12状态间按行Z，C4列以紫框标记。基因、状态、顺序及数值不变。

B：原SPP1、胆固醇代谢、溶酶体及抗原呈递四个局部平均表达/分数网格，四个数值矩阵与phase193逐元素一致。色阶仍使用各场P2–P98，灰色为等值线。三条彩色群位置轮廓与Fig2D一致，不能把C4轮廓当SPP1最高表达区域。RNA单位为log1p(CP10K)，程序使用原分数口径。

C：仍为原C4对其他成熟巨噬细胞GSEA，26名配对患者、12,048个连续配对t统计量排序基因，正排序表示C4较高；BH校正在原299条合格KEGG通路家族内进行，展示原预设四功能通路。完整NES、P和FDR见下方统计表，未改成三群联合分析。三轨分别为富集运行曲线、命中基因位置、排序统计量；未重新运行GSEA。

D：原五研究、原分区髓系轨迹，紫色轮廓为C4，星号表示原单核根，灰色为无有限伪时间。颜色为各分区内归一化伪时间；保留未连通分区，未强制连接C4。研究间不能当作同一连续时间轴。

E：原五类动态表达曲线，实线为主分析14例，虚线为敏感性24例。横轴为分区内伪时间秩0–1，纵轴为原平均log1p(CP10K)。保留单核、SPP1、脂质处理、溶酶体和抗原呈递原基因集及40点原曲线，不重新平滑或拟合。小样本结果按探索性定位。

## 排版说明

图内仅保留面板编号、短标题、必要坐标/刻度、群名、图例和关键数值。方法、n口径、P/FDR和解释集中在此图注与统计附件。Fig2保留原B/C编号，因此右上为C，右中为B；所有面板均以字母引用。本轮只调整图形呈现，C4继续作为后续主分析，旧报告正文和原始结果不覆盖。
'''
# Source statistics retained, not recomputed.
statdir=out/'统计与索引';statdir.mkdir(exist_ok=True)
for src,name in [(r/'phase199_fig3_marker_GSEA_20260918_v5/tables/GSEA_FOUR_PROGRAMS_primary.tsv','Fig3C_GSEA统计.tsv'),(r/'phase202_candidate_labels_comparison_20260918_v6/tables/CANDIDATE_VS_WHOLE_SUMMARY.tsv','Fig2E_双环统计.tsv'),(r/'phase198_fig2_contrast_corr_revision_20260918_v5/tables/TOP20_DISPLAY_VALUES.tsv','Fig2A_完整通路名与数值.tsv')]:shutil.copy2(src,statdir/name)
a=pd.read_csv(r/'phase195_mainfigures_compact_20260917_v2/tables/ORIGINAL_PROTEIN_STATISTICS.tsv',sep='\t');a[a.protein.isin(['CASP9','CHMP1A','PRSS8','PXN'])].to_csv(statdir/'Fig1G_原患者统计.tsv',sep='\t',index=False)
src=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/09_metabolism_trajectory/metabolism_v1/tables/four_pathway_AUCell_correlation_summary_v1.tsv');tab=pd.read_csv(src,sep='\t');tab[tab.predictor=='F1'].to_csv(statdir/'Fig2C_F1原P与FDR.tsv',sep='\t',index=False)
keys=json.loads((d/'tables/PATIENT_LABEL_KEY.json').read_text());assert len(keys)==31
keymd='# Fig1D 患者×程序profile索引\n\nP01–P31表示31个profile，涉及19名独立患者，不能以31作为生物学样本量。\n\n| 索引 | 原研究、患者、程序 |\n|---|---|\n'+'\n'.join(f"| {x['display']} | {x['original']} |" for x in keys)
(statdir/'Fig1D_profile索引.md').write_text(keymd)
g=pd.read_csv(statdir/'Fig3C_GSEA统计.tsv',sep='\t');notes+='\n## Fig3C 完整GSEA统计\n\n| 通路 | NES | P | FDR |\n|---|---:|---:|---:|\n'+'\n'.join(f"| {z['ID']} | {z['NES']:.6g} | {z['pvalue']:.8g} | {z['p.adjust']:.8g} |" for z in g.to_dict('records'))+'\n'
(out/'三张主图图注与排版说明_v10.md').write_text(notes)
checks=json.loads((d/'qc/FONT_PAGE_AUDIT.json').read_text());assert all(z['min_font']>=5.95 and not z['small'] and not z['outside_margin'] for z in checks)
for n in [1,2,3]:
 p=d/'figures'/f'Main_Figure_{n}_A4_v10.pdf';doc=fitz.open(p);assert len(doc)==1;assert np.allclose([doc[0].rect.width/72*25.4,doc[0].rect.height/72*25.4],[210,297],atol=.01)

def export(n):
 p=d/'figures'/f'Main_Figure_{n}_A4_v10.pdf';subprocess.run(['pdftocairo','-svg',str(p),str(p.with_suffix('.svg'))],check=True,stderr=subprocess.PIPE);subprocess.run(['pdftoppm','-r','600','-singlefile','-png',str(p),str(p.with_suffix(''))],check=True,stderr=subprocess.PIPE)
 for ext in ['pdf','svg','png']:shutil.copy2(p.with_suffix('.'+ext),out/p.with_suffix('.'+ext).name)
 return n
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 for n in pool.map(export,[1,2,3]):print('EXPORTED',n,flush=True)
manifest=[dict(path=str(p.relative_to(out)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in out.rglob('*') if p.is_file() and p.name!='FILE_MANIFEST.json'];(out/'FILE_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2));print('PACKAGE_READY',len(manifest),flush=True)
