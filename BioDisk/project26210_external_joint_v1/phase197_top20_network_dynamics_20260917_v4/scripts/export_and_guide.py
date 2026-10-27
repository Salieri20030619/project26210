from pathlib import Path
import subprocess
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase197_top20_network_dynamics_20260917_v4'
guide="""# Top20通路与Fig3非热图版 v4

本次修改Fig2与Fig3，Fig1完整沿用v3。继续独立交付，不插入正文。PDF/SVG用于排版，PNG为600 dpi，JPG用于预览。

## Fig2A：Top20选择规则

从原85通路的患者×程序×TYMP分组AUCell均值出发：

1. 在同一患者、同一cNMF程序中计算TYMP+减TYMP-的AUCell差值；沿用每组至少10个细胞的原纳入条件。
2. 同一患者如有多个合格程序，先对这些程序的差值等权平均，确保患者只贡献一个值。
3. 对16名患者的差值取中位数，按中位数绝对值从高到低选前20。
4. 保留实际正负方向，不按P值筛选，不强制加入预设通路。该次Top20实际均为正向。

点图仍采用原患者等权平均AUCell着色、AUC>0比例决定点面积；各程序合格配对患者数仍是6、9、0、2、3、2、2。排名的16名患者与每个程序的样本数不能混为一谈。F3无合格配对，以灰叉表示。各程序面板为描述性显示，不产生新的P值或FDR。

精氨酸/脯氨酸代谢列第9，嘧啶代谢第16；半胱氨酸/蛋氨酸代谢未进入此Top20，完整85条排名与方向保留于TOP20_RANKING_ALL85.tsv。AUCell差值作为本次展示筛选指标，不等同于代谢通量。

Fig2B/C/D保留原恶性上皮轨迹、MTA生化背景及C4发现。E用原255例组织数据重绘为短横版，患者顺序与数值不变，避免组织热图占据半张主图。MTA—组织SPP1原P=0.02926、FDR=0.09755，仍按名义显著候选联系解释。

## Fig3：原色块热图全部替换

### A｜基因—功能归属网络

四个KEGG程序为胆固醇代谢、溶酶体、吞噬体、抗原呈递。连线仅表示冻结KEGG基因集的成员关系，不表示调控方向、因果作用或基因之间相互作用。17个功能相关基因按原基因名单显示；另外7个身份/巨噬骨架基因独立列在下方，不为没有成员关系的基因强加连线。

共24个基因全部沿用v2患者内配对基因面板。节点颜色为26名合格患者中C4减其他成熟巨噬细胞的Δlog2(CPM+1)中位数，固定-3至3色阶，端色饱和由色条三角标示。该图替代原功能评分及患者基因色块热图，原26名患者的完整变化矩阵保留在phase195。四个程序是功能注释，不能据连线数量称其统一上调；吞噬体原配对检验不支持一致增强，抗原呈递原配对方向以负向为主。

### B｜功能表达等高线

保留既定SPP1、胆固醇代谢、溶酶体及抗原呈递四面板。紫色虚线为冻结C4边界，灰线为等表达线。色标独立，图为转录嵌入而非组织空间。

### C｜C4数量与RNA来源比例

原16名同文库患者的C4细胞占比中位数为7.4%，C4占捕获SPP1 RNA份额中位数为54.6%。用统一比例带显示，紫色为C4，灰色为其余巨噬细胞。这两个值分别是患者中位数，不是合并细胞后的比例；比例带不表示细胞或RNA从一种类型流向另一种类型，也不表示蛋白分泌量。逐患者贡献表仍在原分析中完整保留。

### D｜原始髓系轨迹

保留5研究原始图、根节点、C4定位及没有有效伪时间的灰色区域。主分析C4有效覆盖442/2674，敏感性1743/2674。没有重新拟合或补连轨迹。

### E｜表达动态曲线

由原冻结DYNAMIC_DISPLAY.npz直接读取40个秩位置的患者平衡表达均值，不重新拟合曲线。展示单核细胞基因组、SPP1、脂质处理基因组、溶酶体基因组、抗原呈递基因组；具体基因保存在DYNAMIC_GENE_GROUPS.json，共涉及原35基因中的24个。组曲线是在每个秩位置对所列基因的原平均log1p(CP10K)再等权平均，SPP1保留单基因值。直线连接原40个显示位置，没有加入新的平滑。

紫实线为主分析14名患者，青虚线为敏感性分析24名患者；纵轴为实际平均log1p(CP10K)，各小图独立刻度，无新置信区间或P值。主/敏感性人群不同，曲线间距不能解释为同一批患者的干预效应。原SPP1动态FDR分别为0.07177和0.03069。髓系排序中的动态变化与C4相对其他成熟巨噬群的配对差异是不同对比；例如抗原呈递沿单核至成熟髓系排序上升，并不否定C4相对其他成熟巨噬细胞的较低表达。

## 文件与来源

Fig1文件为phase196的相同版本。Fig2/Fig3分别为Main_Figure_2_top20_v4及Main_Figure_3_network_curves_v4。旧版热图、85通路全图和完整矩阵保留，不删除原始证据。服务器phase197包含本次选择规则、绘图脚本、基因集连线、曲线值及校验记录。正文未变。
"""
(d/'REVISION_GUIDE_CN.md').write_text(guide,encoding='utf8')
selection="""from pathlib import Path
import pandas as pd,json
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
d=r/'phase197_top20_network_dynamics_20260917_v4'
a=pd.read_csv(r/'phase131_sections2_5_2_6_reference_correction_remote_v1/data/PATHWAY_PATIENT_GROUP_DOT_VALUES.tsv.gz',sep='\\t')
v=a.pivot(index=['patient','factor','pathway'],columns='group',values='mean_AUC')
delta=(v.TYMP_positive-v.TYMP_negative).rename('delta').reset_index()
patient=delta.groupby(['patient','pathway']).delta.mean().unstack()
tab=pd.DataFrame({'median_patient_delta':patient.median(),'patients':patient.notna().sum(),'positive':(patient>0).sum(),'negative':(patient<0).sum()})
tab['absolute_median_delta']=abs(tab.median_patient_delta)
tab=tab.sort_values(['absolute_median_delta','median_patient_delta'],ascending=[False,False],kind='stable')
tab.insert(0,'rank',range(1,len(tab)+1));tab['selected']=tab['rank']<=20
tab.to_csv(d/'tables/TOP20_RANKING_ALL85.tsv',sep='\\t')
patient.to_csv(d/'tables/PATIENT_BALANCED_PATHWAY_DELTA.tsv',sep='\\t')
"""
(d/'scripts/select_top20.py').write_text(selection)
for f in sorted((d/'figures').glob('*.pdf')):
 subprocess.run(['pdftoppm','-r','600','-singlefile','-png',str(f),str(f.with_suffix(''))],check=True)
 subprocess.run(['pdftocairo','-svg',str(f),str(f.with_suffix('.svg'))],check=True)
 print('EXPORTED',f.name,flush=True)
