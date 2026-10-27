from pathlib import Path
import json,subprocess,shutil
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase196_selected_panels_20260917_v3'
text="""# 三张主图逻辑重排 v3

本次按用户指定图件修订Fig1和Fig2，继续独立输出，未加入正文。Fig3保留phase195扩充版，文件内容不变。

## Fig1：细胞背景、程序结构与调控关系

A，细胞图谱；B，代表患者CNV；C，cNMF程序、F1和TYMP；D，用户指定的患者—程序表达分块图；E，患者平衡调控子活性；F，用户此前指定的调控子—靶基因网络；G，紧凑蛋白活性密度图。

D直接取自H_patient_program_expression_blocks_v2.pdf，替换原来的基因贡献/F1排序热图组合。保留树状图、分块、相关色标、队列标签和程序标签。原图含19名患者的31个患者×程序表达谱，每个表达谱至少20个细胞，以2,500个可变基因构建标准化平均RNA表达关系。31个表达谱不是31名独立患者；显示簇C1—C7不是重新定义的单细胞群。未重排原矩阵或重新聚类。原图H面板字母由新主图D替代。

F仍为71节点、71有向边的原网络，G仍为原患者等权KDE矩阵。网络为预测调控边，蛋白活性为推断值。F1主导细胞原来源集中、PXN未达显著等关键结果仍见原报告及phase195图注。

## Fig2：完整代谢背景到SPP1状态发现

### A｜左侧完整一列：85条代谢通路

取用户指定03_all85_pathways_TYMP_dotplot_reference_v2.pdf对应的原始PATHWAY_DOTPLOT_SUMMARY.tsv及原绘图函数，按长列重新设置画布、字体和点面积的统一显示比例。全部85条通路顺序、TYMP-/TYMP+分组、F1—F7列、颜色数值、点面积对应的检测比例均保留，未删选通路或突出单独阳性结果。

颜色为患者等权平均AUCell评分，点面积为AUC>0的细胞比例，灰叉表示没有符合条件的配对样本，点大小不是P值。同一程序的TYMP-/TYMP+使用相同患者，要求每名患者每组至少10个细胞；原F1—F7合格患者数分别为6、9、0、2、3、2、2，故此图为描述性完整通路背景，不把细胞数或点数作为独立重复。

### B—E｜右侧按逻辑顺序排列

B，恶性上皮原始Monocle3轨迹，承接上游转录状态排序。保留原11,950细胞嵌入、原主图连线及色标，不重拟合轨迹，不将横断面排序视为实际谱系转化。

C，MTA生化反应背景，说明蛋氨酸/SAM和多胺合成与MTA的联系。实线仍为已知生化反应，不新增TYMP直接催化MTA的因果箭头。

D，成熟巨噬细胞图谱及C4定位，先展示独立的状态发现。C4由冻结无监督群确定，不能将其等同于全部SPP1阳性细胞。

E，4队列255例实测MTA与组织表达图谱，在状态发现之后呈现组织层面的联系。保留全部病例、队列大小以及按MTA排序的原顺序；原MTA—SPP1 P=0.02926、FDR=0.09755，为名义显著候选联系。配对组织数据不等同于C4内的MTA测量。

为减少重复热图，原Fig2两种恶性上皮表达动态热图及C4标志热图不再占用本主图。原件和方法间差别仍保留在phase192、phase194及原报告，未删除分析结果；C4身份基因的患者证据继续见Fig3配对基因图。

## Fig3：沿用扩充版

保留phase195的六组面板：功能程序、26名患者24个功能相关基因的配对变化、功能表达等高线、16名患者的C4细胞占比与SPP1 RNA贡献、原始髓系分区轨迹、患者平衡表达动态。文件按SHA256核对为完全一致，不重生成。

## 交付

目录包含Fig1/Fig2新版PDF、SVG、600 dpi PNG及JPG预览，以及完全沿用的Fig3 v2同格式文件。Fig1宽340 mm；Fig2为适合长列排布的420 mm宽矢量母版，保留通路文字以便放大审阅。源PDF、矩阵、原绘图函数引用、显示参数和面板来源表留在服务器phase196目录。
"""
(d/'LAYOUT_GUIDE_CN.md').write_text(text,encoding='utf8')
for f in sorted((d/'figures').glob('*.pdf')):
 subprocess.run(['pdftoppm','-r','600','-singlefile','-png',str(f),str(f.with_suffix(''))],check=True)
 subprocess.run(['pdftocairo','-svg',str(f),str(f.with_suffix('.svg'))],check=True)
 print('EXPORTED',f.name,flush=True)
