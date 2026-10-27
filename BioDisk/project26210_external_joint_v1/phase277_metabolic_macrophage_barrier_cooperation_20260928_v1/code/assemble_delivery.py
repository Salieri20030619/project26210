from pathlib import Path
import pandas as pd,numpy as np,json,datetime,hashlib,shutil,zipfile
P=Path(__file__).resolve().parents[1];H=P/'harmonized_platform_v1';D=P/'delivery'
s=pd.read_csv(P/'tables/COHORT_SUMMARY.tsv',sep='\t');h=pd.read_csv(H/'tables/COHORT_SUMMARY.tsv',sep='\t');identity=pd.read_csv(P/'tables/C4_IDENTITY_SUMMARY.tsv',sep='\t')
def fmt(x):
 if pd.isna(x):return '—'
 return x if isinstance(x,str) else f'{x:.6g}'
def table(d,cols,labels):
 out=['| '+' | '.join(labels)+' |','| '+' | '.join(['---']*len(cols))+' |']
 for _,r in d.iterrows():out.append('| '+' | '.join(fmt(r[k]) for k in cols)+' |')
 return '\n'.join(out)
cols=['endpoint','model','n','positive','negative','median','P','FDR']
labels=['端点','模型','患者n','正向','负向','中位效应','主P','FDR']
text='''# 26210：代谢相关巨噬状态与免疫屏障协同检验（Phase277）

2026-09-28。研究主线按用户最新意见调整：屏障过程做有限文献对齐，重点评价代谢相关巨噬状态参与并调节免疫屏障。不要求SPP1单独构成或独立解释全部屏障。本轮完成统计与验证，不绘图、不修改当前Phase273正式稿。

## 本轮形成的主要认识

**C4与APOE高表达、溶酶体增强的代谢相关巨噬状态具有明确对应；现有空间数据尚未支持C4进一步增强CAF/ECM相关CD8排斥。** 这两个层次应分开陈述。

1. 单细胞端24位配对患者中，C4的APOE及去共享基因溶酶体程序均24/24高于其他巨噬。该对应为本项目状态与文献APOE巨噬轴提供衔接，不等于两种状态完全相同。
2. C4胆固醇稳态程序在去共享后22/24降低，与原Hallmark结果方向一致。应描述“脂质处理/溶酶体相关状态重塑”，不能写成“胆固醇代谢或所有代谢整体增强”。
3. 跨平台统一后，C4状态与TYMP、多胺背景、胆固醇程序存在未调整的空间正相关，但深度/组成调整后减弱，尚不能归因为独立代谢驱动。
4. C4×CAF/ECM交互没有通过预设原P<0.05标准，尤其CD8结果方向不稳定。文献过程对照中的APOE×POSTN总T负交互为名义显著，但条件斜率仍多为正，不能直接解释为已复现协同屏障。

## 从状态身份到代谢特征

固定原C4成员，使用原24位肿瘤患者C4/其他巨噬配对伪bulk；未根据空间结果重分群或选基因。APOE单基因、文献APOE巨噬程序去除C4共享基因MARCO/S100A9，以及去共享后的胆固醇稳态/溶酶体程序为4个冻结端点。

功能程序来自缓存的Hallmark与KEGG。为避免同一基因同时定义状态和功能，从功能程序剔除C4、APOE巨噬、POSTN-CAF、ECM及T/CD8读出的共享基因，另排除SPP1、TYMP。方向为C4减其他巨噬的平均log2(CPM+1)差。患者为独立n，精确符号检验为主、Wilcoxon为伴随，4项一组BH。

'''+table(identity,['program','n','positive','negative','median','P','FDR','P_wilcoxon','FDR_wilcoxon'],['读出','患者n','升高','降低','中位差','主P','FDR','Wilcoxon P','对应FDR'])+'''

APOE与溶酶体的升高在5个研究中方向均一致；胆固醇稳态降低也在5研究均保持。APOE整套去共享程序仅16/24位升高：主P=0.1516，伴随Wilcoxon P/FDR=0.02486，说明整体文献状态对应受统计方法影响，不能据APOE单基因高表达直接把C4改名为文献APOE巨噬。

原Phase200完整基因家族统计随INHERITED_*表保留：C4肿瘤配对Hallmark胆固醇稳态原中位差−0.03031，P=0.03832、FDR=0.04077；Lysosome中位差+0.33488，原P约1.94e−5、FDR约4.6e−5。本轮去共享检验是对已有状态解释的针对性核查，不包装成独立新队列发现。

## 一个已查明并处理的平台问题

完整C4状态签名采用Phase271冻结的19基因，不含SPP1/TYMP。在GLMF两片与GSE175540 12张FFPE切片中，CES1和FABP3不存在于原始feature列表；不是低表达被补零，也不是QC阶段丢失。已逐一核对14个原始文件及SHA。

因此保留原19基因分析全部结果，另在模型运行前冻结一个按38片共同测量基因确定的技术扩展。去除这两个未测基因后，所有平台统一使用17基因；不根据效果保留不同基因。四种评分定义如下：

| 版本 | 定义 | 用途 |
|---|---|---|
| C4_19 | 原固定19基因 | 原方案主分析，FFPE不具备完整覆盖 |
| C4_17 | 原19基因去MARCO/S100A9 | 原方案去APOE共享核验，仍受CES1/FABP3缺失影响 |
| C4_platform17 | 原19基因去CES1/FABP3 | 平台一致性技术扩展 |
| C4_platform15_noAPOEshared | 平台17再去MARCO/S100A9 | 技术扩展的共享基因敏感性 |

统一17/15基因均在38片满足评分资格；单细胞配对中均24/24位C4评分更高，17基因中位差+2.6093。参考集的5个研究中，17基因评分的细胞类型中位数均以C4最高。该核查基于既有发现资料，不能冒充新的留出验证；参考中成纤维细胞覆盖有限，评分仍不能当作空间C4细胞比例。

原方案102项空间端点完整保留；技术扩展新增66项，以独立登记的技术扩展家族校正。原方案与技术扩展均披露，不用扩展结果替换原主分析或隐藏基因覆盖失败。技术附录见harmonized_platform17/。

## 空间检验如何对应“参与并调节”

沿用38片/37患者的冻结表达分区及Phase276的400 µm边界弧段，不重选边界、切片或阈值。每个边界段用200 µm源带均值表示C4、巨噬和基质RNA程序，要求源带至少10个spot且平均总巨噬比例≥0.05；免疫结局是400 µm邻域内侧减外侧CD8谱系RNA或总T log2映射丰度。

连续交互模型：
结局 = 状态主效应 + 基质主效应 + 状态×基质 + 协变量。
C4/巨噬状态和基质变量在片内按合格边界段标准化。交互β代表超过两项相加关系的效应变化，目标不是任一协变量的最小P值。负交互只有在对应条件斜率/预测方向一致时才支持排斥增强，不能仅凭乘积项负号宣称机制成立。

核心调整固定包括深度均值与内外差、总巨噬丰度、恶性比例均值与内外差、内外spot数量比。C4状态调整再加入SPP1份额、去共享APOE程序及内皮比例，评价是否超出一般巨噬状态；文献过程扩展调整只增加内皮比例。至少50个合格边界段、四个中位数象限每个至少5段、满秩和标准化设计条件数<30，不对边界段作独立生物学P检验。

同患者多切片先合并；每个队列单独按患者方向做主检验，n<10仅描述。整个过程保留缺失、资格不足、反向及全部模型。没有因FDR未通过删除原P<0.05信号。

## 主体：C4是否改变CAF/ECM相关的免疫内外差

平台17基因版本的GSE175540结果：

'''+table(h[(h.dataset=='GSE175540')&(h.family=='C4_primary_interaction')],cols,labels)+'''

C4×POSTN-CAF对总T的核心调整交互，19位患者14负5正，中位β=−0.01960，P=0.06357/FDR=0.3850；尚未通过原P<0.05标准，不纳入名义显著候选。进一步控制SPP1/APOE/内皮后12负7正，中位β=−0.02606，P=0.3593/FDR=0.9581。去共享15基因版本两种调整均为12负7正、P=0.3593。

CD8结局核心模型9负8正，状态调整后7负10正，没有稳定预期方向。ECM替代读出也未获得支持。TJ-RCC有效n9，保持描述性；GLMF1因合格源段28<50不可估计，GLMF2合格50段可描述，但不单独进行人群推断。

留一患者中，GSE175540总T核心交互中位系数保持负值，但P在0.03088—0.09625之间变化；完整样本主P仍为0.06357，不挑选剔除后的P值作为阳性。当前数据不支持“C4增强CAF/ECM相关CD8排斥”这一统计结论。

## 代谢相关背景与C4的空间对应

'''+table(h[(h.dataset=='GSE175540')&(h.family=='metabolic_context')],cols,labels)+'''

GSE175540中，未调整的TYMP、去共享胆固醇、多胺/MTA上下文与C4程序均为15/19位正向，P=0.01921/FDR=0.03842；中位rho依次0.2220、0.4866、0.3640。控制深度和组成后中位rho降至0.06593、0.06343、0.11454，原P分别0.1671、0.3593、0.1671。保留未调整的FDR显著对应，同时明确不能称独立驱动或完整代谢因果链。

TJ-RCC调整后TYMP8/9位正向、中位rho0.1621，多胺6/9位正向、中位0.1136，仅描述。GLMF2多胺正向、胆固醇负向，同样保留。

单细胞“C4相对其他巨噬的胆固醇程序较低”与混合spot“C4评分与胆固醇程序同变”是不同层次的比较；后者受细胞组成及其他细胞来源影响，不能据此推翻前者或说明C4内胆固醇代谢增强。MAT2A/MAT2B/AMD1/SRM/SMS只是转录上下文，不是实测MTA或巨噬细胞内代谢通量。

## 辅助：文献过程对照及其解释

本轮只做有限过程对齐：APOE巨噬/SPP1份额与POSTN-CAF的对应及其连续交互。采用作者冻结程序去共享版本和本项目统一边界，不是作者全部病理、蛋白和动物实验的精确复现。

APOE×POSTN-CAF对总T的血管调整模型，GSE175540 19位中15位负交互，中位β=−0.01448，P=0.01921、FDR=0.15369，为名义显著。核心模型P=0.1671，伴随Wilcoxon P=0.04937/FDR=0.1975，完整保留。

**条件斜率核对限制了这项信号的解释。** 血管调整后，在APOE低/高状态时，POSTN-CAF与总T内外差的条件斜率中位数分别+0.08394/+0.08097，高APOE条件只有5/19位斜率为负。因此这是正向关系减弱的交互，不足以称“二者共同造成T细胞排斥”。这一解释比看到负β后直接称屏障增强更准确。

同队列APOE–POSTN源带相关核心调整rho0.02477/P0.6476；SPP1–POSTN为rho0.2153/P0.06357，未通过主原P标准。没有为了完成复现而把所有过程端点归为阳性。

## 主线组织与证据边界

研究重点已经转到“代谢相关巨噬状态参与并调节免疫屏障”，没有继续要求巨噬细胞单独构成屏障。现阶段证据可组织为：
- 已有上游代谢发现与客户实验背景，各自按原资料来源保留；
- 本项目C4的SPP1富集、APOE高表达和脂质处理/溶酶体状态；
- 已确立的CD8内外空间分区；
- 状态与代谢背景的空间对应，以及合作效应的探索性检验。

前述各层不能直接拼成已经验证的“TYMP/MTA→C4→CAF屏障→CD8阻挡”因果箭头。本轮没有支持最后一段稳定调节关系；后续若以此为核心机制，应直接评价代谢条件改变是否改变巨噬对CAF/基质和T细胞进入的作用。此处是研究缺口定位，不启动或冒充新湿实验。

边界仍由表达分区定义，独立病理/多重定位缺口继承Phase276。C4评分是状态相关转录程序，不是独立细胞丰度；APOE程序不能被直接更名为C4。各研究、平台、病例重复和临床数据来源保持分开。

文献背景：肾癌APOE–POSTN屏障研究（Science Advances2026，https://pubmed.ncbi.nlm.nih.gov/42566529/）与本项目部分空间样本重叠，因此不算完全独立的空间队列。其干预与多重定位证据不能直接移用于本项目。

## 验收与交付

原方案1292条切片效应、1258条患者效应、102队列端点，加4个单细胞身份/功能端点；技术扩展836/814条切片/患者效应、66队列端点。两轮计算均exit0。

原方案720项模型/相关独立复算最大误差2.23e−15、24595项源带聚合核对和76项独立几何核对通过；技术扩展618项模型复算最大误差1.98e−15、11970项聚合核对通过。P值用独立组合公式、BH手工复算、交互系数与条件斜率用statsmodels验证。所有原始H5AD元数据、原派生输入SHA及冻结方案记录保持。

交付目录分original19与harmonized_platform17，避免混淆技术版本和推断家族；全部正负统计、模型资格、条件斜率、基因列表、脚本和验收随附。新图0，正式稿修改0。
'''
(D/'代谢巨噬状态与屏障协同检验报告.md').write_text(text,encoding='utf-8')
for src,label in [(P,'original19'),(H,'harmonized_platform17')]:
 target=D/label
 for sub in ['tables','code','qc']:(target/sub).mkdir(parents=True,exist_ok=True)
 for f in (src/'tables').glob('*.tsv'):shutil.copy2(f,target/'tables'/f.name)
 for f in (src/'code').glob('*.py'):
  if f.name!='assemble_delivery.py':shutil.copy2(f,target/'code'/f.name)
 for f in (src/'qc').glob('*.json'):
  if f.name not in ['LOCAL_DELIVERY_VERIFIED.json']:shutil.copy2(f,target/'qc'/f.name)
for f in ['PROTOCOL_BEFORE_TESTING.md','PROGRAMS_FROZEN.json','FREEZE.json']:shutil.copy2(P/f,D/f)
shutil.copy2(H/'ADDENDUM_BEFORE_HARMONIZED_MODELS.json',D/'harmonized_platform17/ADDENDUM_BEFORE_HARMONIZED_MODELS.json')
v=json.loads((P/'qc/FINAL_VALIDATION.json').read_text());hv=json.loads((H/'qc/FINAL_VALIDATION.json').read_text());assert v['status']==hv['status']=='PASS'
combined=dict(status='PASS',at=datetime.datetime.now().astimezone().isoformat(),original_validation=v,harmonized_validation=hv,primary_original_preserved=True,technical_extension_not_independent_validation=True,conditional_slope_interpretation_checked=True,figures=0,manuscript_changes=0)
(P/'qc/COMBINED_VALIDATION.json').write_text(json.dumps(combined,indent=2));(D/'COMBINED_VALIDATION.json').write_text(json.dumps(combined,indent=2))
files=sorted(f for f in D.rglob('*') if f.is_file() and f.name!='DELIVERY_SHA256.tsv')
pd.DataFrame([dict(file=f.relative_to(D).as_posix(),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files]).to_csv(D/'DELIVERY_SHA256.tsv',sep='\t',index=False)
with zipfile.ZipFile(P/'phase277_delivery.zip','w',compression=zipfile.ZIP_DEFLATED) as z:
 for f in sorted(D.rglob('*')):
  if f.is_file():z.write(f,arcname=f.relative_to(D).as_posix())
(P/'STATUS.json').write_text(json.dumps(dict(status='VALIDATED_AWAITING_DELIVERY',at=datetime.datetime.now().astimezone().isoformat(),files=len(files),figures=0,background_jobs=0),indent=2))
print('DELIVERY_FILES',len(files),'ZIP_BYTES',(P/'phase277_delivery.zip').stat().st_size)
