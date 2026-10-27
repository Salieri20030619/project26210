from pathlib import Path
import pandas as pd,numpy as np,json
p=Path(__file__).resolve().parents[1]
def read(n):return pd.read_csv(p/'tables'/n,sep='\t')
s=read('C4_SURVIVAL_MODELS.tsv');tv=read('C4_PFS_TIME_VARYING_MODELS.tsv');hr=read('C4_PFS_TIME_SPECIFIC_HR.tsv');se=read('C4_TIME_VARIATION_TRIAL_SENSITIVITY.tsv');loo=read('C4_TIME_VARIATION_LEAVE_ONE_PATIENT.tsv')
def one(df,**kw):
 for k,v in kw.items():df=df[df[k]==v]
 assert len(df)==1
 return df.iloc[0]
def fmt(x):return f'{x:.4g}'
def ci(q):return f'{q.HR:.3f}（{q.CI_low:.3f}–{q.CI_high:.3f}）'
q=one(tv,model_id=12);h6=one(hr,model_id=12,month=6);h12=one(hr,model_id=12,month=12)
st=one(se,parent_model_id=12,sensitivity='stratify_trial');cm=one(se,parent_model_id=12,sensitivity='CM025_only')
def md(df):
 return '| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(map(str,row))+' |' for row in df.itertuples(index=False,name=None))
t=[]
for i,name in [(1,'Nivolumab OS：临床调整'),(4,'Nivolumab OS：再调文献评分/总巨噬'),(9,'Nivolumab PFS：临床调整'),(12,'Nivolumab PFS：再调文献评分/总巨噬'),(20,'JAVELIN联合组PFS：完整调整'),(24,'Sunitinib PFS：完整调整')]:
 z=one(s,model_id=i);t.append({'模型':name,'n/事件':f'{int(z.n)}/{int(z.events)}','HR（95%CI）':ci(z),'P':fmt(z.P),'FDR':fmt(z.FDR),'PH检验P':fmt(z.PH_P)})
ht=[]
for tm in [3,6,12,24]:
 z=one(hr,model_id=12,month=tm);ht.append({'月数':tm,'HR（95%CI）':ci(z),'在险人数':int(z.at_risk)})
genes=json.loads((p/'qc/CLINICAL_SCORE_FREEZE.json').read_text())['genes_used']
report=f"""# 26210：C4状态临床延伸分析
Phase271 v1，2026-09-28。内部完整分析与写作建议；主稿仍为Phase267。

## 结论
本轮已从“复现SPP1文献”推进到项目自身C4状态的临床分析。新建的C4状态评分具有跨研究可识别性，并在CheckMate nivolumab组显示独立于文献五基因评分和总巨噬估计的时间依赖性PFS关联。它可以支撑一个探索性结果小节，主体是本项目定义的状态及其临床异质性。

当前数据不支持把C4直接写成统一的免疫耐药指标：独立OS和治疗获益交互未获支持，JAVELIN联合组也没有相同的恒定PFS关联。不能把已复现的SPP1文献风险直接移植到C4。

## 1. 状态定义和验证
保持原C4的2674个细胞成员，不重新聚类。基因发现使用肿瘤组织中24位独立患者、5个研究的配对C4/其他成熟巨噬伪bulk，每组至少20细胞。

事前登记筛选：配对Wilcoxon原始P<0.05；中位log2CPM差异≥0.5；至少70%患者同向；至少80%研究的中位差异为正；C4中位检测比例≥0.1。按研究平衡效应及患者效应排序固定前20基因。完整原始P/FDR保留，不以FDR筛除名义信号。

事先排除SPP1、TYMP、FN1、CLEC5A、OLR1、IL1RN及线粒体/核糖体基因；临床结局未用于基因选择。CheckMate未覆盖CCL18，因此两临床队列统一使用剩余19基因，不补选：
{', '.join(genes)}。

留一研究时每折只用训练研究选基因。24/24位留出患者的C4分数高于同患者其他巨噬，五研究分别4/4、1/1、2/2、12/12、5/5。去除CCL18后仍24/24，中位差异2.376 log2CPM单位。留出研究的跨细胞类型比较中，C4的评分中位数在五研究均最高，其他成熟巨噬/单核次之；成纤维参考覆盖不足等既有边界保留。

这是基于既有冻结标签的研究隔离内部验证，没有每折重新聚类，也不是新增外部单细胞队列。评分表示C4相关转录状态，不能直接视作C4细胞比例；旧NNLS混合恢复未通过的问题仍保留。

## 2. 与文献及总巨噬指标的区分
C4与SPP1单基因的Spearman rho为0.123/0.135（CheckMate/JAVELIN）；与文献五基因为0.426/0.476；与TYMP为0.155/0.189；与官方总巨噬估计为0.073/0.265。均为描述性比较，没有用相关P选择指标。这些结果支持指标之间有区别，但不能据此确认TYMP—C4因果链。

## 3. 恒定HR模型与完整结果
CheckMate调整年龄、统一性别、试验和MSKCC；JAVELIN调整年龄、性别。分别追加SPP1、文献五基因，以及文献评分+总巨噬估计；在每一模型的同一批完整病例上做嵌套增量检验。评分HR按治疗臂内1 SD。24项模型全部保留，以下为关键行。

{md(pd.DataFrame(t))}

C4独立OS关联当前未获支持。Sunitinib完整调整模型呈较低进展风险关联，必须保留；未调整文献评分时P=0.1267，因此不能将这个条件关联写成无条件稳健的TKI疗效预测。

Nivolumab PFS的比例风险假设明显不满足，单个恒定HR不足以描述关联。后续时间模型由这一诊断触发，基因与评分保持不变；不是扫描截点寻找阳性。

## 4. 时间依赖性PFS关联
诊断后另行登记β(t)=β0+β1×log(1+t)，t为治疗后月数。时间交互在各事件风险集中计算。基线RNA评分固定，随时间变化的是风险关联系数，不能推断C4在体内发生功能逆转。

完整调整nivolumab模型n=164，145个事件。加入C4和C4×log(1+t)的两自由度增量似然比为{q.joint_LR:.3f}，P={fmt(q.joint_P)}，FDR={fmt(q.joint_FDR)}。时间变化项相对恒定C4模型P={fmt(q.time_dependence_P)}，FDR={fmt(q.time_dependence_FDR)}。

{md(pd.DataFrame(ht))}

时点区间是同一模型的点估计，不是多次独立验证。t=0截距也完整保存，但不拿它证明“早期耐药”。晚期在险人数减少且风险集构成改变，较低条件风险不等于干预C4可延长生存。

仅临床调整模型也有时间关联（联合P={fmt(one(tv,model_id=9).joint_P)}）。Everolimus对应两种模型的联合P=0.1666/0.1401；不能以两臂显著性不同宣称免疫治疗特异性。

## 5. 稳健性与治疗交互
允许各试验具有不同基线风险后，完整调整模型联合P={fmt(st.joint_P)}，FDR={fmt(st.joint_FDR)}。仅CM025 nivolumab患者n=119、104事件，联合P={fmt(cm.joint_P)}，FDR={fmt(cm.joint_FDR)}。

逐一去掉164位患者，时间系数均保持负向；6个月HR范围{loo.HR_month6.min():.3f}–{loo.HR_month6.max():.3f}，12个月HR范围{loo.HR_month12.min():.3f}–{loo.HR_month12.max():.3f}。正式估计未排除任何患者。这些是同队列的敏感性分析，不冒充独立验证。

CM025完整调整的C4恒定OS交互P=0.7479；时变PFS的两自由度治疗交互P=0.4180、FDR=0.5886；JAVELIN完整调整PFS交互P=0.1427。因此当前数据不支持C4预测相对治疗获益。

## 6. 统计口径与实现修正
基因筛选和临床模型分别登记；时间函数及试验敏感性在各自计算前登记。临床恒定模型BH按队列×终点×调整方式的两治疗臂家族；时变主模型按两治疗臂，时变治疗交互按两调整模型，试验敏感性按两种固定检查分别计算。原始P、FDR及家族均保留。时变分析由诊断触发，整体为回顾性探索；FDR显著不改变其研究设计性质。

独立检查发现CheckMate原合并表同时含F/M与Female/Male，性别四水平与试验存在冗余。已统一为两个类别后重算受影响模型，并增加设计矩阵满秩及独立P有限性检查。初版完整留在attempt_raw_sex_labels，不用于正式结论。

Phase270受影响的8项调整后比较也单列修正：文献五基因nivolumab OS现为HR=1.270（1.045–1.543）、P=0.01611、FDR=0.06444，仍属名义显著，替代此前调整后HR1.256/P0.02080。未调整文献复现、JAVELIN和单独CM025交互不受标签混用影响；旧文件未覆盖。

24项恒定模型、6项交互及其增量检验用Python独立复算；4项时变模型和2项时变交互用显式事件风险集展开再次核对，与R的tt实现一致。模型、输入患者、完整系数、PH检验、19基因覆盖及训练折基因均保留。4项试验敏感性和164次逐患者遗漏另表。没有把细胞数或重复拟合次数作为独立样本量。

方法参考：[survival官方时间依赖系数说明](https://stat.ethz.ch/R-manual/R-devel/RHOME/library/survival/doc/timedep.pdf)。未将患者最终随访时间直接乘入基线评分。

## 7. 文章定位
建议标题：**C4相关脂质处理状态与进展风险的时间依赖性关联**。

正文主线为“项目定义状态—跨研究可识别—与SPP1负荷区分—调整后临床关联及稳健性”。文献复现作为参照，不占据新增章节主体。OS阴性、治疗交互阴性、跨方案差异及回顾性边界在Discussion集中说明；Sunitinib条件关联保留，不把反向结果隐藏。

尚未建立独立临床队列复现、预测准确性或临床决策效用。不能称C4耐药预测器，也不能把TYMP、候选代谢物、C4与临床结局写成已证实的完整因果通路。本轮未改稿、未绘图；另附可整合的结果段。
"""
(p/'delivery/26210_C4状态临床延伸分析报告.md').write_text(report,encoding='utf-8')
draft=f"""# 拟新增结果小节与图表安排
供统一稿整合；编号继承现有顺序。本轮未写入Phase267。

## 结果标题
C4相关脂质处理状态与进展风险的时间依赖性关联

基于24位患者的肿瘤巨噬细胞配对表达谱，我们建立了不含SPP1、TYMP及文献五基因成员的C4状态评分。留一研究内部验证中，全部24位留出患者的C4评分均高于同患者其他成熟巨噬细胞，五研究方向一致。临床分析采用两队列共同覆盖的19基因。该评分与SPP1单基因的相关性较弱（CheckMate和JAVELIN分别ρ=0.123和0.135），与文献五基因评分呈中等相关（ρ=0.426和0.476），支持其反映不同的巨噬转录状态信息。

在接受nivolumab治疗的患者中，C4评分与PFS的关联随随访时间变化。调整年龄、性别、试验、MSKCC风险、文献五基因评分及总巨噬估计后，C4时间模型提供了增量关联信息（n=164，145个事件；两自由度似然比={q.joint_LR:.2f}，P={fmt(q.joint_P)}，FDR={fmt(q.joint_FDR)}）。每升高1 SD评分的条件进展风险HR在6个月为{ci(h6)}，12个月为{ci(h12)}。按试验分层基线风险及仅保留CM025患者时，关联仍获支持（P分别为{fmt(st.joint_P)}和{fmt(cm.joint_P)}），逐患者遗漏后估计方向保持一致。

在JAVELIN的sunitinib组中，调整文献评分和总巨噬估计后，较高C4评分与较低进展风险相关（HR=0.799，95%CI 0.684–0.934，P=0.00491，FDR=0.00983）。这些结果显示C4状态与SPP1总体转录负荷具有不同的临床关联特征。

## Methods必须补充
报告19基因完整名单、24患者配对发现、仅训练研究选基因的留一研究验证、临床共同覆盖和评分公式。CheckMate采用已标准化发布值均值，JAVELIN为log2(TPM+1)均值；统一性别标签；明确事件编码与完整病例标准。报告Cox协变量、PH诊断、诊断后固定β(t)=β0+β1log(1+t)、两自由度增量检验及各BH家族。详细参数见本轮完整报告和冻结计划。

## Discussion集中说明
C4尚未建立独立OS关联或相对治疗获益预测。JAVELIN联合组未显示相同恒定PFS关联；sunitinib较低风险为控制文献评分后的条件关联。时间变化是基线评分与条件风险的关系变化，不证明C4功能逆转；12/24个月在险人数仅27/17。该扩展由PH诊断触发，属回顾性探索；试验及患者敏感性是内部稳健性，尚无独立临床队列验证。评分不能直接等于C4比例，也不能闭合TYMP—代谢物—C4的因果链。

## 建议面板（尚未绘制）
1. 24患者、5研究的留出C4/其他巨噬评分差异。
2. C4与SPP1单基因、文献五基因及总巨噬估计的相关性。
3. 完整调整后的C4-PFS时间HR曲线，95%CI及固定时点在险人数。
4. 分试验基线风险及CM025敏感性，配合完整统计表。
文献OS复现作为对照或补充；OS、治疗交互及所有队列完整结果放补充表，关键限制在Discussion明确。
"""
(p/'delivery/拟新增结果小节与图表安排.md').write_text(draft,encoding='utf-8')
print('REPORTS_WRITTEN',len(report),len(draft))
