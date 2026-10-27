from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,shutil,datetime
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=P/'phase162_MTA_polyamine_part1_reassessment_20260915_v1';Q=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1')
D=pd.read_csv(O/'PART1_GENE_CONTRAST_SUMMARY.tsv',sep='\t');C=pd.read_csv(O/'PART1_TYMP_GENE_CORRELATION_SUMMARY.tsv',sep='\t');M=pd.read_csv(O/'PUBLIC_META_AXIS_ASSOCIATIONS.tsv',sep='\t');R=pd.read_csv(O/'PUBLIC_COHORT_AXIS_ASSOCIATIONS.tsv',sep='\t');TC=pd.read_csv(O/'TCGA_MTAP_SHALLOW_DEEP_DISTINCTION.tsv',sep='\t');I=pd.read_csv(O/'PUBLIC_META_MTAP_INTERACTION.tsv',sep='\t');EXP=pd.read_csv(O/'PART1_PATIENT_TARGET_EXPRESSION.tsv',sep='\t')
def md(df):
 out=['| '+' | '.join(df.columns)+' |','|'+'|'.join(['---']*len(df.columns))+'|']
 for row in df.itertuples(index=False,name=None):out.append('| '+' | '.join(f'{x:.6g}' if isinstance(x,(float,np.floating)) else str(x) for x in row)+' |')
 return '\n'.join(out)
module_file=Path('/mnt/BioDisk/ludingyu/project26210/envs/part2_6_2_7_20260726_v1/sources/scFEA/data/module_gene_m168.csv')
mg=pd.read_csv(module_file,index_col=0);cross=[]
for gene in ['ODC1','SRM','SMS','MTAP','AMD1','MAT2A','TYMP']:
 cross.append(dict(gene=gene,scFEA_modules=','.join(mg.index[(mg.astype(str)==gene).any(axis=1)].astype(str))))
pd.DataFrame(cross).to_csv(O/'FIRST_PART_scFEA_GENE_CROSSWALK.tsv',sep='\t',index=False)
summary="""# MTA–多胺轴可与第一部分对上，应区分合成端和MTAP清除背景

本轮复核第一部分冻结11950恶性细胞、原14患者配对、四公共配对代谢组队列255患者及TCGA-KIRC 507患者。MTA不放弃，SLC16A3上游追加也保留。

**与第一部分最清楚的对应是SRM及多胺合成背景。**
- 原半胱氨酸/甲硫氨酸AUCell通路在TYMP阳性恶性细胞中14/14患者升高，P=0.00110、FDR=0.00226；scFEA相关轴也14/14同向。但模型没有显式MTA，不能把该分数叫实测MTA或MTA净生成通量。
- 新拆解SRM：TYMP阳性相对阴性恶性细胞11/14患者升高，中位log2CPM差+0.290、P=0.00671、FDR=0.0403；F1高克隆中11/14升高，+0.518、P=0.0134、FDR=0.0806。
- ODC1在TYMP阳性组10/14升高，+0.228、P=0.0419、FDR=0.1256，为名义支持；SMS、MTAP没有一致升高/降低对应，完整保留。

**公共数据把SRM与SMS的角色区分开。**
- TYMP—SRM四队列原始均正，随机合并rho=0.386、P=0.000132、FDR=0.000791；控制髓系RNA代理后rho=0.218、P=0.0343、FDR=0.2056。
- SMS—实测MTA四队列均正；控制TYMP、MTAP和髓系代理后rho=0.250、P=0.000578、FDR=0.00289。SMS—组织SPP1的相应TYMP/髓系调整关联也为正，rho=0.354、P=7.76e-7、FDR=4.66e-6。
- 这些结果支持合成轴的分段联系，不等于已确认TYMP—单个基因—MTA—SPP1完整因果中介。

**MTAP应作为清除能力/遗传背景因素保留，不能预设普遍失活。**
- 第一部分TYMP阳性组MTAP没有一致下降（8/14升高，中位差+0.077，P=0.463）；F1高克隆有9/14降低的弱趋势（差-0.402，P=0.194）。
- TCGA历史148“丢失”中，142为浅层-1、仅6为深度-2；深度调用占507患者的1.18%。不能把148都称为完全缺失。
- 公共MTAP—MTA原始四队列均负，rho=-0.250、P=0.00957、FDR=0.0287，支持清除端联系；但低MTAP放大TYMP—MTA关系没有得到支持，连续交互项反为负（髓系调整beta=-0.194、P=0.0335）。

另有直接相关的公开786-O研究：作者在有MTAP表达的亲本786-O上用CRISPR构建MTAP敲除株，因此不能默认786-O天然MTAP缺失；客户这批细胞的实际MTAP状态没有由现有图证明。来源：https://www.nature.com/articles/s41392-019-0035-z

客户图确认shTYMP后MTA下降；44倍属于用户提及的实验幅度，现有图没有被独立核算为该比值，不能归因为MTAP缺陷。无需客户原始数据即可继续公共调查。

**SLC16A3上游补充：** 第一部分完整恶性细胞中，16合格患者深度调整rho中位0.071，13/16正，P=0.0443、追加五项FDR=0.0580；再控制周期P=0.0833。与此前人群重叠，属于追加方向性支持，不是独立验证。

建议保持MTA–多胺合成轴与客户实验的对应主线，MTAP单独作为清除背景，SLC16A3并行补充；整体巨噬承载SPP1联系，C4保留胆固醇处理状态定位。702项核查通过，无新图、无客户数据入模。
"""
(O/'SUMMARY_CN.md').write_text(summary,encoding='utf-8')
report="""# MTA–多胺轴与第一部分的重新对照

## 结论：多胺合成端可与第一部分对应，MTAP定位为清除与遗传背景因素

"""+summary.split('\n\n',1)[1]+"""

## 生物学环节与本次检验分工

ODC1为多胺合成提供腐胺；SRM和SMS利用脱羧SAM进行氨丙基转移，并产生MTA。MAT2A和AMD1分别对应SAM供给及脱羧SAM供给。MTAP参与MTA的清除和甲硫氨酸补救循环。因此，合成相关基因与MTAP不应相加后笼统解释为“MTA更高”。[SRM反应](https://www.reactome.org/content/detail/R-HSA-351215)、[SMS反应](https://reactome.org/content/detail/R-HSA-351210)、[MTAP反应](https://www.reactome.org/content/detail/R-HSA-1237160)

本轮分别检验合成基因、SAM供给和MTAP清除，不把MTAP表达低、基因深度缺失、蛋白缺失、酶活降低当成同一概念。

## Methods

第一部分采用冻结11950恶性细胞原始计数，来自26患者；使用原有TYMP阳性/阴性及F1高克隆/其他克隆标签，并复用原通路比较的14患者配对清单。基因伪bulk表达为log2(CPM+1)。重算TYMP逐细胞标准化表达与原第一部分数值最大绝对差5.33e-15，确认口径一致。

轴基因为ODC1、SRM、SMS、MTAP、AMD1、MAT2A。每种原有对比分别进行患者配对Wilcoxon，六基因构成同类家族BH校正。逐患者细胞共表达使用log1p(CP10K)秩，报告原始、深度及深度+G2-M调整；患者至少50细胞，TYMP及目标基因各至少10检出细胞。各患者rho作为描述性效应，以患者为独立单位进行符号秩检验。细胞周期分数排除全部目标基因。SLC16A3作为独立上游追加目标，两个状态比较及三个连续模型额外提供合计五项的保守补充FDR；该追加家族不替代之前的发现家族。

公共数据沿用RC18、RC20、CPTAC、CPTAC_val的255名已配对肿瘤患者及实测MTA；新增读取六基因的log2(TPM+1)。ODC1/SRM/SMS均值z分数为生产端RNA模块，AMD1/MAT2A均值z分数为SAM供给RNA模块。分别计算TYMP—X、X—MTA/SPP1及控制TYMP的条件关联，报告原始与髓系RNA代理调整；生产端另控制MTAP。六基因与两个模块分别按终点、边和调整口径校正。合并采用Fisher z、有效协变量数修正的近似方差及DerSimonian–Laird随机效应，保留单队列及I²。

额外检验MTA标准化秩与TYMP标准化秩、负MTAP标准化秩及其乘积的关系，使用HC3稳健标准误，并报告有/无髓系代理调整。该交互为预列的清除背景修饰检验，不作二分阈值搜索。四队列交互系数按逆方差随机效应合并。

遗传证据使用已有TCGA-KIRC 507患者MTAP GISTIC调用与RNA结果，将深度缺失-2、浅层丢失-1、中性0及增益1分开。旧MTAP-loss混合分组、旧SPP1程序交互结果保留。未将单细胞RNA零表达或全局CNV负荷当作MTAP基因缺失。

## Results

### 第一部分支持甲硫氨酸与多胺相关代谢状态

原第一部分TYMP阳性恶性细胞的半胱氨酸/甲硫氨酸AUCell通路14/14升高，P=0.00110、FDR=0.00226；精氨酸/脯氨酸通路13/14升高，P=0.00137、FDR=0.00226。F1高克隆的这两条通路均12/14升高。scFEA相关轴也有对应方向。

模型映射显示ODC1进入M62，SRM/SMS进入M31/M69/M70，MTAP进入M31，AMD1进入M31/M70，MAT2A进入M30/M31/M70。这说明第一部分模型包含轴内基因，却将生成与补救循环放在相邻或共同模块中；其原范围表明确MTA不是显式代谢物，属于间接Cys/Met与Arg/Pro代理。该结果支持代谢网络背景，而非直接测量MTA积累。

### SRM是与第一部分TYMP状态最清楚对应的合成基因

SRM在TYMP阳性相对阴性恶性细胞中11/14患者升高，中位log2CPM差+0.290，P=0.00671、FDR=0.0403；在F1高克隆中11/14升高，+0.518，P=0.0134、FDR=0.0806，后者为名义支持。患者内TYMP—SRM在深度调整后11/14正，中位rho=0.0332、P=0.0203、FDR=0.0608；追加细胞周期后仍11/14正，中位rho=0.0311、P=0.0295、FDR=0.1006。

ODC1在TYMP阳性组10/14升高，中位差+0.228、P=0.0419、FDR=0.1256；但细胞内深度调整关联接近零，保留为状态层名义线索。SMS在TYMP阳性组6/14升高、8/14降低，差-0.0766、P=0.326，未显示同步上调。不能写成ODC1/SRM/SMS全部随TYMP升高。

AMD1细胞内TYMP关联在深度调整后12/16正、中位rho=0.0316、P=0.0110、FDR=0.0608；MAT2A较弱。SAM供给背景有部分对应，不能替代实测SAM或脱羧SAM通量。

### 公共组织数据将TYMP—SRM与SMS—MTA两段区分开

TYMP—SRM四队列原始均正，随机合并rho=0.386、P=0.000132、FDR=0.000791；控制髓系代理后仍四队列正，rho=0.218、P=0.0343、FDR=0.2056。该组织层对应与第一部分SRM状态结果同向，调整后证据为名义支持。

SMS—实测MTA原始四队列均正，合并rho=0.278、P=0.00105、FDR=0.00631；控制TYMP和MTAP后rho=0.232、P=0.000295、FDR=0.00148；再控制髓系代理后rho=0.250、P=0.000578、FDR=0.00289，仍四队列正。SMS对应MTA的证据较SRM更强，但TYMP—SMS并没有稳定强上游，原始合并rho=0.258、P=0.0815、FDR=0.1629。

生产端RNA模块与TYMP原始四队列均正，rho=0.284、P=0.000256、FDR=0.000511；控制TYMP及髓系代理后，模块—MTA rho=0.190、P=0.00323、FDR=0.00645；再加MTAP后降至rho=0.139、P=0.0865、FDR=0.1273。模块反映表达协同，不称MTA净生成量。

在SPP1终点，控制TYMP和髓系代理后，SMS—组织SPP1四队列均正，rho=0.354、P=7.76e-7、FDR=4.66e-6；SRM对应rho=0.210、P=0.00168、FDR=0.00505；生产端模块rho=0.328、P=1.82e-7、FDR=3.65e-7。这些结果支持多胺相关表达与SPP1联系，但不将不同基因的最强边拼接为已验证串联中介。

### MTAP清除端有公共支持，普遍失活与TYMP放大机制未获支持

第一部分TYMP阳性组MTAP有8/14升高，中位差+0.0773、P=0.463；没有一致降低。F1高克隆有9/14降低，差-0.402、P=0.194、FDR=0.311，保留为弱方向线索。细胞数至少50的18患者恶性伪bulk均存在MTAP读数，单细胞检出率差异较大；这不能确诊或排除蛋白/酶功能缺陷。

TCGA-KIRC 507患者中，深度缺失调用6例（1.18%），浅层丢失142例（28.0%），中性348例，增益11例。深度缺失组MTAP表达z中位-3.217，中性组+0.282。旧148个“loss”是142浅层与6深度合并，不能将全部148称作纯合缺失；也没有与第一部分每个患者直接配对的DNA或蛋白证据。

公共MTAP—MTA原始四队列均负，随机合并rho=-0.250、P=0.00957、FDR=0.0287；控制TYMP后rho=-0.239、P=0.0143、FDR=0.0430。控制髓系代理后边际rho=-0.216、P=0.0601、FDR=0.1202，异质性较高。它支持将MTAP保留为清除端变量，但不自动说明它是TYMP作用的中介，因为TYMP—MTAP原始合并rho=0.00507、P=0.937。

低MTAP是否增强TYMP—MTA关联的连续交互，原始合并beta=-0.212、P=0.0130；髓系调整beta=-0.194，95%区间[-0.373,-0.0152]、P=0.0335，四队列只有一个正向。当前数据不支持预期的正向放大，反向交互必须保留。该结果属于低RNA表达的组织层修饰检验，不等价于深度DNA缺失的功能实验。

旧TCGA混合loss对旧SPP1程序的交互beta=0.067、P=0.255、FDR=0.382保留；旧程序不等于当前C4，也不是实测MTA，不据该旧阴性结论排除MTA—整体巨噬SPP1路线。

### 公开786-O模型不支持默认天然MTAP缺失

与客户细胞系直接相关的原始研究检测了MTAP蛋白，并在有MTAP表达的亲本786-O上通过CRISPR构建MTAP敲除株。该结果说明786-O可以是MTAP表达完整的背景，不能因文献讨论MTAP缺陷肾癌就把亲本786-O标成天然缺失。客户这批细胞的具体基因/蛋白状态仍未由其提供图证明，本轮不索要原始数据或自行赋予缺失标签。[786-O MTAP原始研究](https://www.nature.com/articles/s41392-019-0035-z)

### 第一部分追加了SLC16A3上游状态线索，幅度仍小

TYMP阳性组SLC16A3在11/14患者升高，差+0.262、P=0.0295、追加五项FDR=0.0580；F1高克隆中11正、2负、1零，差+0.301、P=0.0464、FDR=0.0580。完整第一部分中16合格患者的深度调整TYMP—SLC16A3为13正、3负，中位rho=0.0705、P=0.0443、FDR=0.0580；深度+周期调整中位rho=0.0573、P=0.0833、FDR=0.0833。

这些是原第一部分人群扩展后的追加探索，与Phase161的同文库配对人群重叠。不能因为P从0.083变为0.044就称独立验证或强上游已确立；方向性支持保留，MTA轴与SLC16A3继续并行。

## Discussion与推荐叙事

第一部分的TYMP相关甲硫氨酸/多胺代谢状态、SRM上调，公共组织中的SMS—MTA和多胺相关SPP1联系，加上客户TYMP干预后MTA下降，构成可进一步衔接的分层证据。更适合当前数据的工作假说是“TYMP相关代谢状态与多胺合成端协同，MTA处于生成与清除共同决定的平衡中；MTAP作为清除与遗传背景变量”。不要将它改写成“TYMP必然造成MTAP失活，导致MTA增加44倍”。

客户提供图确认shTYMP后MTA下降、Thymidine上升，并伴随Methionine/SAM/Spermidine等下降。44倍为用户提及的实验幅度，本轮没有从图中倒算或独立确认具体比值；公共数据不需要复刻客户干预倍数，也不能用RNA较低解释精确倍数。无需客户原始数据即可开展当前公共分析。

MTAP在不同癌种的组织与培养模型中可能表现不同。其他癌种原始研究显示，间质中的MTAP清除可影响组织MTA水平，因此基因缺失不必与整体组织MTA幅度一一对应；这提供背景解释，不能直接替代本项目肾癌证据。[MTAP缺失与组织MTA原始研究](https://www.nature.com/articles/s41467-021-24240-3)

本轮所有表达、代谢相关及F1状态比较为观察性证据；原第一部分F1由同一表达数据发现，其功能对应不是外部独立复现。公共组织SPP1不是巨噬特异表达，当前C4也没有同患者实测MTA。保留Phase157中MTA转录代理在C4较低的结果，但代理不等于实测MTA。用户确认的双层分工不变：整体巨噬承载SPP1联系，C4用于SPP1信号集中和胆固醇处理状态定位，不要求每个MTA相关环节必须先通过C4比例检验。

后续优先核查合成端SRM/SMS与客户MTA变化的对应；MTAP深度缺失亚组及功能缺陷需要单独定位，不能以浅层拷贝减少或单细胞零表达代替。SLC16A3的上游追加继续保留，尚无确认的完整因果中介。

## 验证和完整统计

原始第一部分细胞及TYMP计数一致；标准化表达与旧定义一致；14患者比较名单冻结复用。公共偏相关采用独立精度矩阵公式核对，合并效应与FDR独立复算。702/702检查通过，无新图、无客户实验数据进入统计。单细胞、公共代谢组和TCGA是不同证据层，不将样本量相加为一个实验。

"""
for title,df in [('第一部分基因状态比较',D),('第一部分患者内相关',C),('TCGA MTAP拷贝状态分层',TC),('公共MTA轴全部合并结果',M[M.outcome=='MET:5-methylthioadenosine (MTA)']),('公共SPP1全部合并结果',M[M.outcome=='SPP1']),('MTAP修饰交互',I),('第一部分MTAP患者表达',EXP[EXP.gene=='MTAP'])]:
 report+='\n\n### '+title+'\n\n'+md(df)
report+='\n\n完整单队列效应与患者明细见PUBLIC_COHORT_AXIS_ASSOCIATIONS.tsv、PART1_FROZEN_PATIENT_GENE_CONTRASTS.tsv及覆盖表；未仅保留支持方向。\n'
(O/'MTA_POLYAMINE_PART1_REASSESSMENT_CN.md').write_text(report,encoding='utf-8')
status=json.loads((O/'STATUS.json').read_text());assert status['total_checks']==702
status.update(status='COMPLETE',report=str(O/'MTA_POLYAMINE_PART1_REASSESSMENT_CN.md'),MTA_retained=True,MTAP_role='clearance/genetic-background modifier; not established universal loss mediator',SLC16A3_upstream_supplement=True,confirmed_mediator=False,completion_time=datetime.datetime.now(datetime.timezone.utc).isoformat())
(O/'STATUS.json').write_text(json.dumps(status,indent=2))
sources=json.loads((O/'SOURCE_MANIFEST.json').read_text())
extra=[module_file,Path('/mnt/BioDisk/wangmingyuan/project26210_phase79_mtap_natural_perturbation_v1/results_v1/tcga_kirc_natural_perturbation_cohort_v1.tsv')]
extra+=list((P/'phase160_external_bridge_investigation_20260915_v1').glob('PATIENT_VALUES_*.tsv'))
sources +=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in extra];(O/'SOURCE_MANIFEST.json').write_text(json.dumps(sources,indent=2))
note="""## 2026-09-15 Phase162 MTA多胺轴与第一部分重审完成

用户要求继续追加上游、MTA不放弃，并重审ODC1/SRM/SMS/MTAP与第一部分。完成第一部分冻结11950恶性细胞（26患者，原14患者配对）、公共MTA队列255患者、TCGA507患者的分层核查。报告phase162_MTA_polyamine_part1_reassessment_20260915_v1/MTA_POLYAMINE_PART1_REASSESSMENT_CN.md；702项核查PASS，无图无客户数据入模。
第一部分CysMet AUC在TYMP阳性14/14升高；scFEA包含SRM/SMS/MTAP等模块但MTA不显式，不能称实测MTA/净生成。新SRM在TYMP阳性11/14上升、差0.290/P0.00671/FDR0.0403；F1高克隆11/14升、P0.0134/FDR0.0806。ODC1名义上升，SMS/MTAP未一致对应。
公共TYMP—SRM原始四队列正，rho0.386/P0.000132/FDR0.000791，髓系调整后名义rho0.218/P0.0343/FDR0.206。SMS—MTA四队列正，控制TYMP/MTAP/髓系后rho0.250/P0.000578/FDR0.00289；SMS—组织SPP1控制TYMP/髓系rho0.354/P7.76e-7/FDR4.66e-6。分段支持合成轴，非确认串联中介。
MTAP作为清除/遗传背景保留：TCGA148 loss=142浅层+6深度，深度仅6/507=1.18%；第一部分MTAP无普遍降低或与这些细胞配对的基因组缺失证明。公共MTAP—MTA原始四队列负rho-0.250/P0.00957/FDR0.0287；但低MTAP放大TYMP—MTA交互不支持，髓系调整合并beta-0.194/P0.0335反向。不得写普遍MTAP失活导致44倍。客户图shTYMP后MTA下降已确认，44倍未经独立核算，不索要原始数据也不以倍数筛公共结果。
SLC16A3第一部分扩展上游：深度调整16患者13正，rho中位0.0705/P0.0443/追加五项FDR0.058；周期调整P0.083，属重叠人群追加，不是独立验证。
公开786-O原始研究s41392-019-0035-z在MTAP表达亲本上另外构建KO，不能默认客户786-O天然MTAP缺失；该文献只作模型背景，不代表客户批次鉴定。
推荐MTA合成端SRM/SMS与客户背景主线，MTAP单列清除因素，SLC16A3并行；整体巨噬SPP1联系+C4胆固醇状态定位不变，暂无确认中介。

"""
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 p=P/name;shutil.copy2(p,O/(name+'.before'));p.write_text(note+p.read_text(encoding='utf-8'),encoding='utf-8')
with (P/'RESULTS_MANIFEST.tsv').open('a',encoding='utf-8') as f:f.write('\t'.join(['phase162_MTA_polyamine_part1_reassessment_20260915_v1','Phase162','report',str(O/'MTA_POLYAMINE_PART1_REASSESSMENT_CN.md'),status['completion_time'],'PASS','702 computational checks','Part1 frozen14-patient contrasts; MTA cohorts255; TCGA507; no figures'])+'\n')
p=P/'part2_rebuild_v1/STATUS.json';shutil.copy2(p,O/'PART2_STATUS.before.json');b=json.loads(p.read_text());b.update(latest_phase='Phase162',latest_report=str(O/'MTA_POLYAMINE_PART1_REASSESSMENT_CN.md'),stage='MTA_polyamine_part1_reassessment_complete',MTA_retained=True,C4_retained=True,confirmed_mediator=False);p.write_text(json.dumps(b,ensure_ascii=False,indent=2))
manifest={x:hashlib.sha256((O/x).read_bytes()).hexdigest() for x in ['SUMMARY_CN.md','MTA_POLYAMINE_PART1_REASSESSMENT_CN.md']};(O/'DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(manifest));print('COMPLETE')

(O/'normalization_and_supplement_fdr.py').write_text("from pathlib import Path\nimport pandas as pd,numpy as np,json\nfrom statsmodels.stats.multitest import multipletests\nP=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=P/'phase162_MTA_polyamine_part1_reassessment_20260915_v1'\nC=pd.read_csv(O/'PART1_TARGET_CELL_METADATA.tsv.gz',sep='\\t')\ndiff=np.abs(C['expr:TYMP']-C.TYMP_log1p_CP10K_v1)\nprint('ORIGINAL_TYMP_NORMALIZATION maxabs',float(diff.max()),'median',float(diff.median()))\nchecks=[dict(check='Part1_original_normalization',passed=bool(diff.max()<1e-6))]\nD=pd.read_csv(O/'PART1_GENE_CONTRAST_SUMMARY.tsv',sep='\\t');R=pd.read_csv(O/'PART1_TYMP_GENE_CORRELATION_SUMMARY.tsv',sep='\\t')\ndi=D.index[D.gene=='SLC16A3'];ri=R.index[R.gene=='SLC16A3']\nq=multipletests(np.r_[D.loc[di,'P'],R.loc[ri,'P']],method='fdr_bh')[1]\nD.loc[di,'FDR_SLC16A3_supplement5']=q[:len(di)];R.loc[ri,'FDR_SLC16A3_supplement5']=q[len(di):]\nD.to_csv(O/'PART1_GENE_CONTRAST_SUMMARY.tsv',sep='\\t',index=False);R.to_csv(O/'PART1_TYMP_GENE_CORRELATION_SUMMARY.tsv',sep='\\t',index=False)\npd.DataFrame(checks).to_csv(O/'NORMALIZATION_CHECK.tsv',sep='\\t',index=False)\nassert all(r['passed'] for r in checks)\ns=json.loads((O/'STATUS.json').read_text());s['total_checks']+=1;s['normalization_checks']=1;(O/'STATUS.json').write_text(json.dumps(s,indent=2))\nprint(D[D.gene=='SLC16A3'].to_string(index=False));print(R[R.gene=='SLC16A3'].to_string(index=False))\n",encoding='utf-8')

