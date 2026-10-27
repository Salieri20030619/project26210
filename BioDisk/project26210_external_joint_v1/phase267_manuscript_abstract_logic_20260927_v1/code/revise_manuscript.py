from pathlib import Path
import json,re,copy,hashlib,collections,difflib,shutil
P=Path(__file__).resolve().parents[1];R=P.parent
source=R/'phase266_report_suppfigS1_20260927_v1/report/report_materialized.json'
original=json.loads(source.read_text());b=copy.deepcopy(original)
for i,x in enumerate(b):x['_phase266_index']=i
abstract=[
'肾癌细胞的代谢状态与巨噬细胞异质性及组织分布之间的联系尚需在患者层面系统解析。本研究整合5个单细胞队列的75个文库、40位个体，质控保留258,120个细胞，确定26位患者的11,950个恶性上皮细胞；联合4个配对组织转录组与代谢组队列的255名患者、主要空间队列的37个患者单位及38张切片，并纳入两个独立空间来源的9名病例开展对照。',
'恶性上皮F1程序与TYMP表达正相关（20位患者，rho中位数0.115，17/20同向，P=0.0153，FDR=0.0197），伴随嘧啶、半胱氨酸/蛋氨酸及精氨酸/脯氨酸程序增强。ILF2、HOXD3和HNF4A调控子，以及CASP9、CHMP1A和PRSS8推断活性与F1同向。多胺相关酶和配对代谢组支持MTA与合成背景的联系，SMS与实测MTA的调整后合并rho为0.250（P=0.000578，FDR=0.00289）；TYMP与MTA的直接组织关联未获支持。',
'独立髓系细分在22,444个成熟巨噬细胞中识别出SPP1富集的C4状态，共2,674个细胞；SPP1在24/26名配对患者升高（中位log2差1.449，P=6.68e-5，FDR=5.45e-4）。排除TYMP/SPP1后，胆固醇与溶酶体程序仍在26/26名患者升高。实测MTA与组织SPP1在四队列同向，但仅为名义正关联（rho=0.139，P=0.0293，FDR=0.0976）。',
'空间分析显示，多胺相关转录背景与SPP1相关巨噬丰度在GSE175540的23/23名和TJ-RCC的12/12名患者同向，rho中位数分别为0.226和0.163（P=2.38e-7和0.000488，FDR=7.15e-7和0.00586）。独立病例支持代谢与巨噬状态程序的区域联系；TYMP相关胆固醇程序、SPP1相关区域的溶酶体/脂质处理及免疫调节特征提供功能对应，T细胞绝对丰度未一致下降。空间质谱中精胺质量候选信号与SPP1呈名义正关联（12名患者，rho=0.622，P=0.0347，FDR=0.2081）。',
'这些结果表明，TYMP相关恶性细胞状态伴随核苷、蛋氨酸及多胺代谢改变，并与SPP1相关巨噬的脂质处理功能和组织分布相联系。已有TYMP敲低实验中MTA、蛋氨酸、SAM及亚精胺下降，为这一代谢背景提供方向对应。MTA作为候选联系保留，但现有证据尚未建立TYMP经MTA直接诱导SPP1的完整因果链。'
]
patch={
1:'肾癌单细胞、多组学与空间转录组联合研究　2026-09-27',
8:'本研究采用患者分层的横断面多组学设计，分析恶性上皮细胞的TYMP相关转录状态、核苷及氨基酸代谢、配对组织中的实测MTA、巨噬细胞状态及其空间分布。恶性上皮和髓系细胞分别进行状态识别；组织代谢组、空间转录组和空间质谱按来源独立建模，并在患者层面汇总效应。',
9:'单细胞数据用于表达定量、状态识别和网络构建，独立患者用于效应汇总与统计推断。比较包括连续F1关联、患者内F1主导细胞与其他细胞的比较，以及F1高主克隆与其他主克隆的比较，分别评价程序使用量、主导状态和克隆背景。已有786-O细胞TYMP敲低实验用于代谢方向的对照，其数据未并入公共队列的统计模型。',
35:'独立生物学样本量以患者为单位，细胞数、spot数、切片数及患者×程序表达谱数用于描述观测规模。患者内相关和配对差异汇总后，采用单样本或配对Wilcoxon符号秩检验；各端点使用既定的精确检验或正态近似、连续性校正及零差值处理。F1-CNV的rank筛选端点使用单侧正向检验，其余主要结果使用双侧检验。',
36:'有效统计检验以目标效应的原始P<0.05作为名义显著筛选标准，同时报告效应量、方向一致性及相应FDR。P<0.05而FDR≥0.05称为名义显著；只有通过相应FDR阈值的结果称为校正后显著。BH校正在预先定义的同类假设家族内分别进行，不将差异基因、通路和其他不同类型的检验混入同一校正池；预设的跨方法联合家族仍按原范围校正；单个预设主要端点不另作FDR校正。有效独立患者少于10名时按探索性结果报告方向、效应和一致性；无有效生物学重复的描述性投射不进行推断性检验。',
43:'配对组织转录组与代谢组包括RC18、RC20、CPTAC和CPTAC_val四个队列，分别包含97、52、50和56名肿瘤患者，共255名。依据原始样本映射匹配RNA和实测5′-methylthioadenosine；RNA采用log2(TPM+1)，在队列内计算秩相关和条件关联。模型包括未调整、髓系背景调整及TYMP/MTAP条件模型，分别评价合成与清除背景；跨队列采用随机效应合并估计，并记录各模型的P值及检验家族。',
68:'空间队列先计算片内效应，再在患者层进行精确符号检验及各分析预设的伴随检验；相关和差值的区间采用患者重采样。整片质谱对应采用患者层双侧置换相关及bootstrap。基因排序GSEA的P值描述排序富集统计量，与患者配对程序评分检验分别报告。',
69:'完整统计表记录各端点的有效患者数、效应方向、原始P值、FDR、检验家族及模型协变量；不可估计端点保留其资格说明。不同空间尺度、组成调整和基因集处理的结果分别记录。',
97:'调控子子网络包含71条唯一TF→靶基因边，其中18条由两种算法共同支持。与TYMP的直接预测连接仅见于GRNBoost2的ILF2→TYMP；HOXD3和HNF4A均未检出TYMP边。预设HIF1A在两套pySCENIC调控子矩阵中均不可评估，补充候选评估也未获得该调控子。',
126:'三种轨迹方法对恶性上皮细胞状态变化的评价显示方法间差异。Monocle3提供状态排序的图谱背景（图3B），患者内连续关联与末端比较见表6及附表S6-S8。',
134:original[134]['text'],
136:'TYMP相关恶性上皮状态同时具有半胱氨酸/蛋氨酸、精氨酸/脯氨酸程序增强及ODC1、SRM表达变化（图3A、图3C、表5及附表S5），将核苷代谢状态与SAM供给和多胺合成背景联系起来。',
137:'SRM和SMS分别参与亚精胺及精胺生成，其氨丙基转移反应同时产生5′-甲硫腺苷（5′-methylthioadenosine，MTA）；MTA经MTAP参与的分解过程进入蛋氨酸补救循环[1-3]。配对组织RNA与实测代谢组用于评价这一候选代谢物与合成、清除背景的对应。',
144:'无监督髓系细分在39名患者的22,444个成熟巨噬细胞中识别出12个候选状态（图3D）。这些状态通过多基因表达鉴定，为分析肿瘤代谢背景所对应的髓系异质性提供了细胞层面的依据。',
146:'C4、C10和C15合计包含4,502个细胞。在16名同文库配对患者中，三群联合占成熟巨噬细胞的比例中位数为9.76%，C4单群为7.41%；三群联合贡献捕获SPP1原始RNA的份额中位数为59.63%，C4为54.59%（图3E）。SPP1信号因而集中于少数状态，同时具有跨状态分布。',
149:'C4包含2,674个细胞，分布于37名患者和五项研究，表现为GPNMB、CCL18、CYP27A1、LPL、APOE与SPP1共同富集（图4A）。患者配对分析进一步确认了这一表达组合。',
157:'实测MTA与组织SPP1在四个配对转录组/代谢组队列中均呈正向联系。255名患者的未调整合并rho=0.13894，P=0.02926，原多候选检验家族FDR=0.09755，属于名义正关联。髓系背景调整及加入TYMP条件后的模型结果见附表S11；组织MTA及SPP1的排序分布见图3F。',
180:'可估计轨迹中的C4位于单核之后，但与其他成熟巨噬的相对排序没有统一方向（图4D、图4E）。巨噬状态C4与图2D患者×程序表达谱的C1-C7展示簇属于不同分类体系。',
192:'GLMF1/2显示TYMP相关恶性成分、代谢转录背景与SPP1相关巨噬在组织中的区域性分布。两例病例图用于展示患者层统计关系对应的组织结构。',
194:'GLMF1/2的组织组成、TYMP表达、代谢背景和SPP1-state丰度见图5A-C；代谢与去SPP1状态程序的高值重叠见图5D。四条代谢通路及外部ATAC参考活性的空间投射另见补充图S1。',
198:'这些结果支持SPP1相关巨噬与血管成分的数量共变化，但未显示统一的SPP1状态特异血管富集。成纤维及ECM程序的邻近检验也未获得一致正向支持，当前数据不支持普遍的SPP1-CAF屏障模式。',
205:'由SPP1巨噬和血管成分构成的区域指标与T细胞相对平衡呈负向联系。移除T细胞参照中的SPP1成分后，TJ-RCC仍有10/12名患者同向，两种固定区域尺度的相关系数中位数分别为-0.490和-0.387，P均为0.03857，FDR均为0.05786。该名义关联不完全由参照中的共同SPP1成分解释。',
208:'独立空间来源显示相同的区域性代谢-巨噬程序联系。在E-MTAB-12767的六名患者中，代谢转录背景与去SPP1状态程序在固定区域尺度全部同向，相关系数中位数为0.493。有足够TYMP区域变异的三名患者中，TYMP与代谢背景及状态程序也均同向。GSE294006的三个治疗后病例中，代谢背景与SPP1 RNA全部正向，相关系数中位数为0.559；状态程序对应有一例反向。',
214:'该联系对归一化方式敏感：以总强度之比计算时，方向仍为正，但P=0.05605、FDR=0.3363。局部区域的预设质谱与转录对应未获统一支持，TYMP及SPP1巨噬内部比例也未同步重复。因此，精胺候选信号提供的是整片组织层面的探索性联系。',
233:'TYMP催化胸苷与无机磷酸生成胸腺嘧啶和2-脱氧核糖-1-磷酸的可逆反应[7]。本研究显示，TYMP相关恶性上皮状态同时伴随蛋氨酸及部分多胺程序改变，表明其对应的代谢特征超出单一核苷反应。F1、代谢通路评分及ODC1/SRM表达从转录状态、程序和酶基因三个层面支持核苷与氨基酸/多胺背景的共同变化。',
234:'蛋氨酸/SAM供给与多胺合成构成这一代谢联系的不同环节：前者提供氨丙基转移的代谢背景，SRM/SMS参与的合成反应则伴随MTA生成。MTA因而是连接这些代谢过程的具体候选产物，其组织水平还受到清除能力的影响。',
237:'细胞内、胞外和组织平均MTA反映不同代谢区室。Barekatain等在胶质母细胞瘤模型中观察到MTAP缺失细胞释放MTA，而MTAP表达的基质细胞可代谢胞外MTA[6]。这一结果说明周围细胞的清除能力能够改变局部暴露，为培养体系与组织平均读数的差异提供代谢交换层面的解释。',
238:'直接扰动研究表明MTA能够影响髓系响应。Hansen等通过MTA处理及受体/STAT3干预观察到信号依赖的表型变化[4]；Keyel等发现MTA与TLR刺激背景共同改变炎症和共刺激响应[5]。这些结果为组织MTA与巨噬状态的联系提供生物学依据，也表明其效应依赖接收细胞和刺激环境。',
241:'胆固醇处理与溶酶体程序的一致增强表明，C4具有脂质加工和细胞内降解相关的功能特征。该状态同时呈现抗原呈递、炎症和干扰素程序降低。吞噬体在基因排序富集中增强，却未在患者配对评分中一致升高，显示不同功能读出不能合并为笼统的“巨噬活化”。',
244:'综合肿瘤细胞、配对代谢组、巨噬状态和空间结果，可提出代谢环境与髓系功能相联系的工作模型：TYMP相关恶性细胞状态伴随核苷和蛋氨酸/多胺程序改变；MTA与生成、清除背景及组织SPP1存在分段联系；SPP1富集巨噬具有脂质处理和溶酶体特征，并在组织中与相关代谢区域共同分布。',
245:'这一模型将代谢候选、细胞状态身份及组织位置联系起来。独立空间病例的程序共变化和TYMP相关胆固醇背景进一步支持功能层面的对应；局部代谢暴露与巨噬响应之间的直接关系则需要在同一实验体系中检验。',
247:'本研究的主要限制在于不同层次的关联尚未建立连续因果关系。TYMP经典反应不直接生成MTA；公共TYMP-MTA直接关联及同文库肿瘤TYMP与整体巨噬SPP1的联系均不稳定，串联模型的间接效应区间跨零。相同体系内的TYMP相关代谢改变、胞外MTA暴露和髓系SPP1响应仍需联合测定。',
248:'组织MTA和组织SPP1混合了多种细胞来源，不能分别视为肿瘤细胞分泌量和巨噬专属表达。外部髓系扰动中的ARG1/VEGFA等终点不等同于SPP1诱导，受体或STAT3功能依赖亦不等同于普遍的磷酸化升高。其他癌种的代谢交换结果提供机制背景，并非本队列的局部暴露测量。',
250:'通路、调控子、推断蛋白活性及代谢代理均为模型定义的读出，与酶促速率、代谢物浓度、蛋白分泌和实测免疫功能不同。CNV树描述拷贝数状态关系，并非突变验证的祖先关系；F1选择使用了TYMP相关信息，且研究构成不均衡。主克隆划分对距离阈值的敏感性及各端点的有效患者数限定了结果的适用范围。',
251:'已有TYMP敲低实验覆盖细胞内代谢变化，尚未包含胞外MTA、MTA补回及巨噬SPP1干预终点。实验方向用于与公共队列发现比较，未用于公共数据的候选筛选或效应估计。',
253:'已有786-O细胞实验显示，相对O-NC，shTYMP组的5′-MTA、Methionine、SAM、Thymine和Spermidine下降，Thymidine升高。胸苷（Thymidine）与胸腺嘧啶（Thymine）的相反变化与TYMP对应的核苷反应方向一致；蛋氨酸、SAM和亚精胺下降则与TYMP相关细胞状态中的氨基酸/多胺程序相对应（表12）。',
254:'表 12 TYMP敲低后的代谢变化与公共多组学结果',
256:'实验列依据已提供的TYMP敲低结果记录变化方向；本研究未对该实验重新计算效应量或显著性。',
257:'MTA下降使候选代谢物与TYMP干预后的代谢响应形成方向对应。这一实验结果与公共队列中的程序关联共同支持核苷及蛋氨酸/多胺背景的联系，但两类数据分别回答干预响应和患者间状态关联问题。',
259:'本研究确立了TYMP相关恶性上皮状态的核苷、蛋氨酸及部分多胺代谢特征，识别出SPP1富集且具有胆固醇处理和溶酶体程序的巨噬状态，并在多来源空间组织中观察到代谢背景与SPP1相关区域的共同分布。配对代谢组和实验方向将MTA保留为候选代谢联系。上述结果为肿瘤代谢与巨噬功能适应提供了患者层面和组织层面的证据，并明确了检验直接代谢介导关系的具体对象。',
261:'空间结果将单细胞发现延伸到组织结构：多胺相关转录背景与SPP1相关巨噬丰度及状态程序在区域尺度共同出现，独立病例也显示代谢背景与去SPP1状态程序同向。该联系在多基因程序层面比TYMP或SPP1单基因更一致，支持从代谢环境与巨噬功能组合解释局部分布。',
263:'单细胞C4来自本研究的髓系状态发现，空间SPP1-state来自既有参考映射，两者定义不同。空间关联支持SPP1相关巨噬背景的组织对应，尚不能等同于C4的独立空间定位。',
265:'SPP1相关区域同时具有血管绝对丰度、免疫调节及侵袭程序的变化，提示其处于复杂的组织微环境。受体表达对应、T细胞相对分布与CD8区域距离分别描述信号兼容性、组成和局部结构；T细胞绝对丰度没有一致下降，因而不能将相对分布直接解释为淋巴细胞排斥。',
266:'血管和免疫联系随细胞组成、来源和空间尺度变化。当前数据不支持已确认的CD8排斥、普遍CAF屏障或SPP1特异血管重塑；这些区域特征与临床结局、治疗反应的关系未在本研究中检验。',
268:'空间中的多胺及MTA相关转录程序反映合成背景，而非实测MTA浓度。独立病例的调整结果及配对组织对照显示，状态程序与单基因表达具有不同变化模式；空间共变化不确定TYMP、MTA和SPP1之间的诱导顺序。',
269:'空间映射、混合spot程序和表达区域距离受组成、测序深度及共同基因影响，部分关联在调整后减弱。病理边界和空间配准尚缺少独立审核，C4的空间混合恢复未达到全部预设门槛。补充图S1的ATAC活性来自单个外部供者的chr1四motif试点，属于细胞组成加权的描述性投射，不提供切片内直接ATAC测量。',
273:'这些干预证据来自不同疾病和细胞模型，与本研究的区室内表达关联及空间程序对应共同提示，核苷稳态、胆固醇处理和溶酶体功能可构成肿瘤代谢与巨噬表型之间的功能联系。',
340:'主图1-5分别展示整合细胞图谱、恶性上皮调控状态、代谢及巨噬信号来源、C4功能和髓系排序，以及组织空间对应。补充图S1展示两例空间代谢通路及外部ATAC足迹活性投射。图注说明各面板坐标、颜色和统计单位，完整输入身份与源表见配套数据清单。',
344:'各模型使用同一组255名患者。多候选家族FDR与单个预设MTA端点的背景调整检验分别报告，二者对应不同检验家族。',
346:'成熟巨噬状态比较使用22,444个细胞；含单核的轨迹和调控分析使用33,098个细胞。C4标志比较、肿瘤动态主模型及同文库分析分别包含26、14和16名有效患者，各端点采用相应资格集合。',
353:'公共转录组、配对实测代谢组、空间资料及统计源文件的身份对应见配套来源清单。',
}
for i,t in patch.items():b[i]['text']=t
b[159]['headers'][-1]='证据层次'
b[159]['rows'][3][0]='本研究单细胞发现';b[159]['rows'][3][1]='12个成熟巨噬候选';b[159]['rows'][3][3]='定位SPP1相关巨噬状态'
b[255]['caption']=patch[254];b[255]['headers']=['TYMP敲低后的变化','公共多组学对应','证据层次'];b[255]['rows'][-1][0]='未提供方向的其他实验端点'
# Remove editing-history wording from figure legends while retaining their scientific definitions.
legend_patch={
39:('不是新增聚类边界','不表示离散聚类边界'),
40:('基因、行序与原图一致。',''),41:('不再重复嵌入Cell identities总览。',''),
110:('保留原注释与坐标。',''),111:('主图仅稀疏显示1、5、10、15、20号标签以避免拥挤，不删除其余染色体。','标签标出1、5、10、15、20号染色体。'),
112:('未重新拟合cNMF、UMAP或表达值；F1–F7颜色沿用原图。','颜色分别表示程序身份、使用量和归一化表达。'),
114:('顺序、值及原行树保持。',''),115:('71个节点及71条有向连接全部保留。','共71个节点及71条有向连接。'),116:('未重新估计KDE，',''),
164:('原Top20通路及排序保留','Top20通路按预设顺序展示'),165:('使用原UMAP及原全局伪时间，未重拟合。','使用UMAP坐标及全局伪时间。'),167:('将原两幅巨噬定位展示合成一幅，保留12个冻结候选的密度底图，并叠加','12个巨噬候选的密度底图叠加'),
184:('基因、状态、顺序及数值不变。',''),186:('仍为原C4对其他成熟巨噬细胞GSEA','C4与其他成熟巨噬细胞的GSEA'),188:('不重新平滑或拟合。',''),
}
for i,(a,z) in legend_patch.items():b[i]['text']=b[i]['text'].replace(a,z)
b[167]['text']=b[167]['text'].replace('细胞成员、坐标和密度算法未改变。','')
b[186]['text']=b[186]['text'].replace('，未改成三群联合分析','').replace('；未重新运行GSEA','')
# Sequential section numbers, preserving a mapping to the previous manuscript.
def block(start,end):return b[start:end]
method_blocks=[(7,10),(10,16),(16,19),(19,22),(22,25),(25,29),(29,34),(42,45),(45,50),(50,54),(54,59),(70,73),(73,75),(59,63),(75,77),(63,67)]
method=[];number_map=[]
for n,(a,z) in enumerate(method_blocks,1):
 chunk=block(a,z);old=chunk[0]['text'];chunk[0]['text']=re.sub(r'^1\.\d+',f'1.{n}',old);number_map.append(dict(old=old,new=chunk[0]['text']));method+=chunk
method += [dict(type='h2',text='1.17 ATAC足迹活性的空间投射'),dict(type='p',text='以GSE207493的RCC100 ATAC片段为外部参考，按宽谱系归类后在chr1共同开放峰上构建HINT-ATAC足迹；端点兼容输入保留观测切割位置和片段长度，并与直接片段计数核对。固定HNF4A、CEBPD、JUNB和STAT3四个JASPAR motif，经匹配与偏倚校正获得保护分数和周围切割信号，其和定义为活性。将内皮、恶性、髓系和平滑肌/周细胞的活性按spot内条件丰度加权；四类代表丰度占全部参考丰度至少80%时保留投射，否则标记为缺失。GLMF1/2的四条代谢通路采用既有ssGSEA分数，与细胞丰度和ATAC投射并列展示（补充图S1）。'),dict(type='h2',text='1.18 统计分析'),b[35],b[36],b[68],b[69]]
number_map += [dict(old=original[34]['text'],new='1.18 统计分析'),dict(old=original[67]['text'],new='1.18 统计分析')]
result_spans=[(78,85),(85,91),(91,98),(98,108),(117,125),(125,135),(135,143),(143,148),(148,156),(156,162),(170,176),(176,182),(190,195),(215,220),(195,200),(200,204),(204,207),(220,222),(207,212),(212,215)]
results=[];rm={};remove={90,155,161,181,218}
for n,(a,z) in enumerate(result_spans,1):
 chunk=[x for x in block(a,z) if x['_phase266_index'] not in remove];old=chunk[0]['text'];new=re.sub(r'^2\.\d+',f'2.{n}',old)
 if a==148:new=f'2.{n} C4具有跨患者一致的SPP1与脂质处理基因特征'
 if a==156:new=f'2.{n} 实测MTA与组织SPP1呈名义正向关联'
 chunk[0]['text']=new;number_map.append(dict(old=old,new=new));rm[old.split()[0]]=new.split()[0];results+=chunk
 if a==78:results+=block(37,42)
 if a==98:results+=block(108,117)
 if a==156:results+=block(162,170)
 if a==176:results+=block(182,189)
results+=block(222,231)
# Discussion follows findings, experimental correspondence, integration, limitations, conclusion.
disc=[]
def add_disc(title,indices):
 disc.append(dict(type='h2',text=f'3.{1+sum(x["type"]=="h2" for x in disc)} '+title));disc.extend(b[i] for i in indices)
add_disc('TYMP相关状态连接核苷与蛋氨酸/多胺代谢',[233,234])
add_disc('TYMP敲低的代谢变化与转录程序方向相符',list(range(253,258)))
add_disc('MTA的生成与清除背景为髓系联系提供依据',[236,237,238])
add_disc('SPP1富集巨噬呈现脂质处理与溶酶体特征',[240,241,242])
add_disc('空间共分布连接代谢背景与巨噬相关区域',[261,262,265])
add_disc('核苷稳态与胆固醇和溶酶体功能相联系',[272,273])
add_disc('肿瘤代谢与髓系功能适应的工作模型',[244,245])
add_disc('研究局限与适用范围',[247,248,249,250,251,263,266,268,269,270])
add_disc('结论',[259])
# References follow supplementary tables, not the middle of their index.
appendix=block(274,345)+[dict(type='h2',text='4.12 数据来源与分析集合'),b[346],b[353]]+block(367,370)
refs=[dict(type='h1',text='参考文献')]+[b[i] for i in list(range(348,353))+list(range(354,367))]
new=b[:3]+[dict(type='p',text=t,_new_section='abstract') for t in abstract]+[b[6]]+method+[b[77]]+results+[b[231]]+disc+appendix+refs+[dict(type='h1',text='补充图')]+block(370,len(b))
# Update numeric section cross-references only, never statistical decimal values.
for x in new:
 if 'text' in x:
  x['text']=re.sub(r'(正文\s*)(2\.\d+)',lambda m:m[1]+rm.get(m[2],m[2]),x['text'])
  x['text']=x['text'].replace('本项目','本研究')
  for a,z in [('本轮校正','相应家族校正'),('正式分析集','分析集'),('原冻结纳入范围','预设纳入范围')]:x['text']=x['text'].replace(a,z)
# Editorial audit, quantitative tables and figure identity are checked independently.
oldtables=[x for x in original if x['type']=='table'];newtables=[x for x in new if x['type']=='table'];assert len(oldtables)==len(newtables)==32
for old,nw in zip(oldtables,newtables):
 if old.get('_phase266_index') is not None:pass
 # Table12 moves within Discussion; table order otherwise remains identical.
 assert re.findall(r'[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?',json.dumps(old['rows'],ensure_ascii=False))==re.findall(r'[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?',json.dumps(nw['rows'],ensure_ascii=False)),old['caption']
assert [x['figure'] for x in new if x['type']=='figure']==[1,2,3,4,5]
assert sum(x['type']=='supplementary_figure' for x in new)==1
for prefix,expected in [('1',18),('2',20),('3',9)]:
 ns=[int(re.match(r'^'+prefix+r'\.(\d+) ',x['text'])[1]) for x in new if x['type']=='h2' and re.match(r'^'+prefix+r'\.\d+ ',x['text'])]
 assert ns==list(range(1,expected+1)),ns
for phrase in ['客户','报告由','置于实验对应段','阶段结论','构成当前报告的论证顺序','接下来','下一步']:
 assert phrase not in json.dumps(new,ensure_ascii=False),phrase
changes=[]
for x in new:
 if '_phase266_index' in x:
  i=x['_phase266_index'];prev=original[i]
  if x.get('text')!=prev.get('text') or x.get('headers')!=prev.get('headers') or x.get('rows')!=prev.get('rows'):
   changes.append(dict(source_index=i,before=prev,after=x))
(P/'report/report_materialized.json').write_text(json.dumps(new,ensure_ascii=False,indent=2))
(P/'qc/EDITORIAL_CHANGES.json').write_text(json.dumps(changes,ensure_ascii=False,indent=2))
(P/'qc/SECTION_NUMBER_MAP.json').write_text(json.dumps(number_map,ensure_ascii=False,indent=2))
(P/'qc/STRUCTURAL_VALIDATION.json').write_text(json.dumps(dict(status='PASS',tables=32,table_numeric_values_unchanged=True,main_figures=5,supplementary_figures=1,methods=18,result_sections=20,discussion_sections=9,abstract_characters=sum(map(len,abstract)),new_statistical_tests=0,source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest()),indent=2))
(P/'delivery/26210_revised_abstract_CN.md').write_text('# 摘要\n\n'+'\n\n'.join(abstract))
(P/'delivery/26210_全文逻辑调整说明.md').write_text('# 摘要与全文逻辑修订\n\n摘要改为研究问题、资料与核心发现、空间对应及结论，删除文档编排和客户沟通措辞；补充空间队列、独立病例、胆固醇/溶酶体功能和精胺候选联系。\n\n方法按数据和分析对象排列，统计方法集中于1.18；补充ATAC投射方法。结果保持肿瘤代谢—候选MTA—巨噬发现及联系—空间对应的顺序，把脂质/溶酶体空间功能放在区域共分布后、独立病例和质谱前。原未纳入的scDRS分析没有恢复，结果节重新连续编号。讨论依次为主要代谢发现、敲低实验对应、MTA、巨噬状态、空间及功能、工作模型、局限、结论。完整统计附表后集中排列参考文献，补充图S1及图注保留。\n\n32张表的统计数值、5张主图、补充图S1及原统计家族保留；本轮仅进行文字和结构修订。\n\n## 章节对照\n\n'+'\n'.join('- '+x['old']+' → '+x['new'] for x in number_map))
for n in range(1,6):shutil.copy2(R/f'phase263_spatial_reintegrated_report_20260926_v1/figures/Figure{n}.pdf',P/f'figures/Figure{n}.pdf')
shutil.copy2(R/'phase266_report_suppfigS1_20260927_v1/figures/Supplementary_Figure_S1.pdf',P/'figures/Supplementary_Figure_S1.pdf')
print('EDITORIAL_READY',len(new),'blocks',len(changes),'changed blocks',sum(map(len,abstract)),'abstract characters')
