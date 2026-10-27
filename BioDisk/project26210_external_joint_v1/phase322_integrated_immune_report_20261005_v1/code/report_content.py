from pathlib import Path
import json, pandas as pd, numpy as np
P=Path(__file__).resolve().parents[1]
manifest=json.loads((P/'resources/STATISTICAL_SOURCE_MANIFEST.json').read_text())
blocks=[];used=[]
def para(text):blocks.append({'type':'p','text':text})
def head(text,level=1):blocks.append({'type':'h','text':text,'level':level})
def table(title,headers,rows,widths=None,source=''):
 blocks.append({'type':'table','title':title,'headers':headers,'rows':[[str(v) for v in r] for r in rows],'widths':widths,'source':source})
def data(stage,name):
 rows=[x for x in manifest if x['stage']==stage and Path(x['source']).name==name];assert len(rows)==1,(stage,name,len(rows));used.append(rows[0]);return pd.read_csv(P/rows[0]['relative'],sep='\t')
def f(v,d=3):return 'NA' if pd.isna(v) else f'{float(v):.{d}f}'
def pv(v):return 'NA' if pd.isna(v) else f'{float(v):.5g}'
def one(df,**kwargs):
 for k,v in kwargs.items():df=df[df[k]==v]
 assert len(df)==1,kwargs
 return df.iloc[0]
trans={'Glycolysis':'糖酵解','Interferon Alpha Response':'IFNα反应','Interferon Gamma Response':'IFNγ反应','Inflammatory Response':'炎症反应','TNF-alpha Signaling via NF-kB':'TNFα/NFκB',
'GO:0002474__antigen processing and presentation of peptide antigen via MHC class I':'MHC-I抗原呈递','GO:0002495__antigen processing and presentation of peptide antigen via MHC class II':'MHC-II抗原呈递','GO:0006909__phagocytosis':'吞噬相关程序',
'UV Response Up':'UV响应上调程序','Unfolded Protein Response':'未折叠蛋白反应','mTORC1 Signaling':'mTORC1程序','CD8_share_of_CD8_plus_macrophages':'CD8与巨噬捕获组成'}
blocks.append({'type':'title','text':'TYMP与SPP1相关免疫特征综合报告'})
para('26210项目  |  数据与统计整合截至2026年10月5日')
head('摘要')
para('TYMP与SPP1在巨噬细胞内具有可重复的表达联系，同时对应不同的免疫转录特征。TYMP与IFNγ、MHC-I抗原呈递及吞噬相关程序正向关联，SPP1与IFNγ及MHC-II抗原呈递程序负向关联；主要方向在控制另一基因、严格细胞成员和免疫程序去重叠后保持。两基因与糖酵解程序的正向联系仍在，但效应强度受到少数基因影响。')
para('跨细胞和组织层结果呈现不同含义。44个合格匿名供者的组织bulk中，SPP1同时与细胞毒和功能障碍程序正相关；细胞来源限定的CD8、Treg及NK分析未形成统一免疫抑制关联。部分肾癌空间分析支持肿瘤表达区内侧T/CD8信号较低，以及局部基质状态与相对T细胞不足的方向联系，但未稳定证实SPP1特异的物理排斥。')
para('公开SPP1干预数据提供了补充实验依据：特定模型中CAF减少、CD8数量或增殖升高，另一些模型则表现为基质程序下降但免疫端变化不一致。免疫治疗队列中，SPP1相关转录指标与部分生存结局及治疗交互有关；这些临床指标不能等同于TYMP作用、C4细胞数量或已可用的耐药预测工具。综合证据支持TYMP与SPP1相关巨噬状态及免疫过程的参与性联系，完整肿瘤来源TYMP到免疫屏障的因果链尚未建立。')

head('1 方法')
head('1.1 数据层级与独立样本',2)
para('分析以独立患者、动物或独立混样池为统计单位。细胞、spot、空间bin、组织视野和技术分区用于生成测量值，不替代生物学样本量。同患者多文库先按原方案汇总；基线与进展期分别分析。不同队列和不同分析的人数不相加。')
table('表1 主要数据与有效分析范围',['数据来源','独立单位与有效范围','测量内容'],[
['原75文库单细胞集合','258120细胞；C4冻结2674细胞；巨噬内关联最多31患者，CD8接收端16患者','巨噬状态、来源及免疫接收端'],
['Obradovic未治疗ccRCC','11患者；13个CD45+文库；8患者有CD45−配套','巨噬内联系及同患者跨细胞程序'],
['HCRN GU16-260','13位ccRCC患者；基线8、进展6份文库；有效巨噬配对随资格变化','治疗分层方向对照'],
['公开组织bulk PRJNA1223824','44个合格匿名供者','SPP1及文献五基因与组织免疫程序'],
['肾癌Visium','TJ-RCC 12患者；GSE175540 23患者24切片；GLMF 2患者','界面、邻域、基质及T/CD8空间分布'],
['高分辨率肾癌空间案例','10x数据1供者；另有预选命名区域试点','逐细胞几何、SPP1来源及分割敏感性'],
['外部实验与源表','HCC空间15患者；PDAC每组3混样池；CRC各面板5至11只动物/组','SPP1干预、基质、CD8功能及定位'],
['CheckMate与JAVELIN','RNA匹配311人与726人；具体模型按完整病例确定','免疫治疗相关生存和治疗交互']
],[3.5,6.2,6.3],source='各分析原始资格与模型成员表')
para('原单细胞接收端分析中，CD8两个功能终点为16患者、19个同肿瘤文库、4项研究；Treg/CD4和NK端点各13患者。整体巨噬TYMP—SPP1为29患者，C4内部为15患者。Obradovic的11患者用于后续独立队列验证及免疫程序拓展，这些后续模型共享同一队列。HCRN基线/进展文库有一位患者重复取材，分层结果分别保留。')
head('1.2 细胞身份与指标定义',2)
para('原C4成员及CD8、常规CD4、Treg、NK接收端沿用冻结版本，未将混合T/NK群直接当作NK。Obradovic和HCRN采用原始计数质控、髓系细分和多标记资格；其宏观谱系和后续宽松/严格成员分别记录，不重新指定C4。空间旧SPP1相关状态、项目C4、SPP1单基因和文献五基因评分作为不同指标分析。')
para('文献五基因为SPP1、FN1、CLEC5A、OLR1、IL1RN。组织TPM/FPKM先取log2(x+1)再平均；已标准化的CheckMate表达直接平均。单细胞伪bulk使用相应原计数归一化定义。细胞毒程序采用NKG7、GNLY、GZMB、PRF1、CTSW；功能障碍相关程序采用PDCD1、CTLA4、LAG3、HAVCR2、TIGIT、TOX、CXCL13、ENTPD1。Treg比例分母为冻结常规CD4与Treg之和。')
head('1.3 巨噬细胞内关联及最终敏感性分析',2)
para('巨噬细胞内分析先在每位患者中计算秩残差相关，控制测序量、检出基因数、文库和冻结髓系亚群，再汇总患者方向与中位效应。无方向探索使用固定50个Hallmark和6个GO程序，均移除TYMP/SPP1。宽松和严格细胞成员分别计算，56程序、2个源基因和2种成员构成224项免疫检验家族。')
para('最终定向核查限定10条已有关系，加入另一基因的表达秩作为协变量。条件模型20项单独BH；代谢与免疫去重叠、免疫程序之间去重叠共36项单独BH。逐一剔除程序中每个基因，并依据宽松成员的单基因衰减程度固定选取5个基因联合剔除；严格成员沿用同一名单。联合五基因检查只报告效应和方向。原224项发现家族的统计完整保留。')
head('1.4 跨细胞免疫联系与定向通信',2)
para('原接收端分析使用同一批合格肿瘤文库配对来源和接收指标，再汇总到患者。研究调整采用偏Spearman、99999次研究内置换和5000次患者分层bootstrap。文献五基因适配保持同患者、同接收端，仅改变来源侧表征。')
para('Obradovic跨患者模型控制作者处理批次和pT分层，HCRN按治疗时点分开。定向接收端另测PDCD1/ENTPD1共同检出、TOX检出和IFNG/GZMB程序；CD44敲除响应程序来源于真实纯化CD8实验，经一对一人鼠同源映射后固定。自由探索对完整50个CD8 Hallmark程序及CD8/(CD8+合格巨噬)捕获组成分别建模。')
para('定向通信沿用冻结CellChatDB中的9个SPP1受体映射；复合体要求同一接收细胞联合检出。受体表达容量与来源SPP1分别测量，未将含SPP1的乘积分数再与SPP1相关。NicheNet先验靶表达及间质接收网络按原方案单列。')
head('1.5 组织免疫与空间分析',2)
para('公开bulk的60文库先按固定SNP面板核查供者关系；观察文库去除一份低信息样本并合并重复后为44个匿名供者。组织程序用Spearman、100000次置换和5000次供者bootstrap，两个主要端点与两个文献指标端点分别BH。')
para('空间分析分别沿用表达映射界面、计算肿瘤样区域以及逐细胞癌周区域定义。Visium以切片效应再合并到患者，采用固定宽度、邻域及深度、组成和几何调整；部分分析另用协变量约束空间零模型。高分辨率案例使用作者标签/固定多标记及真实坐标，分别计算最近距离、邻域密度与调整后秩相关。')
head('1.6 文献复算与公开扰动',2)
para('CRC原文复算使用作者逐spot源表、GO统计以及经坐标、标签和计数指标逐项确认的P7镜像对象。公开标记迁移到肾癌时另行固定评分和端点，保留原切片校准结果。外部功能研究区分作者实验描述、可获得源数值复算和原始计数重分析。多视野数据先按动物汇总；全身/骨髓Spp1干预、CD44干预和重组蛋白刺激分别报告。')
head('1.7 免疫治疗临床模型',2)
para('CheckMate使用原PFS/OS事件编码，JAVELIN使用1−PFS_P_CNSR。文献评分按规定标准化和固定中位数分组，连续HR按1 SD报告。Cox模型按预设临床协变量调整；评分与治疗交互仅在相应随机两臂范围内检验。CheckMate性别编码统一后的修正模型用于正式汇总。')
para('C4相关转录评分由24患者、5研究配对伪bulk发现，排除SPP1、TYMP及文献五基因等，临床共同覆盖19基因。发现时不使用临床结局，采用留一研究内部验证。Nivolumab PFS比例风险诊断后登记以随访月数加1取自然对数的线性项描述风险系数变化的模型，以事件风险集拟合时间变化系数，并进行预设试验与留一患者检查。')
head('1.8 统计呈现与版本继承',2)
para('有有效推断检验时，原始P<0.05作为名义候选门槛，同时报告FDR与检验家族；FDR<0.05才称校正后显著。n<10的患者/池层结果按原方案报告方向和效应，不以细胞或区域补n。不同估计对象和不同家族的P/FDR不互相替换，原始值与修订值均保留来源对应。本报告整合既有分析，不重新拟合生物学模型。')

head('2 结果')
head('2.1 TYMP与SPP1的巨噬细胞内联系在独立队列重复',2)
para('原整体巨噬集合中，29患者有24位TYMP—SPP1正向，中位患者内ρ=0.049133，P=0.000546、FDR=0.002184；C4内部15患者中13位正向，ρ=0.134872，P=0.007385、FDR=0.029541。Obradovic独立11患者中10位正向，严格成员核心模型ρ=0.046985，加入干扰素/缺氧后ρ=0.071659，两模型P=FDR=0.01171875。不同模型重复使用同一队列，未合并成更大的统一样本。')
para('TYMP还与输出/处理过程相关。原整体巨噬TYMP—TGFB1在31患者中26位正向，ρ=0.052605，P=0.000192、FDR=0.000769；后续Obradovic宽松成员中，TYMP—CTSB为10/11正，ρ=0.113757，P=0.011719、FDR=0.014423，TYMP—TGFB1两个模型均11/11正。C4溶酶体状态与CTSB、CTSD、LGALS3和APOE呈正向联系；SPP1对应为名义显著，P=0.035156、FDR=0.072327。相关程序已去除这些输出基因。')
para('同患者来源分析未获得相同的肿瘤端支持。共同13患者的肿瘤TYMP—巨噬SPP1调整ρ=0.013862、P=FDR=0.98529；Obradovic的8例肿瘤样细胞配套也未显示明确正向连接。该来源区别与髓系内在阳性结果同时保留。')
head('2.2 TYMP与SPP1对应不同的免疫转录特征',2)
d=data(314,'ALL_EFFECT_SUMMARIES.tsv')
selected=['Glycolysis','Interferon Alpha Response','Interferon Gamma Response','GO:0002474__antigen processing and presentation of peptide antigen via MHC class I','GO:0002495__antigen processing and presentation of peptide antigen via MHC class II','GO:0006909__phagocytosis','Inflammatory Response','TNF-alpha Signaling via NF-kB']
rows=[]
for target in selected:
 for source in ['TYMP','SPP1']:
  z=one(d,branch='immune',condition='L2_broad',source=source,target=target)
  rows.append([source,trans[target],f'{int(z.positive)}/{int(z.patients)}',f(z.median_rho),pv(z.P),pv(z.FDR)])
table('表2 巨噬细胞内宽松成员的主要免疫关系',['来源','程序','正向患者','中位ρ','P','FDR'],rows,[1.4,5.0,2.2,2.0,2.7,2.7],source='原224项免疫发现家族，11患者；Phase314完整表')
para('TYMP与IFNα、IFNγ、MHC-I和吞噬程序均11/11患者正向；SPP1与IFNα、IFNγ及MHC-II均11/11负向，主要方向在严格成员保持。两者均与糖酵解正向，10/11同向，原P=0.01171875、FDR=0.0254854。TYMP与总体炎症程序在宽松版本为正向，但严格版本仅9/11正、P=0.065430、FDR=0.116319；其TNFα/NFκB程序反而10/11负向且FDR显著。')
para('其余固定程序也完整保留。TYMP同时与T细胞正调节和负调节GO集合正相关，而两集合重叠53基因；SPP1与促进T细胞活化程序负向。现有结果显示具体程序的方向分化，不能归并为两个基因都代表统一免疫抑制或广泛炎症增强。')
head('2.3 互相控制后主要关系保持而部分程序证据收窄',2)
d=data(321,'FORMAL_SUMMARIES.tsv');rows=[]
for _,z in d[(d.condition=='L2_broad')&(d.variant=='conditional')].iterrows():
 rows.append([z.source,trans[z.target],f(z.median_rho),f'{int(z.positive)}正/{int(z.negative)}负',pv(z.P),pv(z.FDR)])
table('表3 控制另一基因后的巨噬细胞内关联',['来源','程序','中位ρ','11患者方向','P','FDR'],rows,[1.4,5.0,2.0,2.4,2.6,2.6],source='20项条件模型家族；原224项发现FDR不被替换')
para('条件模型中，两条糖酵解关系仍10/11正，ρ分别为0.102671和0.102083，FDR=0.012336；其余主要干扰素、MHC及吞噬方向均保持，除严格炎症外均通过本轮条件家族FDR。')
rows=[]
for _,z in d[(d.condition=='L2_broad')&(d.variant=='immune_program_nonoverlap')].iterrows():
 rows.append([z.source,trans[z.target],f(z.median_rho),f'{int(z.positive)}正/{int(z.negative)}负',pv(z.P),pv(z.FDR)])
table('表4 删除与其他免疫程序重叠基因后的结果',['来源','程序','中位ρ','11患者方向','P','FDR'],rows,[1.4,5.0,2.0,2.4,2.6,2.6],source='36项去重叠家族；所有条目另控制另一基因')
para('去重叠后，TYMP—IFNγ、MHC-I及吞噬仍11/11正；SPP1—IFNγ及MHC-II仍11/11负，P=0.00097656、FDR=0.00146484，严格成员同向。IFNα只剩17/97个非重叠基因，宽松两目标均仅9/11保持原方向，P=0.065430；炎症程序效应接近零，宽松ρ=0.003014，严格ρ=−0.000692。原完整程序关联保留，但IFNα和笼统炎症不再作为额外独立机制证据。')
r=data(321,'EVIDENCE_RETENTION.tsv');rows=[]
for source,target in [('TYMP','Glycolysis'),('SPP1','Glycolysis'),('TYMP','Interferon Gamma Response'),('SPP1','Interferon Gamma Response'),('TYMP',selected[3]),('TYMP',selected[5]),('SPP1',selected[4])]:
 z=one(r,condition='L2_broad',source=source,target=target)
 rows.append([source+' '+trans[target],f(z.conditional_rho),f(z.five_gene_rho),f'{int(z.five_gene_positive)}正/{int(z.five_gene_negative)}负'])
table('表5 五基因联合剔除的效应检查',['关系','条件模型ρ','剔除后ρ','患者方向'],rows,[7.0,3.0,3.0,3.0],source='宽松成员选择五基因后固定用于严格；描述性检查，无新增P')
para('3036次单基因剔除汇总中，中位方向均未翻转。联合删除影响较大的5个基因后，TYMP—糖酵解由0.102671降至0.021528，SPP1由0.102083降至0.059951。TYMP名单为ENO1、PKM、PPIA、SLC16A3、ISG20；SPP1名单为PKM、SLC16A3、CD44、ENO1、LDHA。糖酵解删除与所选免疫程序共有的10个基因后，两关系仍FDR显著，因此“免疫程序基因重叠”与“少数关键基因影响幅度”是不同检查。')

head('2.4 细胞来源限定的CD8 Treg和NK分析未形成统一抑制关系',2)
d=data(213,'PRIMARY_12_HYPOTHESES.tsv');pred={'mature_SPP1_pseudobulk':'整体巨噬SPP1','C4_fraction':'C4比例','C4_SPP1_pseudobulk':'C4内SPP1'};ep={'CD8_cytotoxicity':'CD8细胞毒','CD8_dysfunction':'CD8功能障碍','Treg_CD4':'Treg/CD4比例','Treg_CD4_ratio':'Treg/CD4','NK_cytotoxicity':'NK细胞毒'}
rows=[]
for _,z in d.iterrows():rows.append([pred.get(z.predictor,z.predictor),ep.get(z.endpoint,z.endpoint),str(int(z.n_patients)),f(z.rho),pv(z.p_value),pv(z.FDR)])
table('表6 原接收端12项患者级关联',['来源指标','接收端','n','调整ρ','P','FDR'],rows,[3.3,4.0,1.0,2.2,2.7,2.8],source='研究调整置换P，按原预设家族校正')
para('表6中整体巨噬SPP1对应Treg和NK的FDR为NA，表示单个预设终点未作多重校正。12项研究调整后的P均大于0.05。C4比例—CD8功能障碍在未调整时ρ=0.59706、P=0.01461，调整研究来源后ρ=0.32526、P=0.35629。文献五基因替换来源指标后，成熟巨噬内或全捕获细胞伪bulk对CD8、Treg及NK的8项检验也均未达到原P<0.05。')
d=data(305,'ALL_PATIENT_ASSOCIATIONS.tsv')
d=d[(d.dataset=='OBRADOVIC')&(d.scope=='Untreated')&(d.model=='prespecified_adjusted')&d.source.isin(['Macrophage:TYMP','Macrophage:SPP1'])&d.target.str.startswith('CD8:')];assert len(d)==4
table('表7 独立未治疗队列的四项CD8功能联系',['巨噬来源','CD8程序','n','调整ρ','P','FDR'],[[z.source.split(':')[1],ep.get(z.target.split(':')[1],z.target),str(int(z.patients)),f(z.rho),pv(z.P),pv(z.FDR)] for z in d.itertuples()],[2.0,5.0,1.0,2.0,3.0,3.0],source='Obradovic 11患者，批次与pT调整；原来源×两程序×两模型家族')
para('独立队列也未支持这四项关联。文献定向端点中，SPP1与CD8的PDCD1/ENTPD1共同检出呈反向名义联系，ρ=−0.727451，P=0.026336、FDR=0.158014；增加CD8深度调整后ρ=−0.327832、P=0.427937。TYMP在深度模型中也有反向名义结果，ρ=−0.710583，P=0.048212、FDR=0.289270。TOX及IFNG/GZMB没有支持性名义信号。')
para('真实CD44敲除响应程序的肾癌检验同样未获名义支持。巨噬SPP1与共培养CD44缺失程序，批次/分期调整ρ=−0.488355、P=0.182251；加入CD8深度后ρ=+0.532710、P=0.174047。静态程序及TYMP相关模型均P≥0.05。小样本HCRN各治疗层只保留方向，不将其计作新的显著验证。')
head('2.5 受体测量与候选输出网络提供有限的通信线索',2)
para('9个预设SPP1受体在CD8、CD4和NK中形成27项登记问题，仅少数达到完整复合体检出与患者资格。CD8的SPP1—ITGA4/ITGB1接收容量为名义反向联系，n=12，ρ=−0.60132、P=0.02413、FDR=0.65151；CD44及其余可评估项未形成校正后支持。NicheNet的SPP1先验靶表达未建立稳定的来源—接收联系。')
para('间质接收网络中，7患者的描述性活化目标集对TGFB1的先验匹配排名为1/797，在巨噬表达合格配体中1/190；SPP1分别为692/797和162/190。该结果支持TGFB1作为SPP1相关状态的候选输出之一，网络分值不是患者干预效应P。C4溶酶体相关LGALS3、CTSB、CTSD和APOE的来源端联系，与其同文库CD8下游未稳定支持的结果共同保留。')

head('2.6 组织免疫程序共同升高与细胞特异阴性并存',2)
d=data(243,'IMMUNE_ASSOCIATIONS_ALL.tsv');d=d[d.specification=='primary_44_donors'];assert len(d)==4
table('表8 44供者组织bulk免疫联系',['来源指标','组织程序','ρ','95%区间','P','FDR'],[[str(z.exposure),str(z.endpoint),f(z.rho),f'{f(z.ci_low)}至{f(z.ci_high)}',pv(z.p),pv(z.FDR)] for z in d.itertuples()],[2.4,3.2,1.8,4.0,2.3,2.3],source='SPP1两端点及文献五基因两端点分别校正')
para('SPP1与组织细胞毒程序ρ=0.670754，与功能障碍程序ρ=0.491896，两者均FDR显著；重复供者代表替换后方向保持。文献五基因与两个程序也均正向。该结果显示组织层免疫程序共同变化，不能改写为细胞毒程序下降。组织bulk与谱系内RNA状态、细胞比例和功能实验测量不同。')

head('2.7 CD8自由探索发现名义候选但尚未跨队列重复',2)
d=data(320,'ADJUSTED_NOMINAL_CANDIDATES.tsv')
table('表9 CD8程序及捕获组成的名义候选',['巨噬来源','接收端读数','调整ρ','P','FDR'],[[z.source.split(':')[1],trans.get(z.endpoint,z.endpoint),f(z.rho),pv(z.P),pv(z.FDR)] for z in d.itertuples()],[1.8,6.4,2.1,2.8,2.9],source='全部n=11；Hallmark与组成家族分别校正')
para('SPP1—CD8 mTORC1和未折叠蛋白反应在研究内留一患者、严格资格和深度检查中保留负向；TYMP—CD8糖酵解为调整后的正向名义候选，但严格程序表达资格下无患者可评估。UV响应程序在严格5例子集中反向，因此不作为主要机制。HCRN进展期5例中，TYMP—糖酵解ρ=+0.6，而SPP1—mTORC1为0，SPP1—未折叠蛋白反应及组成均为+0.3。以上不构成跨队列稳定重复。')
para('这些是来源表达与另一谱系CD8程序的跨患者关系，与表2至表5的巨噬细胞内关系不同。CD8探索未通过FDR，不改变原巨噬细胞内炎症/糖酵解的统计记录；严格资格不足也不能迁移到另一分析层面。')
# PART2_CHUNK_0

head('2.8 空间分区和基质联系提供有限的免疫进入证据',2)
para('在固定400 μm表达界面，GSE175540的CD8谱系RNA内侧低于外侧，深度及总巨噬调整后19/20患者负向，中位差−0.260230，P=0.00004005、FDR=0.00016022；200和800 μm方向保持。总T映射丰度调整后17/23内侧较低，P=0.034690、FDR=0.046253。TJ-RCC调整后CD8为9/9内侧较低，GLMF两例同向，按各自样本量作描述。')
para('界面SPP1富集未解释这种分区。GSE175540中，SPP1局部巨噬占比与CD8内外差的调整中位ρ=+0.049298，P=0.115318、FDR=0.461273；与总T内外差ρ=+0.041952，P=0.404873、FDR=0.809746。较早的SPP1高区到CD8富集区距离增加为名义线索（n=16，中位差+300 μm，P=0.03516、FDR=0.05273），不能替代不同终点的局部耦合结果。')
para('基质端在TJ-RCC提供了患者方向线索：200 μm主尺度中，外侧ECM组织程序与相对肿瘤内T信号在7/9患者负向；加入血管相关组成后9/9负向，中位ρ=−0.144059。因n=9不作患者P，9张切片的严格空间零模型均未达原P<0.05。GSE175540相应总T终点ρ=−0.017816、P=0.890381，CD8终点ρ=+0.107930、P=0.414307，未形成稳定跨队列重复。')
para('在相同9患者、相同区域中对齐SPP1—ECM—T三条边，血管调整后SPP1—ECM为7/9正、ECM—T为9/9负、SPP1—T为6/9负，三条预期方向同时出现于5/9患者。因此分段联系保留，但没有把不同患者或不同模型的阳性边拼成完整链条。较早TLS外围SPP1相对占比较低的线索主要来自一个研究背景，也未被改写为TLS抑制或CD8排斥。')
head('2.9 高分辨率病例与公开标记迁移未稳定确认完整屏障',2)
para('一个肾癌高分辨率供者中，成纤维候选丰富邻域的CD8距癌细胞更远：高、低组距离中位数63.23与37.33 μm，20 μm内存在癌细胞的CD8比例9.79%与29.06%。加入局部癌细胞、总细胞、血管相关细胞密度、RNA量和位置趋势后，ρ从0.205639降至0.022175；236个空间块方向不一致。该面板未测TYMP/SPP1，故仅为末端空间案例。')
para('另一预选区域经Proseg重分割后，SPP1计数主要归于推定肾上皮群（96.69%），而非成熟髓系。放宽身份和数量门槛后可计算部分髓系关系，但SPP1—成纤维占比仅弱正，成纤维—最近CD8距离在加入密度后转向；原分割、核内计数与Proseg之间也不一致。区域与技术重复未计为独立患者。')
d=data(317,'COHORT_SUMMARY_ALL.tsv');rows=[]
for ds,label in [('TJRCC','TJ-RCC'),('GSE175540','GSE175540')]:
 for endpoint,short in [('FAP_CAF25_disjoint__SPP1_MAC25_disjoint','FAP与SPP1状态程序'),('SPP1_MAC25_disjoint__ECM_org_heldout','SPP1与ECM程序')]:
  z=one(d,dataset=ds,stratum='overall',endpoint=endpoint,model='depth_adjusted')
  rows.append([label,short,str(int(z.n_patients)),f(z.median_effect),pv(z.p_raw),pv(z.FDR)])
table('表10 原文公开标记的肾癌程序迁移',['队列','关系','n','中位ρ','P','FDR'],rows,[2.6,4.3,1.0,2.1,3.0,3.0],source='去状态共享基因，并从ECM中移除状态定义基因；探索性迁移')
para('迁移分析支持程序联系，但原P7评分校准未通过，不能称原作者精确模型验证。TJ-RCC计算肿瘤样区CD3D平均表达仅6/10患者更低，P=0.753906；GSE175540为13/22，P=0.523467。TJ-RCC的CD3D检出比例反而9/10更高，差约+1.91个百分点，P=0.021484、FDR=0.085938。')
para('单基因和映射对照未重复程序关联：GSE175540的FAP—SPP1单基因调整ρ中位−0.055257，17/23负，P=0.034690、FDR=0.055504。区域内外UMI覆盖也不同。这里的整体计算区域与前述局部界面使用不同定义、终点和调整，因此保留各自结论，不能用某一版本替换其他空间证据。')

head('2.10 文献复算与真实扰动支持部分免疫机制并显示情境差异',2)
para('Qi等CRC研究的作者源数据复算显示，4位患者的FAP/SPP1联合区评分均正相关；40项GO原P复算一致；从公开镜像补齐P7后，4位患者的恶性区域CD3D/CD8A平均表达及检出比例方向均更低。P7论文图值与Pearson图注不一致、原免疫源表重复，以及确切评分签名/对照基因缺失均已记录。原版评分函数和过滤/归一化检查仍未恢复作者全部分值，故定位为关系与源结果复算。[2]')
para('肾癌离体SPP1阻断文献报告CD8中IFNγ、GZMB上升及部分功能障碍蛋白下降；其主要流式面板每组16例，干预转录组为7患者。逐患者流式源值和可靠处理映射未恢复，这部分仍是作者实验依据，不能列为本项目已独立复算的处理效应。[4]')
para('Trehan等研究的公开源值中，小鼠CD8耗竭评分均值由IgG的−1.842升至Spp1刺激的1.275，加入抗CD44后降至0.567，只抵消刺激增加量约22.7%。每组仅2个源评分值，独立供者对应不明，仅报告效应。全身Spp1缺失的肿瘤抗原阳性CD8中，CD39和TIM3比例分别下降18.82与13.03个百分点。[5]')
para('同研究15位HCC患者的配对空间源表中，αSMA富集区的巨噬—总CD8距离在14/15患者更短，距离比中位0.531，P=0.00005555、FDR=0.00049997；巨噬—CD39阳性CD8在15/15更短，距离比0.519，P=0.00042095、FDR=0.00189428。空间接近提供基质环境中的细胞关系，含义不同于免疫进入减少。[5]')
para('GSE274609的PDAC模型为每组3个独立混样池。宿主Spp1缺失后，成纤维黏着斑、收缩装置和肌成纤维程序的3个KO池均低于3个WT池；胶原成熟为7/9次跨组方向比较降低。成纤维捕获比例却增加3.362倍，CD8比例变化不一致，巨噬CXCL9/10程序降低。这里支持基质状态改变，不支持稳定的CD8招募恢复；3×3次比较不是9个独立重复。[6]')
g=data(319,'ALL_GROUP_SUMMARIES.tsv');m=data(319,'QUALIFIED_MODEL_CHECKS.tsv');rows=[]
for panel,endpoint,wt,ko,name in [('Ext Fig10H','Ext Fig10H','wt','Spp1-/-','重量归一化CAF计数'),('Fig7F','CD8+','WT','KO','重量归一化CD8数量'),('Ext Fig10C','Ext Fig10C','WT','KO','CD8中Ki67阳性比例（%）')]:
 a=one(g,panel=panel,endpoint=endpoint,group=wt);b=one(g,panel=panel,endpoint=endpoint,group=ko)
 mm=m[m.panel==panel];pp=pv(mm.iloc[0].P) if len(mm) else '未复算';ff=pv(mm.iloc[0].FDR_two_model_audit) if len(mm) else 'NA'
# PART2_CHUNK_1
 change=f'{(b["mean"]/a["mean"]-1)*100:+.1f}%' if panel!='Ext Fig10C' else f'{b["mean"]-a["mean"]:+.2f}个百分点'
 rows.append([name,f'{int(a.n_units)}只  {a["mean"]:.2f}',f'{int(b.n_units)}只  {b["mean"]:.2f}',change,pp,ff])
table('表11 CRC宿主Spp1干预源表复算',['端点','WT均值及n','KO均值及n','变化','P','本次FDR'],rows,[3.7,2.8,2.8,2.3,2.2,2.2],source='Henriques等；CAF按原线性模型，Ki67为动物随机截距重建；两模型复核家族')
para('CRC肝转移模型中，CAF数量降低44.6%，CD8数量增加53.0%，CD8增殖比例由30.563%升至50.249%。Ki67重建混合模型效应为+19.590个百分点，95%CI 10.650至28.530；它与表中原始动物组均值差是不同统计量。部分骨髓移植源表有重复动物编号和分组冲突，胶原匿名列在D14未降低，完整链条未作为端到端复现成功。[3]')
para('同一CRC研究的纯化T细胞实验中，CD3/28对照分裂读数均值72.833%，三个重组OPN剂量约75.067%、74.833%、75.833%，未见下降。该结果与其他刺激/阻断实验并列，说明不能把SPP1的效应概括为任何条件下直接抑制纯化T细胞增殖。')
para('空间CRISPR尝试中，纯T细胞Cd44测量资格不足；放宽为guide信号区域后，4张可评估切片的SPP1巨噬距离1增、2减、1平，独立动物对应未核实，不计算格点P。外源TYMP刺激的非肿瘤破骨分化背景GSE171542中，SPP1/C4程序反向下降，同样保留为情境不同的关键结果。[7]')

head('2.11 SPP1相关指标具有部分免疫治疗临床联系',2)
d=data(270,'ALL_SURVIVAL_MODELS.tsv');correction=data(271,'PHASE270_ADJUSTED_SEX_CORRECTION_AUDIT.tsv')
clinical=[]
for dataset,endpoint,score,arm,model,label in [
('CheckMate','OS','Jiang5','NIVOLUMAB','median_cohort','Nivo五基因 高低组OS'),
('CheckMate','OS','Jiang5','NIVOLUMAB','continuous','Nivo五基因 连续OS'),
('CheckMate','OS','Jiang5','NIVOLUMAB','continuous_adjusted','Nivo五基因 调整OS'),
('CheckMate','PFS','Jiang5','NIVOLUMAB','median_cohort','Nivo五基因 高低组PFS'),
('JAVELIN','PFS','SPP1','Avelumab+Axitinib','continuous_adjusted','联合治疗SPP1 调整PFS'),
('JAVELIN','PFS','Jiang5','Avelumab+Axitinib','continuous_adjusted','联合治疗五基因 调整PFS')]:
 z=one(d,dataset=dataset,endpoint=endpoint,score=score,arm=arm,model=model)
 if dataset=='CheckMate' and model=='continuous_adjusted':z=one(correction,endpoint=endpoint,score=score,arm=arm)
 clinical.append([label,f'{int(z.n)}/{int(z.events)}',f'{f(z.HR)} [{f(z.CI_low)}, {f(z.CI_high)}]',pv(z.P),pv(z.FDR)])
table('表12 免疫治疗队列的主要生存关联',['指标与模型','n及事件','HR及95%CI','P','FDR'],clinical,[5.2,2.0,4.6,2.1,2.1],source='Cox原家族；CheckMate调整模型使用性别编码修正版')
para('CheckMate的nivolumab组中，文献五基因评分与较短OS相关，连续未调整模型FDR显著；临床调整后HR=1.269884，P=0.016110、FDR=0.064440，为名义显著。固定高低组OS的log-rank P=0.01484，文献该结果及PFS阴性在原精度下复现。JAVELIN联合组的连续SPP1与较短PFS名义相关，但原文二分组效应未精确复现。[4,10]')
d=data(270,'TREATMENT_INTERACTIONS.tsv');z=d[(d.dataset=='CheckMate')&(d.endpoint=='OS')]
table('表13 CM025随机两臂范围内的OS治疗交互',['指标','n及事件','交互HR比值及95%CI','P','FDR'],[[x.score,f'{int(x.n)}/{int(x.events)}',f'{f(x.HR_ratio)} [{f(x.CI_low)}, {f(x.CI_high)}]',pv(x.P),pv(x.FDR)] for x in z.itertuples()],[3.5,2.1,5.6,2.4,2.4],source='目标检验为评分×治疗项；2指标家族')
para('交互HR比值描述两个治疗臂的评分斜率之比，不是患者接受nivolumab的治疗HR。该探索性信号指向评分较高时nivolumab相对everolimus的OS优势减弱。CM025的PFS交互及JAVELIN的PFS交互未获支持，因此不概括为跨方案通用耐药预测。')
head('2.12 C4临床转录状态呈时间依赖性而非统一耐药方向',2)
para('排除TYMP/SPP1及文献五基因后建立的C4状态评分，在留一研究内部验证中24/24位患者的C4分数高于同患者其他巨噬；两临床队列统一使用19个可测量基因。评分与SPP1单基因的相关仅约0.123和0.135，表示不同测量。')
tv=data(271,'C4_PFS_TIME_VARYING_MODELS.tsv');z=one(tv,arm='NIVOLUMAB',model='plus_Jiang5_macro')
para(f'在nivolumab完整调整模型中，C4对PFS的时间依赖性增量检验n={int(z.n)}、事件{int(z.events)}，联合P={pv(z.joint_P)}、FDR={pv(z.joint_FDR)}。该模型控制临床因素、文献五基因及总巨噬估计。')
hr=data(271,'C4_PFS_TIME_SPECIFIC_HR.tsv');hr=hr[(hr.arm=='NIVOLUMAB')&(hr.model=='plus_Jiang5_macro')&hr.month.isin([3,6,12,24])]
table('表14 同一时变模型的特定时点估计',['月数','HR','95%CI','在险人数'],[[str(int(x.month)),f(x.HR),f'{f(x.CI_low)}至{f(x.CI_high)}',str(int(x.at_risk))] for x in hr.itertuples()],[3.0,3.0,6.0,4.0],source='同一模型的点估计，不作多次独立验证')
para('C4独立OS、恒定或时变治疗获益交互及JAVELIN联合组的相同PFS联系均未获支持。Sunitinib完整调整PFS模型HR=0.799，P=0.004913、FDR=0.009826，较低风险方向也完整保留。时变系数描述基线RNA评分与风险集的关联随时间变化，不表示体内C4功能发生逆转。')

head('3 讨论')
head('3.1 最可保留的是来源明确的免疫转录联系',2)
para('综合结果中，最直接贴近TYMP和SPP1两个研究对象的依据来自巨噬细胞内。TYMP—SPP1表达联系在原集合和独立免疫富集队列中重复；两者的共同糖酵解联系与相反IFNγ方向可以同时存在。控制彼此并去除免疫程序共享基因后，IFNγ、抗原呈递与吞噬相关特征仍有患者一致性，支持区分两基因相关免疫状态。')
para('需要收窄的是效应范围。IFNα和IFNγ程序高度重叠，不能将其作为两份相互独立的证据；炎症程序还受成员和基因定义影响。五基因联合剔除使糖酵解效应减弱，尤其TYMP，但这种删除改变了程序本身，不能反推原关联是假象或证明被剔除基因具有因果驱动作用。')
head('3.2 不同层级的阳性与阴性回答不同问题',2)
para('组织bulk同时呈现细胞毒和功能障碍程序升高，描述的是组织构成及转录共变化；谱系限定的细胞内状态、同患者跨细胞联系、蛋白功能和静态空间位置则是不同估计对象。当前细胞来源限定的CD8/Treg/NK模型未形成统一抑制关系，这一结果不能否定巨噬细胞内联系，也不能被组织bulk阳性替代。')
# PART2_CHUNK_2
para('空间同样需要区分局部界面、整块计算区域和逐细胞癌周几何。某个区域内侧T/CD8较低、SPP1状态与ECM程序相关或基质丰富处CD8离癌细胞更远，并不自动证明SPP1介导物理阻挡。深度、局部密度、区域组成及边界定义影响这些估计；调整后变化保留为同一数据的不同条件问题，不选择最有利版本。')
head('3.3 外部干预补充机制依据但不替代本项目的来源链',2)
para('外部SPP1干预使免疫机制具有可核查的实验依据，其中CAF数量、CD8数量及增殖的部分源值已复算。不同癌种、宿主或骨髓来源干预、纯化细胞刺激和混样池实验并不等价；PDAC免疫端异质、纯化T细胞未见分裂下降及肾癌处理映射未解决的记录共同界定了适用范围。')
para('本项目的主要参与性阳性来自髓系来源TYMP。肿瘤来源TYMP到巨噬SPP1、接收端程序及完整免疫进入过程尚未获得同等支持；不能将这些分段相关拼接成已经证实的肿瘤TYMP经SPP1形成屏障。客户既有湿实验只作为独立背景，未被用作公共数据发现、筛选方向或复算效应的来源。')
head('3.4 临床延伸与统计边界',2)
para('SPP1文献指标的OS关联和CM025交互，以及C4状态的时间依赖性PFS关联，是不同指标、不同终点和不同模型的发现。当前没有完成可用于选药的阈值、预测准确性或临床效用验证，不能将其称为TYMP或C4耐药预测器。较晚时点在险人数减少、治疗与取材背景差异、跨方案不一致等限制集中保留。')
para('所有关联和定向复评均维持独立生物学单位。n不足时只描述效应，缺少完整标签或测量资格时记为不可评估。最新条件模型使用较小假设家族，FDR数字更低不代表生物证据自动增强；原发现家族和版本修正均在附表中原样保存。相同患者的多模型、基因剔除、技术分区及文献特征投射不作为新增独立验证。')

head('4 综合结论')
para('TYMP与SPP1相关巨噬状态具有可重复的表达联系，并呈现可区分的免疫转录特征。IFNγ、抗原呈递及吞噬相关方向通过了关键定向复评；糖酵解联系仍在，但对少数基因的联合剔除敏感。细胞特异接收端、空间进入和免疫治疗结果提供各自层面的补充，未建立统一的TYMP经SPP1导致免疫抑制或物理屏障的完整因果链。研究主线可采用“本项目TYMP参与性关联＋可复算的SPP1下游实验依据＋边界明确的免疫与临床延伸”进行收尾。')

head('参考资料')
refs=[
'[1] Liberzon等 MSigDB Hallmark基因集。基因集解释：https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp 。实际程序版本沿用冻结Hallmark2020及既有GO集合。',
'[2] Qi等 2022 Nature Communications 13 1742。Single-cell and spatial analysis reveal interaction of FAP+ fibroblasts and SPP1+ macrophages in colorectal cancer。https://doi.org/10.1038/s41467-022-29366-6',
'[3] Henriques等 2025 Nature Genetics 57 3050–3065。TGF-β builds a dual immune barrier in colorectal cancer by impairing T cell recruitment and instructing immunosuppressive SPP1+ macrophages。https://doi.org/10.1038/s41588-025-02380-2',
'[4] Jiang等 2025 Cancer Immunology Research。SPP1相关巨噬细胞及肾癌免疫治疗研究。https://pmc.ncbi.nlm.nih.gov/articles/PMC12485372/ ，DOI 10.1158/2326-6066.CIR-24-1146。',
'[5] Trehan等 2025 Nature Communications。SPP1与CD44功能及空间源数据。https://doi.org/10.1038/s41467-025-59529-0',
'[6] GSE274609 宿主Spp1缺失PDAC单细胞。https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274609 ；对应全文：https://pmc.ncbi.nlm.nih.gov/articles/PMC12479095/ 。',
'[7] SPAC-seq研究及作者资源。https://doi.org/10.1016/j.cell.2026.04.049 ；https://spac.pku-genomics.org/ ；https://github.com/zenglab-pku/SPAC-seq 。',
'[8] Obradovic肾癌单细胞原始数据与注释。https://data.mendeley.com/datasets/nc9bc8dn4m/1 ；https://github.com/Aleksobrad/single-cell-rcc-pipeline 。',
'[9] HCRN GU16-260数据及原文。https://datadryad.org/dataset/doi:10.5061/dryad.dv41ns28r ；https://pubmed.ncbi.nlm.nih.gov/39992403/ 。',
'[10] Xu等 2024 Human Vaccines & Immunotherapeutics。SPP1单基因与免疫治疗结局。https://pmc.ncbi.nlm.nih.gov/articles/PMC11093034/ 。',
'[11] Salmon等 基质与T细胞运动的实验依据。https://doi.org/10.1172/JCI45817',
'[12] 其余公共数据及技术来源：PRJNA1223824、GSE175540、GSE210041、GSE171542、GSE264586、GSE294109；各样本、来源文件及SHA256见配套来源清单。']
for x in refs:para(x)
head('附录A 完整统计与版本索引')
para('所有主要阴性、反向、资格不足及敏感性结果保留在配套统计附表包。P/FDR沿用每个源分析的原家族，未跨分析重新合并校正。报告中的简表用于阅读，完整系数、患者方向、原始P、FDR和检验家族以原表为准。')
table('附表索引',['内容','来源分析及目录','主要完整表'],[
['细胞身份与接收端','211、213、214','PRIMARY_12_HYPOTHESES；ALL_DIRECTED_HYPOTHESES；LITERATURE_SIGNATURE_ASSOCIATIONS'],
['组织bulk免疫','243','IMMUNE_ASSOCIATIONS_ALL及供者资格'],
['来源及候选输出','278、280、284、289、290、292、293、295、296','来源、过程、条件关联、接收网络及全部阴性'],
['独立队列与定向CD8','305、306、307','ALL_PATIENT_ASSOCIATIONS；RECIPIENT_ASSOCIATIONS_ALL；CD44伪bulk v2'],
['巨噬免疫程序及最终核查','314、321','ALL_EFFECT_SUMMARIES；FORMAL_SUMMARIES；EVIDENCE_RETENTION及所有删基因汇总'],
['额外CD8自由探索','320','ALL_ASSOCIATIONS及全部敏感性、留一患者结果'],
['空间分布与进入','222、234、246、275、291、298、300、309、311、312、313、317','各定义下切片、患者、模型、方向与资格'],
['文献复现与扰动','306、308、315、316、318、319','源值、复算效应、动物身份差异、不可完整复现记录'],
['免疫治疗及C4临床','270、271','生存与治疗交互；PHASE270_ADJUSTED_SEX_CORRECTION_AUDIT；时变模型']
],[3.0,4.2,8.8],source='附表文件逐项保留源路径及SHA256')
# PART2_CHUNK_3
para('版本优先规则：接收端身份采用211的freeze_v2；公开bulk采用243的44供者结果；CD44响应采用307的伪bulk v2；CheckMate受影响调整模型采用271的性别编码修正；原文P7采用316恢复数据；317为探索性迁移；314的发现统计与321的定向复评分别保留。旧失败、旧定义和重复拟合不被计算为独立证据。')
para('统计来源索引包含429份原始表或配套元数据及47条既有文字报告来源记录（44个文件）。它们按原目录结构保存；部分旧版本仅作追溯。STATISTICAL_SOURCE_MANIFEST与REPORT_SOURCE_MANIFEST提供来源、大小、SHA256及相对路径。')
for b in blocks:
 if b['type']=='table':
  b['rows']=[[{'literature5':'文献五基因','Jiang5':'文献五基因','cytotoxic':'细胞毒程序','dysfunction':'功能障碍程序'}.get(v,v) for v in row] for row in b['rows']]
for b in blocks:
 if b['type']=='table' and b['title'].startswith('表11 '):
  for row in b['rows']:
   for j in [4,5]:
    try:row[j]=format(float(row[j]),'.2e')
    except ValueError:pass
 if b['type']=='table' and b['title']=='附表索引':
  for row in b['rows']:row[1]=row[1].replace('、',', ')
(P/'inputs/DOCUMENT_CONTENT.json').write_text(json.dumps(blocks,ensure_ascii=False,indent=2))
(P/'inputs/DOCUMENT_USED_TABLES.json').write_text(json.dumps(list({x['source']:x for x in used}.values()),ensure_ascii=False,indent=2))
print('CONTENT_READY',len(blocks),'TABLES',sum(x['type']=='table' for x in blocks))

