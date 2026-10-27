from common import *
import json,shutil,hashlib,datetime
read=lambda p:pd.read_csv(D/p,sep='\t')
validation=json.loads((D/'qc/VALIDATION_STATUS.json').read_text());assert validation['status']=='PASS'
assert json.loads((D/'qc/ATTRIBUTION_STATUS.json').read_text())['status']=='PASS'
def fmt(v):
 if pd.isna(v):return '不检验/不可估计'
 if isinstance(v,(float,np.floating)):return f'{v:.4g}'
 return str(v)
def table(df,cols=None,rename=None):
 if cols:df=df[cols]
 if rename:df=df.rename(columns=rename)
 lines=['| '+' | '.join(map(str,df.columns))+' |','|'+'|'.join(['---']*len(df.columns))+'|']
 for row in df.itertuples(index=False,name=None):lines.append('| '+' | '.join(fmt(v).replace('|','/') for v in row)+' |')
 return '\n'.join(lines)
s=json.loads((D/'qc/RESULTS_FOR_REVIEW.json').read_text())
reg=read('core/CONTRAST_REGISTRY.tsv');markers=read('core/KEY_MARKER_COMPARISONS.tsv')
fun=read('functional/ALL_PROGRAM_RESULTS.tsv.gz');paired=read('functional/U_MINUS_C4_SHARED_PATIENT_PROGRAMS.tsv.gz')
bridges=read('bridge/ALL_FIXED_CANDIDATE_BRIDGES.tsv');receiver=read('receiver/WITHIN_CELL_SUMMARY.tsv')
rhof=read('bridge/U_MINUS_C4_RHO_BOOTSTRAP.tsv')
stat=read('receiver/STAT3_ACTIVITY_SUMMARY.tsv');statdiff=read('receiver/STAT3_SHARED_PATIENT_CHANGE.tsv')
attrib=read('receiver/STAT3_ATTRIBUTION_SUMMARY.tsv')
terms=['Cholesterol metabolism','Lysosome','Phagosome','Antigen processing and presentation','Glycolysis','E2F Targets','G2-M Checkpoint']
mod=[
('群体/marker/SPP1来源','原始计数重新汇总；联合、C4、C10、C15及固定剔除敏感性','完成','core/'),
('功能评分与ORA','3套冻结基因集；完整及去除TYMP/SPP1；同患者、共同基因和共同背景','完成','functional/'),
('KEGG GSEA','联合/C4、共同患者、去C10主导患者；各自完整KEGG家族BH','完成','gsea/'),
('pySCENIC与三TF','原冻结网络；更新群间对比、轨迹相关、TYMP/SPP1联系及靶基因重叠','完成','regulons/'),
('原图轨迹/动态','联合及C10/C15在原分区图上的覆盖、排序及动态统计','完成','trajectory/'),
('TYMP/候选X桥接','同16文库；原116候选、38代谢成分、4个MTA相关代理；共632个比例终点模型','完成','bridge/'),
('群内SPP1表达终点','各自合格与共同12患者；全部原候选家族，2784个模型','完成','bridge/MATCHED_GROUP_SPP1_FIXED_CANDIDATES.tsv.gz'),
('MTA/受体/接收端','MTA响应代理、受体、SLC16A3、IL6/IL10/LGALS3、STAT3；更新组成调整','完成','receiver/'),
('外部实测MTA及独立文献','没有C4标签的组织/细胞系/扰动数据，保留原结果和定义','不受本次分群改变影响','原Phase157/160/162/163/174/175/176/178/179/180/183–187'),
('旧空间与免疫模型','原SPP1-state不是C4，更不是三群联合；不把旧状态直接改名','原结果保留；无既有联合群空间映射','原Phase164–170及下游资格记录'),
('第一部分及客户湿实验','无C4依赖的既有分析及客户数据','保持原结果；未进入新模型','原项目报告与来源记录')]
pd.DataFrame(mod,columns=['模块','处理内容','状态','结果位置']).to_csv(D/'SCOPE_AND_DEPENDENCIES.tsv',sep='\t',index=False)
# Concise client-readable outcome summary.
summary_text="""# C4+C10+C15联合复评：覆盖更完整，并出现新的功能线索，但不是所有关联都增强

## 结论
保留三群联合研究是有价值的。本轮显示：联合群扩大了SPP1来源覆盖，揭示了C15参与的增殖、糖酵解和STAT3转录活性代理信号；原C4的脂质/溶酶体特征仍在，但被混合状态稀释。TYMP直接关联、MTA响应代理及轨迹覆盖没有同步改善。

## 最重要的对照
| 项目 | 原C4 | C4+C10+C15 | 判断 |
|---|---:|---:|---|
| 冻结细胞数 | 2,674 | 4,502 | 覆盖扩大，三个原簇身份保留 |
| all-tissue有效配对患者 | 26 | 28 | 增加2名；比较优劣另用共同26名 |
| 同16患者的细胞占比中位数 | 7.41% | 9.76% | 增加 |
| 同16患者的巨噬SPP1原始RNA贡献中位数 | 54.59% | 59.63% | 增加约5.04个百分点 |
| 糖酵解程序配对中位差 | -0.133 | +0.123 | 由负转正；联合FDR=0.00543 |
| 恶性TYMP—群体比例rho（同16患者） | 0.468 | 0.291 | 点估计减弱；联合P=0.274，FDR=0.411 |
| 原主轨迹有效覆盖 | 442/2674，16.53% | 463/4502，10.28% | 有效细胞仅增加21，整体覆盖率下降 |

SPP1在联合群中仍富集：28名患者25名正向，配对log2(CPM+1)差中位1.821，P=2.93e-5，FDR=4.17e-4。
但在共同26名患者中，逐患者比较“联合对比效应减C4对比效应”，仅10名增加、16名减少，中位增量-0.052，P=0.980；因此不能把两个总体中位数的差直接说成每位患者均增强。

## 新增的可用线索：C15参与的糖酵解、增殖及STAT3状态
联合群的G2-M/E2F程序变为正向。去掉C10后方向保留，去掉C15后糖酵解和G2-M回到负向，说明新增功能主要与纳入C15有关，不是原C4的所有细胞都变强。

CollecTRI STAT3转录活性代理由C4的无明确群间升高，变为联合群26名患者中20名升高：中位差+0.0231，P=0.0111，原主比较家族FDR=0.0333。
同24名患者、同一个non-U参照背景下，联合减C4的STAT3代理中位增量+0.0186，24/24正向，bootstrap95%区间0.0109–0.0290。去掉C15后不再见这一整体升高，去掉C10后仍保留方向。该结果支持保留C15作独立解释。

恶性细胞糖酵解程序与联合群比例的原始rho从0.565升至0.682；同患者分层bootstrap的rho差95%区间0.0089–0.275。研究调整后的关联仅为名义显著，原50程序家族FDR未通过，保留为探索性线索。

## 未被合并改善的关键结果
- TYMP—群体比例没有增强；研究调整后联合rho=0.107，置换P=0.725。
- 原公共成纤维MTA响应转录代理在联合群中仍以负向差异为主；肿瘤26患者18名负向。它没有因合并变成正向MTA响应群。
- 脂质和溶酶体功能仍有支持，但GPNMB/CCL18等原C4特征减弱；吞噬体全程序的患者级增强仍不支持。
- 主轨迹对联合群覆盖不足；排除应激/增殖的旧敏感性模型本身没有C15，不能当作完整三群联合的验证。
- 现有外部组织实测MTA、旧SPP1空间状态不能被直接重新命名为联合群证据。

## 建议采用的研究分工
联合群用于总体SPP1来源、丰度及上游联系；三群分别解释状态：C4偏脂质/溶酶体，C15偏增殖并贡献新增STAT3/糖酵解线索，C10保留为高度患者集中的SPP1高表达候选。
因此适合“双层保留”，不以联合群取代C4，也不再把C4当作唯一SPP1研究对象。

本轮为统计复评，没有改图或插入正文。完整方法、正负结果、敏感性和旧分析依赖关系见《完整复评报告.md》及统计表；项目原始数据和旧结果均保留。
"""
(D/'SUMMARY_CN.md').write_text(summary_text,encoding='utf8')
# Complete numerical report, including opposing and non-estimable findings.
parts=["# C4∪C10∪C15联合候选群的系统复评",summary_text.split("## 结论")[1].split("## 最重要")[0].strip(),
"## 数据与方法",
"分析沿用75文库项目的既有成熟巨噬细胞集合：22,444细胞、39名患者。联合群固定为原C4/C10/C15的集合，不重新聚类。全组织marker/功能主比较保持原分析人群；接收端使用肿瘤组织，同文库上游关联使用原16名恶性细胞与巨噬同时合格患者。原始细胞计数按患者、文库、组织及原簇汇总；全部旧患者×簇计数逐项复现。",
"表达比较要求同一患者目标群与参照群各≥20细胞。全基因表达为原始counts汇总后log2(CPM+1)，先计算每位患者的配对差再汇总。分别提供各定义自身的合格人群、共同患者，以及以非联合群为共同背景的比较。功能的直接U−C4比较另外使用两边共同可测基因。",
"基因过滤沿用原count覆盖条件；配对统计为双侧Wilcoxon近似检验、连续性校正，少于3名患者不作推断。3–9名患者按探索性呈现。FDR分别在每个对比的基因家族、基因集库/评分版本、regulon家族、候选类别/终点/边/调整方式内计算。",
"功能库沿用Hallmark_2020、KEGG_2021_Human及Reactome_2022；完整与去除TYMP/SPP1两版独立报告。程序效应为患者内成员基因log表达差的平均值，gene-set交集10–500；ORA沿用名义P<0.05且绝对中位差≥0.25的上/下基因及原测量背景。KEGG GSEA使用患者配对t统计量排序，完整合格KEGG家族BH，原对象和完整结果保留。",
"同16患者TYMP关系采用Spearman及研究内秩残差相关/9,999次研究内置换。候选X沿用原116项初筛、38项代谢成分及4个MTA相关代理；标准化秩关联分解与2,000次分层bootstrap沿用原方法。原始rho差的2,000次bootstrap在每次抽样后重新排名。群内表达终点另分析各自合格与共同12名患者，保留全部候选家族。",
"调控使用已冻结pySCENIC新网络与逐细胞AUCell，以及既有STAT3目标集/评分缓存。细胞内相关先在患者内计算效应，达到≥50细胞和相关RNA各≥10检出细胞时纳入；随后以患者效应统计，报告深度、检出基因数、缺氧及群体组成调整。轨迹沿用原分区图、根与伪时间，分区效应在患者内归并；动态和调控轨迹均更新联合群相关分区的资格。",
"固定敏感性包括去掉C10主导患者、去研究、去C10/去C15成分、共同患者、共同背景、功能去TYMP/SPP1，以及原无应激/增殖轨迹。源文件、脚本、失败日志、版本和校验均在服务器本阶段目录。",
"## 联合群增加SPP1来源覆盖，患者内富集增量并不一致",
table(read('core/COMPOSITION.tsv'),rename={'group':'群体','cells':'细胞数','patients':'出现患者数','studies':'出现研究数','dominant_patient_fraction':'最大单患者细胞占比','dominant_study_fraction':'最大单研究细胞占比'}),
table(pd.DataFrame(s['SPP1']),rename={'contrast':'比较','n':'患者n','median_delta':'SPP1配对中位差','positive':'正向患者','P':'P','FDR':'FDR'}),
"同16患者的细胞占比和原始SPP1 RNA份额中位数如下。各群中位数不可直接相加；联合份额由每位患者先合计原始counts再求中位数。",
table(read('bridge/SAME_LIBRARY_GROUP_METRICS.tsv').groupby('group')[['fraction','SPP1_RNA_share']].median().reset_index()),
"两个群体各自的效应中位数之差不等于逐患者增量的中位数。下面直接配对检验，common为各自原参照但相同患者，shared_background为两者都对非联合群。",
table(read('core/DIRECT_U_MINUS_C4_MARKER_EFFECTS.tsv').query("gene in ['SPP1','GPNMB','CCL18']"),['scope','reference','gene','n','U_median','C4_median','median','positive','negative','P','FDR','CI_low','CI_high']),
"## 脂质与溶酶体特征保留，新增糖酵解及增殖轴主要由C15参与",
table(fun[(fun.contrast.isin(['all_U','all_C4']))&(fun.version=='without_SPP1_TYMP')&fun.term.isin(terms)],['contrast','term','n','median','positive','P','FDR']),
"全部共同患者、共同背景、共同基因集合的直接功能增量见统计表；以下为核心四个KEGG程序和糖酵解。",
table(paired[(paired.scope=='all')&(paired.reference=='shared_background')&(paired.version=='without_SPP1_TYMP')&paired.term.isin(terms)],['term','identical_genes','n','U_median','C4_median','median','positive','P','FDR']),
"预先固定的去成分诊断如下。去C15与去C10是归因诊断，不用于从多个组合中挑选显著性最优定义。",
table(read('functional/CHANGED_AXIS_COMPONENT_ATTRIBUTION.tsv'),['contrast','term','n','median','P','FDR']),
"## 排序富集与患者整体程序评价分开呈现",
table(read('gsea/FOCUS_GSEA_RESULTS.tsv').query("model in ['all_U_common','all_C4_common'] and version == 'complete'"),['model','ID','patients','intersection_size','NES','pvalue','p.adjust']),
"吞噬体的GSEA与患者整体程序不是同一终点：前者评价基因在排序中的集中，后者评价全部成员基因的患者配对平均变化。保留完整二者，不用前者替代后者的阴性。",
"## STAT3转录活性代理出现新的联合群差异，C15是主要成分",
table(stat,['resource','group','n','median','positive','P','FDR']),
"同一患者、同一参照下的直接变化如下；nonU比较排除了参照组变动的影响。",
table(statdiff),
table(attrib.query("resource == 'collectri' and reference == 'own_complement'"),['group','population','n','median','positive','P','FDR']),
"## TYMP直接联系没有增强，糖酵解与部分代谢候选保留探索性线索",
table(pd.DataFrame(s['TYMP_direct'])),
table(rhof),
"研究调整后的预设重点候选列于下表。糖酵解及SLC16A3部分边为名义显著，但原候选家族FDR未过；关联分解区间跨零时不称中介成立。",
table(bridges[bridges.group.eq('U')&bridges.adjustment.eq('study_adjusted')&bridges.candidate.isin(['ligand:IGFBP3','program:Cholesterol Homeostasis','program:Glycolysis','gene:SLC16A3','gene:LDHA','direct_MTA_response_RNA_proxy','polyamine_synthesis_core_RNA','MTAP_clearance_RNA'])],['candidate','n_patients','a_rho','a_P','a_FDR','XY_rho','XY_P','XY_FDR','b_partial_rho','b_P','b_FDR','ab_low','ab_high']),
"群内SPP1表达/检出/贡献使用共同12患者直接对照：",
table(read('bridge/MATCHED_GROUP_SPP1_DIRECT_TYMP.tsv').query("population == 'common_eligible'")),
"## 接收端关联并非全面增强，原MTA响应代理仍非正向富集",
table(read('receiver/MTA_RESPONSE_SUMMARY.tsv')),
"以下为患者内相关后按患者汇总的深度调整结果。C10合格患者不足时保留资格表，不用其细胞数替代样本数。n<10的结果按探索性阅读。",
table(receiver[receiver.group.isin(['U','C4','C10','C15'])&receiver.adjustment.eq('depth')&receiver.module.isin(['same_library_lactate','STAT3_frozen_collectri','RNA'])],['module','group','edge','n','median','positive','P','FDR']),
"主要受体的RNA检测及配对表达情况如下：",
table(read('receiver/RECEPTOR_GENE_PANELS.tsv').query("group in ['U','C4'] and family == 'MTA_receptor'"),['group','gene','n','target_median_detection','target_patients_ge10pct','median','positive','P','FDR']),
"原116候选分泌配体对应受体的各亚基可测性完整保存于服务器receiver/FIXED_LIGAND_RECEPTOR_AVAILABILITY.tsv；整体巨噬的C4组成调整也已与联合组成调整并列计算。",
"## 原三TF结论按新定义复核，未预设统一上调",
table(read('regulons/ALL_GROUP_AUC_SUMMARY.tsv').query("group in ['U','C4'] and scope == 'tumor' and TF in ['JUN','JUNB','CEBPD']"),['group','TF','n','median','positive','P','FDR']),
"逐患者regulon、TYMP/SPP1 RNA相关、目标重叠与组成调整均在完整统计目录。调控目标包含结局基因的情况单列审计，不把分数与其成员基因相关当作独立机制验证。",
table(read('regulons/THREE_TF_OUTCOME_MEMBERSHIP_AUDIT.tsv')),
"## 联合群在原轨迹中的有效覆盖率下降",
table(read('trajectory/COVERAGE_SUMMARY.tsv').query("pseudotime == 'pt_primary'")),
"联合群主分析有效细胞463，比C4的442多21；但原始联合群有4502细胞，因此覆盖率降为10.28%。无应激/增殖模型不含C15：其联合剩余成员有效率43.95%，相对原4502细胞仅38.94%，不能解释成完整三群联合65%以上覆盖。",
"原图上的群间排序、表达动态和regulon趋势完整保存；不强连分区，不据少量有伪时间细胞推断整个联合群的终末走向。",
"## 历史分析逐项依赖处理",
table(pd.DataFrame(mod,columns=['模块','处理','状态','位置'])),
"## 讨论与研究分工",
"联合群是根据已有SPP1分布形成的事后组合，不是一个重新发现的均一亚型；SPP1富集及同数据差异属于条件性证据。C10的97.64%细胞来自同一患者，其all-tissue配对n=2，肿瘤接收端部分分析仅n=1；C15的配对n=9或更少，保持探索性定位。患者内对照和固定剔除检查用于识别这些组成影响。",
"STAT3为转录活性代理，程序分数/GSEA为转录层结果；不直接等同蛋白磷酸化、代谢通量或细胞功能。原MTA响应签名来自公共成纤维扰动，不能等同肾癌巨噬实测MTA效应。当前没有同患者实测MTA与三群联合的对应资料，不能由合并后的STAT3阳性确认MTA→STAT3→SPP1链条。临床及免疫/空间后果继续按独立终点处理。",
"建议将联合群用于总体SPP1来源与上游关联，保留C4脂质/溶酶体、C15增殖及新增STAT3/糖酵解、C10患者集中候选三个分支。主要阴性/反向结果仍进入总判断；不因联合分析产生阳性而替换所有原C4结论。",
"## 校验与文件",
f"原始counts与旧患者×簇矩阵逐项一致；原C4 marker、功能和regulon基线复现。独立QR、BH、GSEA运行分数及原始SPP1汇总核验共{validation['checks']}项通过。STAT3共享背景及去成分检查另列。无新图，既有报告和图件未改。",
f"权威目录：{D}。所有完整原始统计、脚本、R对象、失败日志和敏感性均在此目录；本地交付包含可读说明与聚合统计表。"]
(D/'FULL_REASSESSMENT_CN.md').write_text('\n\n'.join(parts),encoding='utf8')
# Explicit new-authorized scope recorded without changing the original C4 membership.
agreement="""
## 2026-09-18 用户授权C4+C10+C15联合及拆分复评
用户同意保留C4、C10、C15，并授权从联合群角度重新评估此前依赖C4定义的分析。固定U=C4∪C10∪C15原成员，不重新分群；原C4及三群拆分均保留，不以显著性最优组合替代预先定义联合群。采用患者单位及相同患者/背景比较，保留增强、减弱、阴性和研究集中结果。与C4定义无关的外部实测MTA、旧空间状态及第一部分不因本次任务改名或重算。当前先交付文字统计，不自动改现有三张主图。
"""
for p in [R/'AGENTS.md',R/'part2_rebuild_v1/AGENTS.md']:
 backup=D/'record_backup'/('root_AGENTS_before_union.md' if p.parent==R else 'part2_AGENTS_before_union.md')
 if not backup.exists():shutil.copy2(p,backup)
 if '用户授权C4+C10+C15联合及拆分复评' not in p.read_text():
  with p.open('a',encoding='utf8') as f:f.write(agreement)
# All source-derived attribution summaries reproduce the primary values.
a=attrib.query("reference == 'own_complement' and population == 'own' and group in ['U','C4']")
j=a.merge(stat,on=['resource','group'])
assert np.allclose(j['median_x'],j['median_y']) and np.allclose(j.P_x,j.P_y)
for keys,z in attrib.groupby(['resource','reference','population']):assert np.allclose(bh(z.P),z.FDR,equal_nan=True)
for keys,z in statdiff.groupby('reference'):assert np.allclose(bh(z.P),z.FDR,equal_nan=True)
(D/'qc/ATTRIBUTION_VALIDATION.json').write_text(json.dumps({'status':'PASS','primary_effects_and_P_reproduced':True,'attribution_BH_verified':True,'shared_background_BH_verified':True},indent=2))
delivery=D/'delivery';(delivery/'统计表').mkdir(exist_ok=True)
for src,dst in [(D/'SUMMARY_CN.md',delivery/'联合群复评_结论摘要.md'),(D/'FULL_REASSESSMENT_CN.md',delivery/'完整复评报告.md'),(D/'SCOPE_AND_DEPENDENCIES.tsv',delivery/'统计表/分析范围与依赖.tsv')]:shutil.copy2(src,dst)
files={
'01_群体与SPP1':D/'core/KEY_MARKER_COMPARISONS.tsv',
'02_同患者marker增量':D/'core/DIRECT_U_MINUS_C4_MARKER_EFFECTS.tsv',
'03_完整功能评分':D/'functional/ALL_PROGRAM_RESULTS.tsv.gz',
'04_共同基因功能增量':D/'functional/U_MINUS_C4_SHARED_PATIENT_PROGRAMS.tsv.gz',
'05_完整GSEA':D/'gsea/ALL_GSEA_RESULTS.tsv',
'06_TYMP直接关联':D/'bridge/TYMP_DIRECT_RESULTS.tsv',
'07_固定候选桥接':D/'bridge/ALL_FIXED_CANDIDATE_BRIDGES.tsv',
'08_群内SPP1候选关联':D/'bridge/MATCHED_GROUP_SPP1_FIXED_CANDIDATES.tsv.gz',
'09_接收端关联':D/'receiver/WITHIN_CELL_SUMMARY.tsv',
'10_STAT3成分归因':D/'receiver/STAT3_ATTRIBUTION_SUMMARY.tsv',
'11_STAT3同患者增量':D/'receiver/STAT3_SHARED_PATIENT_CHANGE.tsv',
'12_全部regulon':D/'regulons/ALL_GROUP_AUC_SUMMARY.tsv',
'13_轨迹覆盖':D/'trajectory/COVERAGE_SUMMARY.tsv',
'14_候选相关系数差':D/'bridge/U_MINUS_C4_RHO_BOOTSTRAP.tsv',
'15_原始同文库SPP1来源':D/'bridge/SAME_LIBRARY_GROUP_METRICS.tsv',
'16_MTA响应代理':D/'receiver/MTA_RESPONSE_SUMMARY.tsv',
'17_受体表达':D/'receiver/RECEPTOR_GENE_PANELS.tsv',
'18_增强功能成分归因':D/'functional/CHANGED_AXIS_COMPONENT_ATTRIBUTION.tsv'}
for name,p in files.items():shutil.copy2(p,delivery/'统计表'/(name+('.tsv.gz' if p.name.endswith('.gz') else '.tsv')))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
files_to_hash=list((D/'scripts').glob('*'))+[D/'PLAN_BEFORE_RESULTS.md',D/'inputs/FROZEN_UNION.json']
(D/'qc/CODE_PLAN_SHA256.json').write_text(json.dumps({str(p):sha(p) for p in files_to_hash if p.is_file()},indent=2))
manifest={str(p.relative_to(delivery)):sha(p) for p in delivery.rglob('*') if p.is_file() and p.name!='SHA256.json'}
(delivery/'SHA256.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
(D/'STATUS.json').write_text(json.dumps({'status':'VALIDATED_READY_FOR_DELIVERY','union_cells':4502,'all_tissue_paired_U':28,'all_tissue_paired_C4':26,'same_library_patients':16,'original_results_preserved':True,'validation_checks':validation['checks'],'new_figures':False,'report_inserted':False,'delivery_files':len(manifest)},indent=2))
print('REPORT_AND_PACKAGE_READY',len(manifest),flush=True)
