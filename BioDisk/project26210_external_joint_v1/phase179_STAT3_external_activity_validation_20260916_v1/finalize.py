from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,shutil,datetime,csv
r=Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1");d=r/"phase179_STAT3_external_activity_validation_20260916_v1";v=d/"run_v2"
valid=json.loads((v/"VALIDATION_STATUS.json").read_text());assert valid["checks"]==valid["passed"]
sc=pd.read_csv(v/"SC_ACTIVITY_SPP1_SUMMARY.tsv",sep="\t");em=pd.read_csv(v/"BULK_EDGES_META.tsv",sep="\t");pm=pd.read_csv(v/"BULK_SHORTPATH_META.tsv",sep="\t");pc=pd.read_csv(v/"BULK_SHORTPATH_COHORT.tsv",sep="\t");c4=pd.read_csv(v/"C4_ACTIVITY_SUMMARY.tsv",sep="\t");lo=pd.read_csv(v/"SC_LEAVE_STUDY.tsv",sep="\t")
def table(f):
 def fmt(x):return f"{x:.6g}" if isinstance(x,(float,np.floating)) else str(x)
 return "| "+" | ".join(f.columns)+" |\n|"+"|".join(["---"]*len(f.columns))+"|\n"+"\n".join("| "+" | ".join(fmt(x) for x in row)+" |" for row in f.itertuples(index=False,name=None))
primarysc=sc[(sc.resource=="collectri")&(sc.scoring=="zmean")]
primaryedges=em[(em.resource=="collectri")&(em.scoring=="zmean")&(em.adjustment=="MTAP_myeloid")]
primarypath=pm[(pm.resource=="collectri")&(pm.scoring=="zmean")&(pm.adjustment=="MTAP_myeloid")]
primarycohort=pc[(pc.resource=="collectri")&(pc.scoring=="zmean")&(pc.adjustment=="MTAP_myeloid")]
report="""# MTA—STAT3—SPP1详细验证报告
最终统计版本：run_v2。原始run及独立验证失败记录保留，不能混用。
## 核心判断
**STAT3相关转录活性与巨噬细胞SPP1的联系获得了项目内支持；“MTA提高STAT3活性，再提高SPP1”的完整中介解释当前不成立。**
这一结果推进了SPP1接收端的机制候选定位，同时将缺口具体化到MTA端及细胞区室对应。MTA研究线路继续保留，不能把未获支持的统计短路径写成已确认中介，也无需据此删除外部受体依赖实验依据。

## 结果一：独立外部靶集支持整体巨噬STAT3—SPP1联系
主分析CollecTRI靶集共205个共同可测基因（180个激活方向、25个抑制方向），去除方向不明确边并核查SPP1、TYMP、STAT3均不在评分内。DoRothEA A/B/C资源对照62基因（57正向、5负向）。两资源存在证据重叠，不当作独立队列验证。

肿瘤内成熟巨噬31名有效患者，主评分深度调整后29/31正向，中位rho=0.08087，患者bootstrap中位数95%区间[0.04471,0.11129]，P=4.3194e-06。单一预设主要终点不作FDR。效应较小，但患者方向一致；属于转录活性代理与表达的相关证据。

主评分敏感性：
"""+table(primarysc[["adjustment","n","median_rho","positive","negative","P","median_CI_low","median_CI_high"]])+"""

控制C4组成后仍29/31正向；增加检出基因数后26/31正向；同分区拟时序调整后11/12正向。DoRothEA标准化均值评分27/31正向，rho=0.06396，P=4.0348e-05；CollecTRI秩均值评分29/31正向，rho=0.07403，P=4.7442e-06。以上为预设资源/评分敏感性，不以其最优结果替代主分析。

主评分逐研究留出：
"""+table(lo[(lo.resource=="collectri")&(lo.scoring=="zmean")&(lo.adjustment=="depth")][["omitted_study","n","median_rho","positive"]])+"""

各次留出中位方向仍正。此处为内部稳健性描述，不再增加留研究P值。

## 结果二：STAT3相关信号不等于C4特异富集
24名患者可比较C4与其他成熟巨噬。主CollecTRI标准化均值评分的C4差值中位数-0.00520，10/24正向、P=0.85267；当前不支持C4特异升高。
"""+table(c4)+"""
资源与评分敏感性中的部分C4结果为正，但主评分阴性且定义敏感，因此不写“C4由STAT3高活性驱动”。整体巨噬内部的相关性与C4群间差异回答不同问题；继续沿整体巨噬功能与C4定位的双层分工。

## 结果三：实测MTA端未连接到同向STAT3活性
四队列共255患者，主模型调整MTAP和髓系组成代理。采用与单细胞共同测得、方向一致的外部靶基因定义；各队列独立标准化，不将两类数据拼成同患者。
"""+table(primaryedges[["edge","n","positive_cohorts","rho","P","FDR_three_edges","CI_low","CI_high","I2"]])+"""
MTA—STAT3活性代理仅1/4队列正向，合并rho=-0.12332，P=0.10520、FDR=0.15780；STAT3代理—bulk SPP1条件关联rho=-0.01424，P/FDR=0.82858。当前数据不支持预期的两段正向连接。

增加缺氧调整后，MTA—代理4/4为负，代理—SPP1也4/4为负（后者rho=-0.25994、P=0.0002750、FDR=0.0008250），与巨噬细胞内的正向关系不同。该关键反向结果必须保留。两种负向路径乘积即使为正，也不能解释成“MTA升高STAT3，再升高SPP1”。

MTA—SPP1原始边际rho=0.13894、P=0.02926沿用同一历史255患者数据，不是新增独立重复。本轮三边家族FDR与历史候选家族FDR不同，统计表保留本轮口径，但不据较小校正池升级旧结果。主调整的总关联P=0.10448，不以总关联必须显著作为短路径计算门槛。

## 结果四：完整短路径未获支持
同队列内拟合 MTA→STAT3活性代理→SPP1 的关联分解：
"""+table(primarypath[["n","positive_cohorts","ab","P","CI_low","CI_high","bootstrap_low","bootstrap_high"]])+"""
主路径乘积a×b=0.003792，95%区间[-0.02131,0.02889]，P=0.76715；固定合并权重bootstrap区间[-0.01972,0.03190]同样跨零。全部12种资源×评分×调整组合的合并常规区间和bootstrap区间均跨零。

四队列主模型完整结果：
"""+table(primarycohort[["cohort","n","a_beta","b_beta","ab","bootstrap_low","bootstrap_high"]])+"""
本轮不报告“中介比例”，不将bulk短路径视为细胞间传递或因果中介。

## Methods
### 资源与评分
通过OmniPath下载人类STAT3的CollecTRI及DoRothEA记录，保存URL、原始表、边参考来源及冻结靶集SHA256。仅保留激活/抑制方向唯一且有参考文献的边；同基因重复边确认方向一致后去重。共同基因宇宙为冻结单细胞基因表和四个bulk表达表交集，使用已有单细胞输入过滤，不为结果选择基因。
主评分为log归一化表达的逐基因z分数与调控方向乘积的等权均值；单细胞在每患者整体成熟巨噬内标准化，bulk在队列内标准化。零方差基因不参与评分。CollecTRI最低20、DoRothEA最低10个有效基因；秩均值为敏感性。所有评分固定到12位小数后存储和排名，以保持并列秩可复算。
本评分命名为“外部STAT3有符号靶基因转录活性代理”，不是新SCENIC regulon，不是受体活性或pSTAT3实测，也不是直接套用某软件的标准ULM输出。

### 单细胞统计
固定肿瘤内成熟巨噬范围，排除单核样C1和补充Monocyte，不改变C4成员。每患者≥50细胞且SPP1≥10细胞检出。患者内计算原始和调整测序深度的秩相关；其他敏感性分别增加检出基因数、C4指示变量或同分区拟时序。后者要求≥50细胞、≥10个不同拟时序及≥10个SPP1检出细胞，各分区效应在患者内取中位。
患者层双侧Wilcoxon正态近似、连续性校正；主效应中位数区间来自5000次患者bootstrap。C4与其他巨噬均≥20细胞时另做患者配对评分差。细胞数只用于资格，不作独立推断样本量。

### 公共实测MTA统计
按冻结患者映射连接表达、MTA及协变量，四队列样本量97/52/50/56。变量队列内秩标准化，主模型调整MTAP与髓系代理，raw及增加缺氧为敏感性。三条相关边在每个队列/汇总及每项固定资源×评分×模型内作BH，短路径为预设单一主要终点。
路径a来自活动代理对MTA回归，b来自SPP1对代理及MTA的条件回归；每队列2000次患者bootstrap，重算秩但固定原队列评分校准。48个模型共96000次bootstrap。相关边采用Fisher-z随机效应合并，短路径用bootstrap标准误的随机效应合并；另报告固定合并权重bootstrap区间和留一队列结果。

## 独立验证与版本修正
第一次独立复算发现9项bulk秩评分路径因文本读回后的浮点并列秩差异不一致，最大差0.0006544。该原版及失败表完整保留。run_v2先固定评分精度后完整重跑；主CollecTRI z均值结果未改变，不放宽验证容差。
run_v2从稀疏原始计数逐基因复算两资源评分，以QR独立计算患者相关和bulk短路径，并检查摘要、BH和输入哈希：366/366通过，最大数值误差2.78e-16；另有480个抽样bootstrap逐次独立复算通过。合计846项技术核查，不称为846项独立生物学验证。

## Discussion：如何继续保留MTA逻辑
外部MTA实验显示STAT3功能依赖，但这不等价于MTA一定提高STAT3的转录活性代理；原研究也未稳定检出所测时点pSTAT3增加。由此推理，“A2B相关信号需要既有STAT3功能参与”仍是可检验的机制候选，与此次未获支持的单调上调中介模型不同。
这种解释目前是从分段实验推导的假说，不能用来回避公共MTA端的阴性或bulk反向结果。受体抑制/敲除、STAT3功能抑制以及MTA处理后的SPP1读出才直接对应功能依赖问题；客户尚未提供这些实验，不写成已完成。

靶基因代理也可反映其他调控、炎症或细胞状态，两个资源的方向一致并不保证STAT3特异性。bulk混合细胞来源、组织MTA与胞外暴露、不同评分校准层级均限制跨区室外推；本轮没有独立STAT3扰动队列来校准代理的生物学特异性。单细胞相关不证明STAT3驱动SPP1，bulk统计路径不识别因果。

## 后续取舍
保留：外部MTA—A2B/STAT3功能依赖依据、本轮整体巨噬STAT3代理—SPP1正向联系、既定C4表型。
不写：已证实MTA经STAT3活性上调驱动SPP1、C4特异STAT3激活、跨队列已完成中介。
下一轮如深化，优先获取直接髓系MTA/受体干预后的SPP1读出或合格独立扰动资源，检验功能依赖；不再用换靶集、换评分追求MTA→STAT3正相关。MTA继续保留，免疫和空间下游仍暂缓。

## 资源链接
[CollecTRI原始论文](https://academic.oup.com/nar/article/51/20/10934/7318114)；[官方资源仓库](https://github.com/saezlab/CollecTRI)；[OmniPath接口说明](https://r.omnipathdb.org/reference/collectri.html)。
机制来源及适用对象见Phase178报告；本轮新增证据为项目内统计对应，不将既有文献改称新实验。
"""
(d/"STAT3_DETAILED_VALIDATION_CN.md").write_text(report,encoding="utf8")
(d/"STAT3_MANUSCRIPT_MODULE_CN.md").write_text("""# 可拆入统一正文的STAT3候选调控模块
## 外部靶基因程序支持整体巨噬细胞中STAT3相关转录状态与SPP1的联系
基于独立资源确定STAT3有符号靶基因，剔除SPP1、TYMP及STAT3本身后构建转录活性代理。31名有效肿瘤患者中29名呈正向关联，深度调整后相关系数中位数0.0809（P=4.32×10⁻⁶）；控制检出基因数、C4组成及拟时序后方向保留，另一外部靶集与秩评分支持相同方向。该结果为STAT3相关转录状态参与SPP1髓系表型提供关联依据。

四个实测MTA队列的255例患者未支持MTA与该活性代理的稳定正向联系，MTA—代理—SPP1短路径区间跨零。因此，STAT3保留为SPP1相关调控候选；MTA是否通过A2B相关信号依赖STAT3功能影响SPP1仍作为待检验机制，不写成已经确认的活性上调中介。

定位补充：主评分不支持C4相对其他巨噬普遍升高，故上述结果归于整体巨噬状态分析，不将C4定义为STAT3高活性群体。完整方法、统计及关键反向结果见详细验证报告。
""",encoding="utf8")
appendix="# 完整汇总统计（最终run_v2）\n\n"
for name,title in [("SC_ACTIVITY_SPP1_SUMMARY.tsv","单细胞全部资源/评分/模型"),("C4_ACTIVITY_SUMMARY.tsv","C4定位"),("BULK_EDGES_META.tsv","公共队列三边合并"),("BULK_SHORTPATH_META.tsv","全部短路径合并"),("BULK_SHORTPATH_COHORT.tsv","全部队列短路径")]:
 appendix+="## "+title+"\n\n"+table(pd.read_csv(v/name,sep="\t"))+"\n\n"
(d/"ALL_STATISTICS_CN.md").write_text(appendix,encoding="utf8")
stamp=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
status={"status":"COMPLETE_VALIDATED","authoritative_statistics":"run_v2","primary_targets":205,"sensitivity_targets":62,"sc_n":31,"bulk_n":255,"bootstrap_draws":96000,"numeric_checks":846,"new_figures":0,"STAT3_SPP1_association_supported":True,"MTA_STAT3_SPP1_mediation_supported":False,"updated_at":stamp}
(d/"STATUS.json").write_text(json.dumps(status,indent=2))
bk=d/"root_record_backups";bk.mkdir(exist_ok=True)
entry=f"""
## Phase179 STAT3外部活性详细验证完成（{stamp}）
用户要求详细验证MTA—STAT3—SPP1。冻结外部CollecTRI205靶(180+/25-)及DoRothEA62靶，SPP1/TYMP/STAT3均不参与评分。最终统计run_v2。
整体肿瘤巨噬31患者29正，主深度rho中位0.08087/P4.319e-6；基因数调整26/31正，C4组成29/31正，拟时序11/12正；资源/秩评分及留研究支持方向。主C4群间评分不升高(P0.85267)，不称C4特异激活。
公共255例主MTA—代理rho-0.12332/P0.10520/FDR0.15780；代理—SPP1条件rho-0.01424/P=FDR0.82858；短路径ab0.003792/CI[-0.02131,0.02889]/P0.76715，bootstrap亦跨零。全部12汇总短路径不支持。缺氧调整下两段负向不可包装正向上调中介。
366独立计数/QR/汇总/哈希+480bootstrap复算全通过，96000bootstrap，无图。原run有9个浮点并列秩复算差异，保留失败记录；run_v2评分12位小数后重跑，主结论不变。
报告：{d}/STAT3_DETAILED_VALIDATION_CN.md
正文：{d}/STAT3_MANUSCRIPT_MODULE_CN.md
保留STAT3-SPP1关联与MTA功能依赖候选；不确认MTA提高STAT3活性的中介。继续MTA研究时优先直接髓系MTA/受体扰动读出，不再换评分追正向；免疫空间仍暂缓。
"""
for name in ["PROJECT_STATUS.md","DECISIONS.md","CODEX_HANDOFF.md"]:
 p=r/name
 if not (bk/name).exists():shutil.copy2(p,bk/name)
 old=p.read_text()
 if "## Phase179 STAT3外部活性详细验证完成" not in old:p.write_text(entry+"\n"+old,encoding="utf8")
p=r/"RESULTS_MANIFEST.tsv"
if not (bk/p.name).exists():shutil.copy2(p,bk/p.name)
if "phase179_STAT3_external_activity" not in p.read_text():
 with p.open("a",newline="") as f:csv.writer(f,delimiter="\t").writerow(["phase179_STAT3_external_activity","L2","report",str(d/"STAT3_DETAILED_VALIDATION_CN.md"),stamp,"PASS","run_v2:846 checks;96000 bootstrap","SC association supported;bulk short path unsupported"])
de=d/"delivery";de.mkdir(exist_ok=True)
for n in ["STAT3_DETAILED_VALIDATION_CN.md","STAT3_MANUSCRIPT_MODULE_CN.md","ALL_STATISTICS_CN.md","STATUS.json"]:shutil.copy2(d/n,de/n)
for n in ["SC_ACTIVITY_SPP1_SUMMARY.tsv","BULK_EDGES_META.tsv","BULK_SHORTPATH_META.tsv","FROZEN_STAT3_TARGETS.tsv","VALIDATION_STATUS.json"]:shutil.copy2(v/n,de/n)
hh={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in de.iterdir()};(d/"DELIVERY_SHA256.json").write_text(json.dumps(hh,indent=2));print(json.dumps({"status":status,"hashes":hh},ensure_ascii=False))
