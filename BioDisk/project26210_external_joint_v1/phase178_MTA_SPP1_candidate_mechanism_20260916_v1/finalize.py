from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,shutil,csv,datetime,urllib.request,concurrent.futures
r=Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1");d=r/"phase178_MTA_SPP1_candidate_mechanism_20260916_v1"
# Preserve failed download history; retry without confusing retrieval with source review.
urls={"A2B_OPN_2010.xml":"https://www.ebi.ac.uk/europepmc/webservices/rest/PMC2797041/fullTextXML","STAT3_SPP1_2025.html":"https://link.springer.com/article/10.1186/s13046-025-03303-z"}
def get(kv):
 k,u=kv
 try:
  b=urllib.request.urlopen(u,timeout=15).read();(d/"sources"/k).write_bytes(b)
  return {"file":k,"url":u,"status":"retrieved","bytes":len(b),"sha256":hashlib.sha256(b).hexdigest()}
 except Exception as e:return {"file":k,"url":u,"status":"failed","error":str(e)}
with concurrent.futures.ThreadPoolExecutor(2) as ex:retrieval=list(ex.map(get,urls.items()))
(d/"sources/retrieval_retry.json").write_text(json.dumps(retrieval,indent=2))
report="""# MTA→SPP1继续深化：A2B受体与STAT3候选机制
## 结论：可增加有原始实验依据的中间机制，MTA线路继续保留
最值得优先深化的是“MTA—A2B依赖的髓系信号—SPP1”，STAT3作为候选下游调控环节。现在找到的是不同模型中的机制分段依据，不是本项目已经验证的完整串联中介。免疫与空间分支暂缓，继续优先补MTA接收端。

不能要求MTA提高受体RNA或STAT3 RNA才算通路有作用：受体参与和转录因子功能依赖，与其转录量变化是不同问题。此前C4中STAT3 RNA较低不能单独排除STAT3依赖机制；低ADORA2B检出仍应保留为本项目接收端测量限制。

## 新增的两段直接作用依据
|候选连接|原始研究与实验|可支持的表述|
|---|---|---|
|MTA—A2B/STAT3依赖髓系作用|Hansen 2022：MTA处理，受体拮抗/敲除、STAT3药理和显性负性干预|MTA引发部分髓系响应需要这些环节；研究未直接证明SPP1响应|
|A2B—SPP1/OPN|Schneider，在线2009、卷期2010：原代小鼠肺泡巨噬，NECA刺激；A2B拮抗和敲除阻断OPN诱导|A2B参与巨噬SPP1诱导具有直接实验依据|
|STAT3—SPP1|Pan 2025：肿瘤相关BMDM体系，Stattic抑制降低Spp1；SPP1启动子野生型/位点突变双荧光素酶试验|支持STAT3参与SPP1转录调控；报告基因试验实际在293T中，不能写为巨噬内ChIP已证实|

来源：
- [Hansen 2022](https://www.nature.com/articles/s41598-022-07697-0)
- [Schneider 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2797041/)；[DOI](https://doi.org/10.1096/fj.09-140772)
- [Pan 2025](https://link.springer.com/article/10.1186/s13046-025-03303-z)

Schneider实验使用NECA，不是MTA。Pan研究背景为VSIG4/乳酸化，不是MTA；不能将其VSIG4或乳酸上游自动接到MTA。Hansen未可靠观测到所测时间点pSTAT3 Tyr705增加，因此候选措辞为“STAT3依赖/转录调控”，不先写“MTA升高pSTAT3”。三个实验体系共享节点，是候选路线的理由，不是整链因果证明。

## 本轮新分析：A2B转录背景未显示稳定效应修饰
以四个冻结实测MTA与转录组匹配队列共255患者（97/52/50/56）做定向检验。研究问题为：A2B受体RNA背景越高，MTA与SPP1的正向联系是否更强？这属于效应修饰检验，不是“MTA使A2B表达升高”的中介检验。

预设主模型：SPP1的秩标准化值 ~ MTA + ADORA2B + MTA×ADORA2B + MTAP + 髓系组成代理。单一预设主要合并交互终点不作FDR；未调整、增加缺氧代理均为敏感性，不以其最优结果替换主分析。

主模型3/4队列交互系数为正；随机效应合并系数0.06718，95%区间[-0.08387,0.21824]，P=0.38335，I²=18.79%。固定合并权重的患者bootstrap区间[-0.05958,0.18808]，亦跨零。
RC18方向较强但未达名义显著（系数0.18924，P=0.05993）；RC20为负向。全部12个队列×模型bootstrap区间均跨零。增加缺氧后的合并系数0.08555，P=0.19777。
当前数据不支持稳定的A2B RNA效应修饰，完整结果和反向队列保留。没有宣称此次公共数据验证了该候选路线。

## STAT3活性资格：本项目存在测量缺口
对冻结的新SCENIC调控子注册表及全细胞AUCell列核查：没有STAT3 regulon，也没有STAT3 AUCell列。因此本轮未计算不存在的STAT3活性、未将STAT3 RNA替代活性、未把其他调控子的靶基因中出现STAT3误认为STAT3调控子。
无推断调控子不等于无生物学活性，不能以此终止STAT3路线。既有JUN/JUNB/CEBPD与SPP1负向结果继续保留，不强行恢复三TF共同驱动模型。

## 备选路线优先级
1. A2B受体依赖路径：MTA端和SPP1端各有髓系直接干预依据，当前机制优先级最高；RNA交互检验阴性。
2. STAT3转录调控路径：同样有分段依据，需补独立定义的活性估计；目前缺本项目活性测量。
3. PRMT5/甲基化：MTA抑制甲基转移酶有生化依据，但本轮没有确认方向正确、对象匹配的PRMT5→髓系SPP1直接证据，不提前指定中介或批量筛选。
4. 乳酸/胆固醇状态：项目已有SPP1端联系，MTA端不足，不将原有代谢相关性冒充MTA响应，也不换名为MTA机制。

## Methods与验证
各队列按既有患者映射匹配ADORA2B TPM，log2(TPM+1)后在模型内作秩标准化；MTA、SPP1及协变量采用队列内秩标准化。交互项为两项秩标准化值乘积，模型包含相应主效应。使用OLS及HC3标准误；各队列重采样患者2000次，每次重算秩，共24000次bootstrap。随机效应合并采用DerSimonian–Laird方差估计，另提供留一队列结果；固定权重bootstrap用于敏感性。
本轮仅一个预设受体候选；不存在按显著性改受体、改评分或改模型。132项QR/逐次bootstrap独立复算通过，源文件记录SHA256。未重聚类或绘图，未使用客户湿实验数据。

## Discussion
bulk受体RNA和SPP1混合多种细胞来源，组织实测MTA不等于巨噬受体处的胞外暴露。交互系数描述观测联系，不能识别受体蛋白活性；四队列合并与置信区间仍属探索性机制对应。受体药理/基因干预与bulk表达效应修饰检验回答不同问题，RNA交互阴性既不验证、也不能排除受体依赖机制。外部干预研究的物种、组织和疾病模型不同。
网页原始证据经本轮在线阅读；服务器下载状态单独保存在sources/retrieval*.json，下载失败不记为原文成功归档。

## 下一轮具体执行顺序
先用独立文献/实验支持的STAT3有符号靶基因资源定义活性指标：事先冻结来源、版本和置信门槛，剔除SPP1、TYMP，检查有效靶基因数和测量覆盖。它是外部定义的验证指标，与本项目新SCENIC分开；不重建网络去追求STAT3出现。
资格合格后，在整体巨噬患者内检验活性—SPP1，深度调整为主，组成/拟时序及研究一致性为敏感性；C4仅作定位。在实测MTA匹配队列检查同一指标的MTA端，bulk与单细胞分别报告。
只有两个端点测量含义一致、来源合格时，才在同一实测MTA队列做MTA—活性—SPP1短路径分解；即使总效应不显著也可探索，但区间、反向及混杂必须完整保留，不做跨队列拼接中介。
"""
(d/"MTA_SPP1_MECHANISM_DEEPENING_CN.md").write_text(report,encoding="utf8")
(d/"CANDIDATE_NARRATIVE_CN.md").write_text("""# 候选机制衔接文字
既有MTA扰动研究提示其髓系作用具有A2B受体和STAT3功能依赖性；独立研究分别提供了A2B参与巨噬SPP1诱导及STAT3调控SPP1转录的实验依据。据此，本项目将A2B受体相关信号及STAT3转录调控列为MTA与SPP1之间优先检验的候选中间环节。
本轮四队列bulk分析未显示稳定的A2B RNA效应修饰，且冻结SCENIC缺少STAT3活性指标，故目前定位为有分段实验依据的机制假说，而非已确认中介。MTA线路继续深化，不因该项RNA检验阴性转移研究目标。
原始来源与完整结果见MTA_SPP1_MECHANISM_DEEPENING_CN.md。
""",encoding="utf8")
# Validate immutable inputs after computation.
inps=json.loads((d/"INPUT_MANIFEST.json").read_text());valid=[]
for row in inps:
 h=hashlib.sha256()
 with Path(row["path"]).open("rb") as f:
  for b in iter(lambda:f.read(8*1024*1024),b""):h.update(b)
 valid.append(h.hexdigest()==row["sha256"])
assert all(valid)
stamp=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
status={"status":"CANDIDATE_MECHANISM_AND_RECEPTOR_TEST_COMPLETE","MTA_SPP1_remains_priority":True,"formal_mediation_run":False,"interaction_patient_n":255,"bootstrap_draws":24000,"numeric_checks_passed":132,"input_hash_checks_passed":len(valid),"STAT3_activity_available":False,"new_figures":0,"source_retry":retrieval,"next":"independently defined STAT3 target activity qualification","updated_at":stamp}
(d/"STATUS.json").write_text(json.dumps(status,indent=2))
bk=d/"root_record_backups";bk.mkdir(exist_ok=True)
entry=f"""
## Phase178 MTA→SPP1中间机制继续深化（{stamp}）
用户明确MTA→SPP1不能放弃，可加入中间环节；免疫/空间分支暂缓。新增原始来源：Schneider2010 A2B干预影响原代肺泡巨噬OPN；Pan2025 STAT3抑制影响BMDM Spp1及293T启动子报告支持。与Hansen2022构成分段候选依据，不是同实验整链。
冻结SCENIC无STAT3调控子/活性列，属于测量缺口，不以STAT3 RNA替代。受体RNA更适合作背景而非MTA诱导表达中介。
四队列255例A2B交互主beta0.06718，95%CI[-0.08387,0.21824]，P0.38335，3/4正；bootstrap亦跨零，当前不支持稳定RNA效应修饰。12模型24000bootstrap，132数值+{len(valid)}输入哈希检查通过。单预设合并主要终点不FDR，无图。
报告：{d}/MTA_SPP1_MECHANISM_DEEPENING_CN.md
下一步独立外部STAT3靶集资格与活性估计，再分别检验髓系SPP1端及实测MTA端；不跨队列拼中介，不恢复旧三TF假说。
"""
for name in ["PROJECT_STATUS.md","DECISIONS.md","CODEX_HANDOFF.md"]:
 p=r/name
 if not (bk/name).exists():shutil.copy2(p,bk/name)
 old=p.read_text()
 if "## Phase178 MTA→SPP1中间机制继续深化" not in old:p.write_text(entry+"\n"+old,encoding="utf8")
p=r/"AGENTS.md"
if not (bk/p.name).exists():shutil.copy2(p,bk/p.name)
tag="## 2026-09-16 用户明确继续深化MTA到SPP1"
if tag not in p.read_text():
 with p.open("a",encoding="utf8") as f:f.write("\n\n"+tag+"\nMTA→SPP1继续作为优先问题，可按直接作用依据引入中间机制。免疫与空间下游暂缓，不能因已有弱关联就终止接收端深化。机制节点不强制等同RNA表达中介；新候选先定来源/方向与测量资格，保留阴性，不跨队列拼成同患者中介。\n")
p=r/"RESULTS_MANIFEST.tsv"
if not (bk/p.name).exists():shutil.copy2(p,bk/p.name)
if "phase178_MTA_SPP1_mechanism" not in p.read_text():
 with p.open("a",newline="") as f:csv.writer(f,delimiter="\t").writerow(["phase178_MTA_SPP1_mechanism","L2_candidate_mechanism","report",str(d/"MTA_SPP1_MECHANISM_DEEPENING_CN.md"),stamp,"PASS","132 numeric checks;9 input hashes","A2B RNA interaction not supported; STAT3 activity measurement gap; MTA retained"])
de=d/"delivery";de.mkdir(exist_ok=True)
for name in ["MTA_SPP1_MECHANISM_DEEPENING_CN.md","CANDIDATE_NARRATIVE_CN.md","RECEPTOR_INTERACTION_META.tsv","RECEPTOR_INTERACTION_COHORT.tsv","RECEPTOR_INTERACTION_LEAVE_COHORT.tsv","STAT3_ACTIVITY_GATE.json","STATUS.json"]:shutil.copy2(d/name,de/name)
hh={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in de.iterdir()};(d/"DELIVERY_SHA256.json").write_text(json.dumps(hh,indent=2));print(json.dumps({"status":status,"hashes":hh},ensure_ascii=False))
