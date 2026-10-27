from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,shutil,csv,datetime
r=Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1");d=r/"phase184_mature_macrophage_perturbation_20260916_v1"
va=json.loads((d/"VALIDATION_STATUS.json").read_text());assert va["all_pass"]
R=pd.read_csv(d/"RNA_TARGET_SUMMARY.tsv",sep="\t");A=pd.read_csv(d/"ARRAY_EFFECT_SUMMARY.tsv",sep="\t")
def tab(f):
 def fmt(x):return f"{x:.6g}" if isinstance(x,(float,np.floating)) else str(x)
 return "| "+" | ".join(f.columns)+" |\n|"+"|".join(["---"]*len(f.columns))+"|\n"+"\n".join("| "+" | ".join(fmt(x) for x in row)+" |" for row in f.itertuples(index=False,name=None))
report="""# 成熟人巨噬细胞扰动调查：IL-10添加、中和与STAT3敲低
## 结论
本轮完成三类公开干预读出的定向核查。72小时IL-10添加后SPP1在3名严格配对供体中均下降；在另一M-CSF来源巨噬体系，IL-10中和没有显示SPP1下降，STAT3敲低后的SPP1方向随有无LPS而不同。
这些结果不支持将“IL-10—STAT3—SPP1”作为跨条件普遍成立的中介模型。既有本项目STAT3代理—SPP1相关及Phase183正向实例仍保留，但不能把一般IL-10反应或STAT3依赖的其他基因输出直接当作SPP1响应。所有数据均未含MTA处理，不填补MTA箭头。

## 一、72小时IL-10添加：3名供体SPP1均下降
[GSE46903官方记录](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46903)
采用GM-CSF预分化的人巨噬细胞，只纳入IL10 72h与同供体con72h。BC58、BC76、BC78身份在GEO明确标注，六个样本逐项匹配。
SPP1两条探针ILMN_1651354与ILMN_2374449均在三名供体同向下降，处理与对照的两条探针均通过原作者技术检测背景。按每样本所有合格探针中位数汇总：
"""+tab(A)+"""
SPP1中位log2归一化强度差为-0.4030，供体范围[-1.06055,-0.14590]。每名供体结果独立保留；n=3只作方向/效应描述，无新增P/FDR。
SOCS3和STAT3 RNA在72h也均下降，不将其解释为早期受体信号不存在。该数据没有本轮可用的同供体IL-10早期时间序列，不能据此确定何时发生方向变化。LGALS3为2正1负，不支持一致升高。
GSE46903为广谱巨噬刺激数据，当前只分析预先选定的同培养背景同时间对照，不与M-CSF预分化或其他时间合并。

## 二、IL-10中和：SPP1没有出现预期下降
[GSE181250官方记录](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181250)
该研究将原代单核培养7天，形成M-CSF来源巨噬；以anti-IL10或同型IgG预处理1h，再加LPS 4h。每组3份独立制备，公开元数据没有明确逐供体对应，未凭MSF21–26编号强行配对。
SPP1组中位log2(FPKM+1)差为+0.28864，处理和对照范围明显重叠，不能称为稳定上调，也未显示预期的下降。SOCS3与STAT3 RNA降低且处理组各值均低于对照范围，说明其他转录读出有变化，但不能替代SPP1终点。
"""+tab(R[R.study=="GSE181250"][["gene","n_treated","n_control","median_delta","treated_min","treated_max","control_min","control_max"]])+"""

## 三、STAT3敲低：操作读出明确，SPP1方向仍依赖条件
[GSE180897官方记录](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180897)
同一研究在M-CSF培养第5天转染siSTAT3或siCNT，48h后加LPS 4h或保持未刺激；各组3份独立制备。虽然有#1–3标号，逐供体身份未明确，主报告仍采用组水平描述。
STAT3 RNA在未刺激和LPS背景下分别减少约2.192和2.280个log2(FPKM+1)单位，处理组全部低于对照范围，支持RNA敲低操作有效。SPP1则未刺激组中位差-0.38527，LPS组+0.43665；两者处理/对照范围均重叠，不作显著性或“所有供体变化”判断。
"""+tab(R[R.study=="GSE180897"][["contrast","gene","median_delta","treated_min","treated_max","control_min","control_max"]])+"""
本轮核查的是RNA，不直接证明STAT3磷酸化或蛋白功能完全消失；不将siRNA单实验泛化为所有STAT3作用的否定。

## 四、与原始论文的关系
两组RNA数据来自同一[Cuevas等研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC9149444/)，不称为两个独立外部研究。论文结论涉及一组IL-10/STAT3依赖的基因程序，并不意味着该程序中每个候选、尤其本项目指定的SPP1都必然同向。
本轮从公开FPKM中提取SPP1/LGALS3/SOCS3/STAT3，是额外的固定终点描述，不冒充原文已宣称的SPP1机制。

## Methods
主要终点预先固定SPP1，辅助LGALS3/SOCS3/STAT3。先核查对象、供体、背景和干预，再读取目标表达。RNA资料是RSEM输出FPKM而非原始计数；用log2(FPKM+1)描述组中位差和全部样本范围，不再重算CPM，也不作原始计数模型。
GSE46903使用作者发布的log2转换、分位数归一化强度，六样本均存在完整table_end。GPL6947官方完整注释确认Entrez与基因符号一致，保留所有单一基因注释且六样本共同测得的探针。SPP1两条、LGALS3两条、SOCS3一条、STAT3三条；每样本多探针取中位，逐探针结果另存。另一个SOCS3探针ILMN_2156250虽在平台注释，但未出现在作者发布表达表，明确记录缺失，不填零、不按效应选择探针。
三个小样本数据不进行P/FDR检验，不训练全基因签名，不投射本项目评分。技术detection P为作者微阵列测量指标，不当作生物学显著性。

## 下载与验证
完整series matrix和未压缩平台注释下载超时，部分文件保存但未用于结果。改取六个完整GSM小文件；压缩平台注释通过HTTP范围下载补齐，逐段核对Content-Range、长度并完整gzip解压校验后使用。未经确认完整的文件不用于探针映射。
共173项数值/独立解析检查通过，完整输入SHA256保存在COMPLETE_INPUT_MANIFEST.json。原版因未发布SOCS3探针而在分析阶段中止；修正版先审查共同测量覆盖，保留缺失记录和原脚本。没有新图、无客户实验使用。

## MTA并行检索的结果与边界
本轮继续检索直接MTA或腺苷受体干预的髓系表达资源，未新增合格直接MTA—SPP1转录组。检索命中的GSE272546属于腺苷四磷酸/ATP处理，不能因名称含adenosine就当作腺苷受体或MTA实验。原MTA髓系功能实验及其受体/STAT3依赖文献保留；IL-10体系不能替代它们。
这是一轮有界调查，未证明所有公共资源都不存在相应实验；未联系研究作者、未接管客户湿实验。

## 研究取舍与下一断点
- 本项目TYMP—STAT3肿瘤内/组织层关联、整体巨噬STAT3—SPP1关联、C4表型和所有关键阴性保留。
- IL-10背景可测，但经两轮直接数据核查后，不能将它指定为MTA到SPP1的通用必经节点。LGALS3也不作强制中介。
- Phase183的IL-10正向系列继续保留，本轮成熟巨噬数据增加了不同方向的实例；不能选择性只保留任何一侧。
- 下一步优先回到MTA本身：直接MTA暴露、胞外条件、受体功能依赖及SPP1具体输出。停止用追加一般IL-10刺激数据或更换泛STAT3评分来追求完整线性链。
- 从上述数据尚不能确定为何各体系SPP1方向不同；培养因子、成熟程度、时间、供体和炎症背景是待区分因素。提出这些因素不等于已验证其解释。
"""
(d/"MATURE_MACROPHAGE_PERTURBATION_RESULTS_CN.md").write_text(report,encoding="utf8")
(d/"EVIDENCE_UPDATE_CN.md").write_text("""# 证据更新
独立成熟人巨噬干预没有支持普遍的IL-10—STAT3—SPP1正向链：GM-CSF体系IL10 72h后3名配对供体SPP1均下降，两探针一致；M-CSF体系IL10中和未见SPP1下降；STAT3 RNA敲低有效，但SPP1在静息和LPS背景组间方向不同、范围重叠。
均为小样本探索性描述，没有新增P/FDR，未涉及MTA处理。保留Phase183正向实例及本项目关联，不把一般IL10响应当MTA验证，也不把本轮阴性泛化为STAT3在所有条件下无作用。后续重点回到直接MTA作用及SPP1输出，免疫/空间分支仍暂缓。
""",encoding="utf8")
allstats="# 成熟巨噬干预完整描述统计\n\n## IL10添加配对结果\n\n"+tab(pd.read_csv(d/"ARRAY_DONOR_EFFECTS.tsv",sep="\t"))+"\n\n## 功能阻断组间结果\n\n"+tab(R)+"\n\n## 功能阻断所有样本\n\n"+tab(pd.read_csv(d/"RNA_TARGET_SAMPLE_VALUES.tsv",sep="\t"))
(d/"ALL_TARGET_STATISTICS_CN.md").write_text(allstats,encoding="utf8")
stamp=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
status={"status":"COMPLETE_VALIDATED","studies":["GSE46903","GSE181250","GSE180897"],"array_paired_donors":3,"RNA_preparations_per_group":3,"RNA_donor_pairing_confirmed":False,"new_P":0,"new_FDR":0,"checks":173,"direct_MTA_data_added":False,"new_figures":0,"updated_at":stamp}
(d/"STATUS.json").write_text(json.dumps(status,indent=2))
bk=d/"root_record_backups";bk.mkdir(exist_ok=True)
entry=f"""
## Phase184 成熟人巨噬三类扰动核查完成（{stamp}）
GSE46903 GM-CSF背景IL10 72h，同BC58/76/78配对，SPP1 3/3降低、中位log2差-0.403，两条探针均同向且可检出。六GSM和完整GPL注释校验，损坏大矩阵未使用。
GSE181250 MCSF7天巨噬IL10中和+LPS4h：SPP1组中位差+0.28864、范围重叠，未显示预期下降。GSE180897 STAT3敲低RNA明显降低(约-2.2log2)，SPP1静息-0.38527、LPS+0.43665，范围重叠。RNA各组3独立制备，逐供体映射未知不强配；同论文不称两次独立研究。
共173核查通过，无P/FDR、无图、无签名训练。MTA直接转录组未新增，GSE272546 ATP/四磷酸不是MTA/腺苷受体替代。
报告：{d}/MATURE_MACROPHAGE_PERTURBATION_RESULTS_CN.md
前期MTA与STAT3正负结果全保留。一般IL10路线不作通用必经中介；下一断点回到直接MTA暴露/胞外条件/受体功能与SPP1具体输出，停止泛IL10评分追阳性，免疫空间暂缓。
"""
for n in ["PROJECT_STATUS.md","DECISIONS.md","CODEX_HANDOFF.md"]:
 p=r/n
 if not (bk/n).exists():shutil.copy2(p,bk/n)
 old=p.read_text()
 if "## Phase184 成熟人巨噬三类扰动核查完成" not in old:p.write_text(entry+"\n"+old,encoding="utf8")
p=r/"RESULTS_MANIFEST.tsv"
if not (bk/p.name).exists():shutil.copy2(p,bk/p.name)
if "phase184_mature_macrophage" not in p.read_text():
 with p.open("a",newline="") as f:csv.writer(f,delimiter="\t").writerow(["phase184_mature_macrophage","L2_direct_perturbation","report",str(d/"MATURE_MACROPHAGE_PERTURBATION_RESULTS_CN.md"),stamp,"PASS","173 checks;complete annotation;6 paired arrays","No universal IL10-STAT3-SPP1 support;no MTA intervention"])
de=d/"delivery";de.mkdir(exist_ok=True)
for n in ["MATURE_MACROPHAGE_PERTURBATION_RESULTS_CN.md","EVIDENCE_UPDATE_CN.md","ALL_TARGET_STATISTICS_CN.md","RNA_TARGET_SUMMARY.tsv","ARRAY_EFFECT_SUMMARY.tsv","ARRAY_DONOR_EFFECTS.tsv","VALIDATION_STATUS.json","STATUS.json"]:shutil.copy2(d/n,de/n)
hh={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in de.iterdir()};(d/"DELIVERY_SHA256.json").write_text(json.dumps(hh,indent=2));print(json.dumps({"status":status,"hashes":hh},ensure_ascii=False))
