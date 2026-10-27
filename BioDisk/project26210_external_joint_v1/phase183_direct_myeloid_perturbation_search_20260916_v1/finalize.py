from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,shutil,datetime,csv
r=Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1");d=r/"phase183_direct_myeloid_perturbation_search_20260916_v1"
val=json.loads((d/"VALIDATION_STATUS.json").read_text());assert val["all_pass"]
S=pd.read_csv(d/"TARGET_EFFECT_SUMMARY.tsv",sep="\t")
def md(f):
 def fmt(x):return f"{x:.6g}" if isinstance(x,(float,np.floating)) else str(x)
 return "| "+" | ".join(f.columns)+" |\n|"+"|".join(["---"]*len(f.columns))+"|\n"+"\n".join("| "+" | ".join(fmt(x) for x in z)+" |" for z in f.itertuples(index=False,name=None))
spp=S[(S.gene=="SPP1")&(S.method=="median_ratio")&S.contrast.isin(["IL10_vs_Rest","IL10_plus_stim_vs_stim"])]
report="""# 直接髓系扰动调查：人源IL-10处理提供SPP1方向不一致的独立实例
## 核心进展
本轮取得并复核GSE261243公开人源髓系IL-10扰动数据，获得独立于本项目观察队列的直接处理对照。SPP1在一个实验系列3/3供体上升，另一个系列3/3下降；两种预设归一化方法方向一致。当前支持“响应依赖具体实验条件”的探索性判断，不支持IL-10普遍提高SPP1。
这不是MTA扰动验证。当前有界检索仍未取得合格的“直接MTA处理髓系、可配对独立重复、含SPP1全转录组读出”资源；不宣称全球不存在该数据。既有MTA、TYMP—STAT3、STAT3—SPP1等结果完整保留。

## 一、资源资格与排除理由
|资源|对象/设计|本轮处理|
|---|---|---|
|GSE261243|人原代CD14单核细胞，M-CSF过夜培养，IL-10处理；LPS/TNF两个系列各3供体×4条件|通过定向描述资格；下载矩阵并分析|
|GSE130070|人单核LPS±IL10R阻断；仅一名供体有完整对照，另一名只有阻断组|不把细胞数或第二供体单组当有效配对重复|
|GSE51446|人单核来源巨噬多刺激微阵列，7条件各1样本|缺少有效重复，不做推断验证或签名训练|
|GSE86167|小鼠BMDM的IL-10/LPS RNA-seq，两套重复，属GSE86170|候选补充，未混入人源结果|
|GSE86168/86169|同研究H3K27ac/STAT3 ChIP|不是转录表达矩阵，不代替SPP1表达响应|
|GSE34553|小鼠DSS结肠组织±口服MTA|直接MTA暴露但为混合组织，不作直接巨噬响应签名|
|GSE249082（既有）|MTA处理人成纤维细胞|继续排除直接巨噬签名用途|

上述GEO设计依据来自官方元数据，六份系列SOFT及GSE261243完整family元数据保存在sources。某论文中同时存在RNA-seq和MTA药理实验，并不意味着公开RNA-seq包含MTA组；未证实处理映射的命中不纳入直接MTA数据。

## 二、GSE261243实际对象、供体和处理
[GEO官方记录](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261243)
对应[Mishra等，Nature Immunology 2025](https://www.nature.com/articles/s41590-025-02137-3)，研究IL-10对人单核炎症响应的影响。本轮结论为我们对公开矩阵的预设目标基因复核，不冒充原文已经声称MTA或SPP1机制。

- 对象：人CD14+原代单核细胞，M-CSF 20 ng/ml培养过夜；不是成熟肿瘤巨噬，不等同本项目C4。
- IL-10：100 ng/ml过夜；随后LPS 10 ng/ml处理3小时，或TNF 20 ng/ml处理6小时。
- LPS系列：D562、D563、D565，各有Rest、IL10、LPS、IL10+LPS，共12文库。
- TNF系列：D1、D2、D3，各有Rest、IL10、TNF、IL10+TNF，共12文库。
- 两系列供体身份对应未知，不合计为已确认6名独立供体；不把24文库当24个独立样本。
- 同系列以供体后缀配对。LPS对照元数据名Rest_None_Dxxx与矩阵Rest_Dxxx的命名差异逐项显式映射。
- GEO另有6个IRF基因编辑样本，本轮不纳入，不与IL-10处理混组。
- IL-10单独/Rest具体收样时点未单独在逐样本说明中完全展开，按所属实验系列对照分析，不能假定两系列所有条件除刺激外完全相同。

## 三、SPP1方向在两个实验系列相反
主要效应定义为供体内log2(归一化计数+1)差，不是蛋白分泌量或临床效应。每个系列n=3，仅报告探索性方向和幅度：
"""+md(spp[["background","contrast","n_donors","positive","negative","median_delta","min_delta","max_delta"]])+"""
LPS系列IL10单独处理相对Rest的中位差+0.8385，联合处理相对LPS为+0.9700，均3/3上升。TNF系列对应中位差-3.0197和-1.7391，均3/3下降。
关键点：连“IL10单独vsRest”也已在两系列之间反向，不能把全部差异归因于LPS与TNF。供体、实验系列、取样条件或其他未记录因素均待区分，本轮没有确定差异的原因。
SPP1在全部24个纳入样本有较充足计数，最低原始计数146；不存在因基因未提供而填0的问题。

## 四、归一化敏感性不改变方向
"""+md(S[(S.gene=="SPP1")&S.contrast.isin(["IL10_vs_Rest","IL10_plus_stim_vs_stim"])][["background","method","contrast","n_donors","positive","negative","median_delta"]])+"""
上四分位数方法保留所有主要对比的方向。它与中位比值法都基于同一提供矩阵，不能替代绝对RNA量或外源参照归一化；本轮未因结果更换方法。

## 五、辅助读出显示SPP1不是泛IL-10/STAT3表达响应的同义词
"""+md(S[(S.method=="median_ratio")&S.gene.isin(["LGALS3","SOCS3","STAT3"])&S.contrast.isin(["IL10_vs_Rest","IL10_plus_stim_vs_stim"])][["background","gene","contrast","n_donors","positive","negative","median_delta"]])+"""
SOCS3和STAT3 RNA在两个系列的IL-10相关对比均3/3上升，SPP1却方向相反。因此，出现IL-10相关转录反应不保证SPP1上调。此处STAT3是RNA，不是pSTAT3或Phase179的靶基因活性代理。
LGALS3在两个系列IL10单独处理均3/3下降，LPS背景也3/3下降，TNF背景2/3下降。这是本体系的转录读出，不支持把LGALS3同步上升预设为普遍必经步骤；不能据此否定心肌修复研究中的蛋白功能依赖。

## Methods
先固定资源资格、供体配对、主要SPP1与辅助LGALS3/SOCS3/STAT3，后读取目标结果。LPS矩阵11059基因、TNF矩阵10811基因，均为过滤后的原始整数计数，gene_id与gene_name无重复。本轮不对提供矩阵总和声称全库测序深度。
各系列分别采用中位比值法：所有样本非零的基因计算跨样本几何均值，再对逐样本基因比值取中位，size factor几何均值缩放为1；参与基因数LPS10496、TNF10597。非零计数上四分位数因子为预设敏感性。
主要对比IL10 vs Rest、IL10+刺激 vs刺激；刺激vsRest、联合vsRest与差中差另存诊断表。按供体报告log2归一化计数加1差；不跨系列合并、无新P/FDR、无全基因差异发现、无响应签名训练/转移。

## 独立核查
完整GSM—矩阵列—处理—供体映射24/24；四个目标基因均实际提供。原始计数由gzip文本独立解析，size factor采用等价形式核对，效应由归一化计数比独立复算，摘要与文件SHA256核查。669项技术检查通过，最大数值差1.78e-15。没有新图，也未改动前期患者队列。

## 研究取舍
1. 保留LPS系列的正向干预实例，作为“某些人源髓系条件下IL-10可伴随SPP1增加”的探索性依据；必须同时保留TNF系列反向结果。
2. 不把IL-10响应、M2标签、STAT3 RNA或泛STAT3代理直接命名为SPP1诱导程序；不能用一个正向系列训练签名再在原队列寻找支持。
3. MTA仍需直接端点证据。本轮没有MTA处理，不将上述IL-10结果填入MTA箭头；A2B受体依赖文献保留为另一类证据。
4. 下一断点优先核查人原代单核来源成熟巨噬的独立IL-10时间/条件对照，或真正包含MTA/受体干预及SPP1读出的数据。先辨别是否可重复，随后才考虑背景交互。既有GSE46903等广谱刺激资源需先确认IL10条件、供体与对照，不能只凭数据集名称直接纳入。
5. 当前不追加本项目更多相关性模型追求正向。已有TYMP—STAT3、整体巨噬STAT3—SPP1、C4表型与MTA端未获支持的结果均继续保留，免疫/空间分支仍暂缓。
"""
(d/"DIRECT_PERTURBATION_FINDINGS_CN.md").write_text(report,encoding="utf8")
(d/"EVIDENCE_UPDATE_CN.md").write_text("""# 新证据定位
本轮获得人原代CD14髓系IL-10直接扰动数据的供体配对读出，独立于项目观察队列：LPS实验系列SPP1在IL10单独和联合条件下均3/3上升，TNF实验系列均3/3下降。两种归一化方法一致；每系列仅3供体，按探索性方向和效应报告，无P/FDR。
因此保留IL-10可在部分髓系条件下伴随SPP1上升的直接实例，但不能声称普遍诱导，不能将系列差异简单归因于LPS/TNF，也不能把它写成MTA干预证据。SOCS3/STAT3 RNA上升与SPP1方向并不等价；LGALS3转录同步上升不是本体系普遍现象。MTA以及前期所有正负结果继续保留。
""",encoding="utf8")
(d/"ALL_TARGET_STATISTICS_CN.md").write_text("# 目标基因完整描述统计\n\n无新增P值或FDR；差值为log2(归一化计数+1)。\n\n"+md(S),encoding="utf8")
stamp=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
status={"status":"COMPLETE_VALIDATED_BOUNDED_INVESTIGATION","GSE261243_libraries":24,"donors_per_series":3,"cross_series_donor_independence":"unresolved;not pooled","P_values_generated":0,"FDR_values_generated":0,"technical_checks":669,"direct_MTA_macrophage_transcriptome_qualified":False,"IL10_SPP1_uniform_upregulation_supported":False,"new_figures":0,"updated_at":stamp}
(d/"STATUS.json").write_text(json.dumps(status,indent=2))
bk=d/"root_record_backups";bk.mkdir(exist_ok=True)
entry=f"""
## Phase183 直接髓系扰动调查与人IL10数据核查完成（{stamp}）
本轮筛查官方GEO6系列，取得GSE261243人CD14单核MCSF过夜培养的直接IL10 RNAseq；LPS/TNF各3供体4条件，24库不作24独立n，两系列供体对应未知不合并。
SPP1：LPS系列IL10vsRest中位log2差+0.8385、IL10+LPSvsLPS+0.9700，均3/3升；TNF系列对应-3.0197/-1.7391，均3/3降。上四分位数敏感性同向。IL10单独也反向，不将差异单独归因LPS/TNF。SOCS3/STAT3 RNA两系列均升，LGALS3单独IL10均降。
各n3仅描述，无P/FDR、无签名训练/转移。669独立/计算检查通过。未取得合格直接MTA髓系SPP1转录组，GSE34553结肠/GSE249082成纤维不冒充巨噬，GSE130070配对重复不足/GSE51446每条件1样本。
报告：{d}/DIRECT_PERTURBATION_FINDINGS_CN.md
下一断点：独立成熟人巨噬IL10/时间背景或直接MTA受体扰动数据资格，不选阳性系列训练签名，不以IL10结果填MTA箭头。此前正负结果全保留，无图；免疫空间暂缓。
"""
for n in ["PROJECT_STATUS.md","DECISIONS.md","CODEX_HANDOFF.md"]:
 p=r/n
 if not (bk/n).exists():shutil.copy2(p,bk/n)
 old=p.read_text()
 if "## Phase183 直接髓系扰动调查" not in old:p.write_text(entry+"\n"+old,encoding="utf8")
p=r/"RESULTS_MANIFEST.tsv"
if not (bk/p.name).exists():shutil.copy2(p,bk/p.name)
if "phase183_direct_perturbation" not in p.read_text():
 with p.open("a",newline="") as f:csv.writer(f,delimiter="\t").writerow(["phase183_direct_perturbation","L2_external","report",str(d/"DIRECT_PERTURBATION_FINDINGS_CN.md"),stamp,"PASS","669 checks;24 GSM mapped;no inferential tests","IL10 perturbation series discordant;not MTA evidence"])
de=d/"delivery";de.mkdir(exist_ok=True)
for n in ["DIRECT_PERTURBATION_FINDINGS_CN.md","EVIDENCE_UPDATE_CN.md","ALL_TARGET_STATISTICS_CN.md","TARGET_EFFECT_SUMMARY.tsv","TARGET_DONOR_EFFECTS.tsv","GSE261243_SAMPLE_AUDIT.tsv","VALIDATION_STATUS.json","STATUS.json"]:shutil.copy2(d/n,de/n)
hh={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in de.iterdir()};(d/"DELIVERY_SHA256.json").write_text(json.dumps(hh,indent=2));print(json.dumps({"status":status,"hashes":hh},ensure_ascii=False))
