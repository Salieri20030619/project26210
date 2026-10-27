from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,shutil,csv,datetime
r=Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1");d=r/"phase182_IL10_LGALS3_receiver_audit_20260916_v1"
valid=json.loads((d/"VALIDATION_STATUS.json").read_text());assert valid["all_pass"]
R=pd.read_csv(d/"RECEIVER_COVERAGE_SUMMARY.tsv",sep="\t");S=pd.read_csv(d/"EDGE_SUMMARY.tsv",sep="\t");source=pd.read_csv(d/"SOURCE_PANEL_SUMMARY.tsv",sep="\t");P=json.loads((d/"C4_LGALS3_SUMMARY.json").read_text())
def md(f):
 def fmt(x):return f"{x:.6g}" if isinstance(x,(float,np.floating)) else str(x)
 return "| "+" | ".join(f.columns)+" |\n|"+"|".join(["---"]*len(f.columns))+"|\n"+"\n".join("| "+" | ".join(fmt(x) for x in row)+" |" for row in f.itertuples(index=False,name=None))
src=source[(source.gene=="IL10")&source.source_class.isin(["C4","Monocyte_refined","OtherMacrophage_refined","T cell","Malignant_conservative"])]
primary=S[(S.resource=="collectri")&(S.scoring=="zmean")]
report="""# IL-10背景、巨噬受体与LGALS3定向核查
## 结论
**IL-10接收背景通过预设RNA可测性门槛；LGALS3在C4中富集，但当前数据不支持把它设为STAT3→SPP1的必经正向中间环节。**
本轮新增来源/覆盖与三条关联检查，不新增MTA交互或中介模型。MTA候选、TYMP—STAT3和原STAT3—SPP1正向关联，以及MTA端阴性和跨区室弱关联全部保留。

## 一、IL-10候选RNA来源以可测髓系为重点，不能默认来自肿瘤细胞
在原75文库258120细胞图谱中，按肿瘤内患者和细胞类聚合原始计数。新C4、其他成熟巨噬和单核采用冻结精细成员，其他类沿用原图谱标签。
"""+md(src[["source_class","n_patients","patients_detected","patients_ge10_detected","median_detection","median_observed_RNA_share"]])+"""
表中median_detection为0–1比例；各类纳入患者不同，RNA贡献的中位数不能相加成组成图。
其他成熟巨噬29/34、精细单核22/29、T细胞15/33名患者达到至少10个IL10检出细胞；保守恶性细胞仅3/20达到该门槛，9/20至少有一次检出。故目前不宜将IL10直接指定为TYMP高肿瘤细胞发出的配体。髓系及部分T细胞提供了需要进一步验证的微环境来源背景。

所捕获IL10 RNA不是胞外IL-10浓度，也未证明分泌或自分泌。IL10表达与其受体在同类细胞存在，不能直接构成自分泌回路。旧未纳入精细髓系的细胞单列保留，罕见细胞类和少量检出不作生物学来源定论；完整表保留。

## 二、IL-10受体可用于后续RNA研究，A2B不适合当前单细胞高低分组
34名患者各有≥20个肿瘤内成熟巨噬；按至少10个检出细胞、且至少10名患者可用，定义技术可测性门槛：
"""+md(R[R.feature.isin(["IL10","IL10RA","IL10RB","IL10RA_AND_IL10RB","IL6","IL6R","IL6ST","IL6R_AND_IL6ST","SOCS3","ADORA2B"])])+"""
IL10RA、IL10RB分别33、29名患者达到单链门槛，两链共同检出在23名患者达到门槛，共检出细胞比例中位数4.787%。支持后续针对IL-10响应背景的研究，不证明功能性受体复合物活性。
ADORA2B只有3名患者达到门槛，中位检出比例0.273%。当前不支持据此硬划单细胞A2B高/低组，也不能用RNA低检出否定A2B蛋白或功能。
IL6R/IL6ST两链共同检出在19名患者达到门槛，可作为对照背景；此共同RNA检出指标不覆盖IL-6反式信号。SOCS3是反馈相关候选，不单凭表达升高定义通路更强。

## 三、防循环审查通过，无须更换旧STAT3评分
LGALS3原本不在CollecTRI205靶或DoRothEA62靶中，故本次“剔除LGALS3”审查未改变靶集。TYMP/SPP1/STAT3自身继续排除。
原STAT3—SPP1逐患者效应已与Phase179最终结果复核一致。本轮为预设三条边新增BH家族，旧单主要终点统计定位不被改写。这里的STAT3仍是有符号靶基因转录代理，不是磷酸化实测。

## 四、LGALS3—SPP1存在小幅正向联系，受C4组成影响
主深度调整分析31患者20正、11负，rho中位数0.04996，P/FDR=0.02889。三边家族FDR通过，但中位效应bootstrap区间[-0.00346,0.06238]跨零，患者一致性有限。Wilcoxon秩检验与中位数bootstrap针对的统计量不同，不能只择其有利者。
增加C4指示变量后，rho中位数降至0.01850，P/FDR=0.21336；与群体状态组成影响相容，不把调整前关联全部解释为群内调控。
同分区拟时序可分析子集12名患者中10名正向，rho=0.02852，P=0.02537、FDR=0.03805；子集的正向信号保留，不代替31患者的组成敏感性。
"""+md(primary[primary.edge=="LGALS3_SPP1"][["adjustment","n","median_rho","positive","negative","P","FDR_three_edges","CI_low","CI_high"]])+"""

## 五、STAT3—LGALS3主分析方向相反，正向串联不成立
主深度调整33名患者27负、6正，rho中位数-0.10947，P=2.6817e-05、FDR=4.0225e-05；增加检出基因数或C4组成后主评分仍以负向为主。
"""+md(primary[primary.edge=="STAT3_LGALS3"][["adjustment","n","median_rho","positive","negative","P","FDR_three_edges"]])+"""
DoRothEA标准化均值rho=-0.03545、P/FDR=0.06313；秩均值rho=-0.02701、P/FDR=0.20458，均未支持正向。主评分同分区拟时序子集14患者7正7负、P/FDR=0.37947。
因此，当前不支持将心肌修复文献中的STAT3→LGALS3→SPP1正向顺序直接搬到本项目。负相关也不能反过来写成STAT3直接抑制LGALS3的机制。

## 六、LGALS3是值得保留的C4状态特征
24名可配对患者中LGALS3均在C4较其他成熟巨噬升高，中位log2(CPM+1)差+1.98519，P=1.9422e-05。该单一预设配对终点不另作FDR。
其C4定位价值与“是否是STAT3必经中介”是两个问题。C4标签和成员保持冻结，不根据本轮LGALS3再定义细胞。

## 七、既有STAT3—SPP1结果保留
主深度调整31患者29正，rho=0.08087，P=4.3194e-06，本轮三边家族FDR=1.2958e-05。控制C4组成后仍29/31正，说明该已有联系未因本轮候选核查被替换。
MTA—STAT3代理的既有阴性、TYMP肿瘤内/组织层正向、同文库跨区室弱联系均不变。

## Methods
固定面板IL10/IL10RA/IL10RB/IL6/IL6R/IL6ST/SOCS3/ADORA2B/LGALS3/SPP1。全图谱来源按患者和细胞类聚合，≥20细胞作描述摘要；保留低资格类原始行，无来源差异P。受体共同检出要求同一细胞两基因均有原始计数，按患者总结，不用细胞数作独立重复。
三条关联在肿瘤内成熟巨噬、患者≥50细胞、单基因各≥10检出条件下计算；STAT3代理采用既有有符号靶集，主患者内z均值，资源及秩均值敏感性。测序深度主调整，另控基因数、C4指示或同分区拟时序；分区效应取患者中位。患者层双侧Wilcoxon近似加连续性校正；5000次患者bootstrap给出中位效应区间。BH在同资源/评分/调整的三条关联内计算。各研究留出仅描述效应。
C4/其他成熟巨噬每侧≥20细胞，LGALS3采用原始计数聚合后的log2(CPM+1)配对差。研究单位始终为患者。

## 审计与交付
固定全图谱SHA256与此前来源记录一致；直接回查80个随机细胞的10基因原始计数，检查全部患者基因总量与RNA贡献和为1。关联以QR独立复算、BH及摘要回算，原STAT3—SPP1逐患者结果一致。2974项技术检查全部通过，最大数值差2.22e-16。未生成图、未重聚类、未使用客户实验数据。

## Discussion与下一断点
IL-10候选来源和受体表达为背景研究提供了可测条件，但未测得胞外IL-10、受体活性或MTA调节作用。全图谱RNA来源仍受细胞捕获、检测深度和环境RNA影响；不同类的RNA份额不能当分泌贡献。受体共检出也不证明空间接触或信号传递。
LGALS3富集与主STAT3—LGALS3负向说明简单线性串联不符合当前观测，不能靠改变评分/分群继续追求预期正向。仍保留LGALS3作为C4功能表型线索。
后续优先核对合格的IL-10/腺苷受体髓系扰动资源和实际响应输出，明确“MTA是否影响背景依赖响应”，再决定是否有依据开展MTA×背景检验。不以当前巨噬IL10 RNA直接代替暴露量，不把技术门槛通过称为该机制成立。MTA及此前两段继续保留，免疫与空间分支仍暂缓。
"""
(d/"IL10_LGALS3_AUDIT_RESULTS_CN.md").write_text(report,encoding="utf8")
(d/"MANUSCRIPT_MODULE_CN.md").write_text("""# IL-10接收背景可测，LGALS3定位于C4相关状态
来源与受体核查显示，精细注释髓系细胞及部分T细胞可检出IL10 RNA，巨噬IL10RA/IL10RB具有后续RNA分析覆盖；当前不将IL10默认归为肿瘤细胞分泌信号。LGALS3在24名配对患者的C4中一致富集，中位log2(CPM+1)差+1.985（P=1.94×10⁻⁵）。

整体巨噬LGALS3—SPP1主关联幅度较小，控制C4组成后减弱；STAT3转录活性代理—LGALS3主关联为负向。因此，LGALS3保留为C4状态特征，当前不支持其作为STAT3正向连接SPP1的必经中介。既有STAT3—SPP1关联与MTA候选线路保留，背景依赖的MTA响应仍需另行验证。
""",encoding="utf8")
stats="# 完整统计与来源摘要\n\n"
for name in ["RECEIVER_COVERAGE_SUMMARY.tsv","SOURCE_PANEL_SUMMARY.tsv","EDGE_SUMMARY.tsv","C4_LGALS3_PATIENTS.tsv","LEAVE_STUDY.tsv"]:stats+="## "+name+"\n\n"+md(pd.read_csv(d/name,sep="\t"))+"\n\n"
(d/"ALL_STATISTICS_CN.md").write_text(stats,encoding="utf8")
stamp=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
status={"status":"COMPLETE_VALIDATED","checks":2974,"IL10_receptor_RNA_gate":True,"ADORA2B_RNA_gate":False,"LGALS3_positive_mediator_supported":False,"C4_LGALS3_enrichment_supported":True,"prior_STAT3_SPP1_preserved":True,"new_MTA_interaction_models":0,"new_figures":0,"updated_at":stamp}
(d/"STATUS.json").write_text(json.dumps(status,indent=2))
bk=d/"root_record_backups";bk.mkdir(exist_ok=True)
entry=f"""
## Phase182 IL10来源/受体与LGALS3核查完成（{stamp}）
用户授权下一步核查，固定10基因面板；原75文库258120细胞来源审查，冻结精细髓系统计。IL10RA/RB共同RNA检出23患者达到>=10细胞，A2B仅3患者，后者不宜单细胞高低分组。IL10 RNA不默认肿瘤来源：保守恶性仅3/20患者达10检出，其他成熟巨噬29/34、单核22/29、T15/33；RNA不等于分泌。
LGALS3原不在两套STAT3靶集中，无须改评分。LGALS3—SPP1主31患者20正/rho0.04996/P=FDR0.02889，但中位bootstrap区间跨零；控制C4后rho0.01850/P=FDR0.21336。STAT3—LGALS3主33患者27负/rho-0.10947/P2.68e-5/FDR4.02e-5，DoRothEA亦无正向支持；不支持正向必经中介。
C4 LGALS3配对24/24升高，中位log2CPM差1.98519/P1.94e-5，保留定位价值。旧STAT3—SPP1逐患者结果复现；MTA/原两段全保留。2974核查全通过，无图、无新MTA交互/中介。
报告：{d}/IL10_LGALS3_AUDIT_RESULTS_CN.md
下一断点：背景依赖髓系扰动资源/响应资格，不把IL10 RNA当胞外暴露，不继续换评分追LGALS3中介；免疫空间暂缓。
"""
for n in ["PROJECT_STATUS.md","DECISIONS.md","CODEX_HANDOFF.md"]:
 p=r/n
 if not (bk/n).exists():shutil.copy2(p,bk/n)
 old=p.read_text()
 if "## Phase182 IL10来源/受体与LGALS3核查完成" not in old:p.write_text(entry+"\n"+old,encoding="utf8")
p=r/"RESULTS_MANIFEST.tsv"
if not (bk/p.name).exists():shutil.copy2(p,bk/p.name)
if "phase182_IL10_LGALS3" not in p.read_text():
 with p.open("a",newline="") as f:csv.writer(f,delimiter="\t").writerow(["phase182_IL10_LGALS3","L2_receiver_audit","report",str(d/"IL10_LGALS3_AUDIT_RESULTS_CN.md"),stamp,"PASS","2974 checks;old edge reproduced","IL10R measurable;LGALS3 mediation unsupported;C4 enrichment retained"])
de=d/"delivery";de.mkdir(exist_ok=True)
for n in ["IL10_LGALS3_AUDIT_RESULTS_CN.md","MANUSCRIPT_MODULE_CN.md","ALL_STATISTICS_CN.md","RECEIVER_COVERAGE_SUMMARY.tsv","SOURCE_PANEL_SUMMARY.tsv","EDGE_SUMMARY.tsv","SCORE_CIRCULARITY_AUDIT.json","VALIDATION_STATUS.json","STATUS.json"]:shutil.copy2(d/n,de/n)
hh={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in de.iterdir()};(d/"DELIVERY_SHA256.json").write_text(json.dumps(hh,indent=2));print(json.dumps({"status":status,"hashes":hh},ensure_ascii=False))
