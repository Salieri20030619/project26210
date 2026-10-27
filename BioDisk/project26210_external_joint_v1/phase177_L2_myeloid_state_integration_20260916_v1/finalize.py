from pathlib import Path
import csv,json,hashlib,statistics,math,datetime,shutil
import pandas as pd,numpy as np
r=Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1");d=r/"phase177_L2_myeloid_state_integration_20260916_v1"
shutil.copy2(r/"run_phase177.py",d/"run_analysis.py")
def read(n):
 with (d/n).open() as f:return list(csv.DictReader(f,delimiter="\t"))
checks=[]
def check(name,ok):checks.append({"check":name,"pass":bool(ok)})
for item in json.loads((d/"INPUT_MANIFEST.json").read_text()):
 check("input_unchanged_"+item["key"],hashlib.sha256(Path(item["source"]).read_bytes()).hexdigest()==item["sha256"])
vals=read("COMMON_PATIENT_SENSITIVITY_VALUES.tsv")
for row in read("COMMON_PATIENT_SENSITIVITY_SUMMARY.tsv"):
 z=[v for v in vals if v["feature"]==row["feature"] and v["comparison"]==row["comparison"]]
 label=row["feature"]+"_"+row["comparison"]
 check(label+"_unique_patients",len({v["patient"] for v in z})==len(z)==int(row["n_patients"]))
 for col in ["primary","sensitivity"]:
  vv=[float(v[col]) for v in z]
  check(label+"_"+col+"_median",math.isclose(statistics.median(vv),float(row[col+"_median"]),rel_tol=1e-12,abs_tol=1e-12))
  check(label+"_"+col+"_positive",sum(x>0 for x in vv)==int(row[col+"_positive"]))
 check(label+"_both_positive",sum(float(v["primary"])>0 and float(v["sensitivity"])>0 for v in z)==int(row["both_positive"]))
 check(label+"_same_sign",sum(float(v["primary"])*float(v["sensitivity"])>0 for v in z)==int(row["same_nonzero_sign"]))
for row in read("COMMON_PATIENT_LEAVE_STUDY.tsv"):
 z=[v for v in vals if v["feature"]==row["feature"] and v["comparison"]==row["comparison"] and v["study"]!=row["omitted_study"]]
 check("LOSO_"+row["feature"]+row["comparison"]+row["omitted_study"],len(z)==int(row["remaining_n"]) and all(math.isclose(statistics.median(float(v[c]) for v in z),float(row[c+"_median"]),abs_tol=1e-12) for c in ["primary","sensitivity"]))
assert all(c["pass"] for c in checks)
pd.DataFrame(checks).to_csv(d/"INDEPENDENT_VALIDATION.tsv",sep="\t",index=False)
def md(df):
 def fmt(v):
  if isinstance(v,(float,np.floating)):return f"{v:.6g}"
  return str(v)
 return "| "+" | ".join(df.columns)+" |\n|"+"|".join(["---"]*len(df.columns))+"|\n"+"\n".join("| "+" | ".join(fmt(v) for v in row)+" |" for row in df.itertuples(index=False,name=None))
dy=pd.read_csv(d/"TUMOR_DYNAMIC_EVIDENCE.tsv",sep="\t")
main=dy[(dy["mode"]=="primary")&(dy.pt_version=="pt_primary")&(dy.adjustment=="depth_adjusted_rho")]
sens=dy[(dy["mode"]=="without_stress_cycle")&(dy.adjustment=="depth_adjusted_rho")]
cols=["feature","n_patients","median_effect","positive","P","FDR"]
common=pd.read_csv(d/"COMMON_PATIENT_SENSITIVITY_SUMMARY.tsv",sep="\t")
report="""# L2第二步：SPP1髓系状态、轨迹与患者联系
## 结果一：SPP1与胆固醇处理程序沿可估计轨迹呈正向变化
以肿瘤内单核/巨噬细胞、患者内测序深度调整后的相关系数为重点，主Monocle3轨迹中SPP1为12/14患者正向，中位rho=0.04946，P=0.01435、原家族FDR=0.07177，属于名义显著。其效应较小，但方向具有患者一致性。胆固醇代谢程序12/14正向，rho=0.36499，P=0.01706、FDR=0.08203，同属名义显著；溶酶体程序10/14正向，主分析未达名义显著。

主分析完整重点统计：
"""+md(main[cols])+"""

去应激/周期敏感性中，SPP1为18/24正向，rho=0.05625、P=0.01842、FDR=0.03069；胆固醇代谢与溶酶体程序均获原家族FDR支持。该敏感性涵盖不同患者和细胞，不能直接用较小FDR取代主结果。
"""+md(sens[cols])+"""

## 结果二：固定共同患者后，正向状态变化仍然保留
本轮新增描述性复算，固定患者后比较主Monocle3与另两种排序；未新增推断性检验。
"""+md(common[["feature","comparison","n_patients","n_studies","primary_median","sensitivity_median","primary_positive","sensitivity_positive","both_positive"]])+"""

SPP1：主分析与去应激/周期共同13患者中，11人与9人正向、9人共同正向；主分析与DPT共同14患者中，12人与13人正向、11人共同正向。逐一留出研究后两种方法的SPP1中位方向仍为正；部分留出仅剩5–9患者，仅作描述。
完整DPT肿瘤集合32患者则为19人正向，rho=0.04220、P=0.12291、FDR=0.13657。因此，当前支持有覆盖范围限定的状态趋势，不能宣称全部患者/全部轨迹一致。共同患者也不意味着共同细胞，本轮没有重建相同细胞轨迹。

## 结果三：C4是SPP1集中的状态，不能写成唯一终末巨噬
主轨迹仅442/2674个C4细胞获得有效拟时序（16.53%）；去应激/周期后1743/2674（65.18%）。这是算法覆盖率，细胞数不是独立样本量。
肿瘤内C4相对单核：主分析8/8患者更晚，中位标准化拟时序差+0.55866，按小样本探索性呈现；敏感性19/19更晚，差+0.62201，P=0.0001430、原FDR=0.0002860。
相对其他成熟巨噬，主分析仅4/8正向；敏感性7/19正向、中位差-0.06937、P=0.30481、FDR=0.30481。当前数据不支持C4普遍是所有巨噬的终末阶段。

同文库16患者中，C4细胞比例中位数7.406%，其贡献的成熟巨噬SPP1原始RNA计数比例中位数54.590%。该结果支持C4的定位价值；并非分泌蛋白量、全组织总量，也不能把两个中位数相除作为每细胞倍数。
SPP1富集、胆固醇/溶酶体表型与整体轨迹同向变化是不同层面的证据；程序分别沿轨迹上升，不等于已检验其在每个细胞内与SPP1直接共调控。

## 结果四：同患者发送端桥接尚弱，整体SPP1不宜被C4比例替代
|16名同肿瘤文库患者的肿瘤TYMP对应终点|原始rho|P|原三终点FDR|研究内秩调整rho|置换P|FDR|
|---|---|---|---|---|---|---|
|C4比例|0.46765|0.06776|0.20328|0.30293|0.2851|0.8553|
|整体巨噬SPP1 pseudobulk|0.10588|0.69633|0.69633|-0.10810|0.7947|0.8553|
|整体巨噬SPP1检出比例|0.32941|0.21280|0.31920|0.08012|0.8553|0.8553|

旧17患者的平均SPP1表达关联rho=0.50735、P=0.03763、原FDR=0.18816，研究调整P=0.4324。其患者、表达聚合和文库口径不同，仍可保留为旧口径名义信号，但不能覆盖同文库16患者结果。

既有巨噬侧分析中，控制研究及C4比例后，LDHA—SPP1 rho=0.81778（P=0.001155、原选择性敏感性家族FDR=0.001733），SLC16A3—SPP1 rho=0.84003（P=0.0006266、FDR=0.001733）。这支持整体巨噬SPP1还有不由C4比例完全解释的代谢状态联系。它不是MTA响应证据，也不自动把LDHA/SLC16A3变成TYMP的中介。

## 关键反向结果：不恢复三TF共同驱动SPP1的叙事
31名肿瘤患者中，深度调整后JUN/JUNB/CEBPD活性与SPP1分别28/31、31/31、29/31负向，中位rho=-0.09865/-0.11540/-0.16830，原三TF家族FDR均通过。完整P/FDR见SPP1_TF_EVIDENCE.tsv。
髓系侧TYMP沿轨迹多为负向，不等于肿瘤发送端TYMP降低；两种细胞区室必须分开。轨迹伴随变化与TF—SPP1关系也不得互相代替。

## Methods
复用冻结患者内动态效应、Monocle3分区覆盖、原桥接表、巨噬组成敏感性及C4贡献表，保留原始患者资格和统计家族。主叙事优先肿瘤内、深度调整、pt_primary；原表保存全部组织/未调整/根选择/DPT版本。相同患者交集内计算中位效应、正向人数和共同方向，并按研究逐一留出计算描述性中位数。新计算不生成P值，不重算FDR。
从逐患者表复核400条历史动态摘要的n、中位效应、正向数和负向数，共1600项一致性检查全部通过。10份输入SHA256记录并复查未变；共同患者和留研究结果使用标准库statistics独立复算。程序基因表检查TYMP/SPP1是否在内，见PROGRAM_TARGET_GENE_AUDIT.tsv。
本轮不重建轨迹、不调整分群、不生成新图、不使用客户湿实验数据作公共统计。

## Discussion与正文连接
现有数据可以组织为：第一部分TYMP相关代谢发现引出MTA候选；外部直接髓系扰动提供研究髓系状态的理由；本项目建立SPP1相关状态变化与C4定位；后续分别研究免疫功能和空间结构。无需要求每段均有强因果证据，但不能将“提出待验证联系”写成“MTA已驱动该状态”。

当前最需要实验补充的是MTA暴露与SPP1/对应髓系程序之间的关系，而非再次证明C4存在。客户现有图仅包含肿瘤细胞代谢变化，不能补齐这一箭头。
拟时序是表达排序，不能证明单核谱系转化或真实时间；分群后标志和状态统计属于同数据条件性发现。外部MTA实验与本队列不是同一患者证据，M2标志不等于C4身份。主要阴性已在结果保留。

## 当前完成范围与下一步
完成L2当前可用数据的来源—状态—轨迹—患者关联整合和共同患者敏感性。未建立MTA直接诱导SPP1机制，未运行两条下游终点。
接下来先进入免疫功能分支：核对整体巨噬SPP1与T细胞状态/抗原呈递的患者级关联及细胞组成；空间分支独立核对SPP1相关状态与血管、邻域及组织边界，旧SPP1-state不可直接改名新C4。两条分支分别设终点和校正家族，功能结果不代替空间结果。
"""
(d/"L2_STATE_INTEGRATION_RESULTS_CN.md").write_text(report,encoding="utf8")
module="""# 可拆入统一正文的L2模块草稿
## SPP1相关髓系状态变化与C4定位形成相互补充的证据
在肿瘤内单核/巨噬细胞的可估计轨迹中，SPP1随状态排序呈患者一致的正向变化：深度调整后12/14患者为正，中位rho=0.0495（P=0.0144，FDR=0.0718），为名义显著。去应激/周期敏感性支持相同方向；固定共同患者后正向趋势保留。胆固醇代谢程序也呈正向变化，支持进一步考察SPP1相关的代谢适应状态。

冻结分群中，C4具有SPP1富集及脂质处理/溶酶体相关特征。同文库16患者中，C4占成熟巨噬的比例中位数为7.4%，所捕获SPP1原始RNA计数贡献中位数为54.6%，说明其是SPP1信号集中的状态。可估计轨迹中C4位于单核之后，但当前数据不支持其普遍晚于其他成熟巨噬。

同文库患者层分析尚未显示肿瘤TYMP与整体巨噬SPP1的稳定关联；三TF共同正向驱动SPP1的模型也不受当前数据支持。因此，结合外部MTA髓系作用依据，本研究提出MTA相关代谢环境参与SPP1髓系状态调节的待验证联系，后续分别考察该状态的免疫功能联系与空间组织特征。

使用说明：与第一部分结果拆分整合，共用Methods/Discussion和最终图号。主分析覆盖、效应、名义显著与关键阴性不得删除；全部敏感性和统计保留于本阶段审计报告。
"""
(d/"L2_MANUSCRIPT_MODULE_CN.md").write_text(module,encoding="utf8")
stamp=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
status={"status":"COMPLETE_CURRENT_DATA_L2_STATE_INTEGRATION","historical_summary_rows":400,"historical_numeric_checks":1600,"new_independent_checks":len(checks),"all_checks_pass":True,"new_p_values":0,"new_figures":0,"MTA_SPP1_causality_established":False,"next":"immune endpoint branch; spatial endpoint separately","updated_at":stamp}
(d/"STATUS.json").write_text(json.dumps(status,indent=2),encoding="utf8")
bk=d/"root_record_backups";bk.mkdir(exist_ok=True)
entry=f"""
## Phase177 L2髓系状态整合完成（{stamp}）
肿瘤内主Monocle3深度调整SPP1 12/14正、rho0.04946、P0.01435/FDR0.07177名义；去应激周期18/24正/FDR0.03069。共同13患者主/敏感性11/9正，共同14患者主/DPT12/13正；留研究中位方向保留，均为描述性、不新增P。
C4主轨迹覆盖16.53%、敏感性65.18%；晚于单核不等于所有成熟巨噬终末。C4细胞中位7.406%、SPP1原始RNA贡献54.590%复现。
16同文库TYMP—SPP1表达rho0.10588/P0.69633；旧17例名义信号不能替代。三TF与SPP1主要负向保留。MTA→SPP1仍是待验证联系，不把候选叙事写成机制完成。
400历史摘要/1600数值检查及{len(checks)}新增独立检查全通过，10输入未改，无新图/P/FDR。
报告：{d}/L2_STATE_INTEGRATION_RESULTS_CN.md
正文模块：{d}/L2_MANUSCRIPT_MODULE_CN.md
当前L2可用证据整合完成。下一步免疫功能分支患者级复核，空间分支单独处理，不改旧状态为C4。
"""
for name in ["PROJECT_STATUS.md","DECISIONS.md","CODEX_HANDOFF.md"]:
 p=r/name
 if not (bk/name).exists():shutil.copy2(p,bk/name)
 old=p.read_text()
 if "## Phase177 L2髓系状态整合完成" not in old:p.write_text(entry+"\n"+old,encoding="utf8")
p=r/"RESULTS_MANIFEST.tsv"
if not (bk/p.name).exists():shutil.copy2(p,bk/p.name)
if "phase177_L2_state_integration" not in p.read_text():
 with p.open("a",newline="") as f:csv.writer(f,delimiter="\t").writerow(["phase177_L2_state_integration","L2","report",str(d/"L2_STATE_INTEGRATION_RESULTS_CN.md"),stamp,"PASS","1600 historical + "+str(len(checks))+" independent checks","No new P/FDR; no figures; MTA-SPP1 still hypothetical"])
delivery=d/"delivery";delivery.mkdir(exist_ok=True)
for n in ["L2_STATE_INTEGRATION_RESULTS_CN.md","L2_MANUSCRIPT_MODULE_CN.md","COMMON_PATIENT_SENSITIVITY_SUMMARY.tsv","COMMON_PATIENT_LEAVE_STUDY.tsv","TUMOR_DYNAMIC_EVIDENCE.tsv","SPP1_TF_EVIDENCE.tsv","STATUS.json"]:shutil.copy2(d/n,delivery/n)
hashes={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in delivery.iterdir()}
(d/"DELIVERY_SHA256.json").write_text(json.dumps(hashes,indent=2))
print(json.dumps({"status":status,"hashes":hashes},ensure_ascii=False))
