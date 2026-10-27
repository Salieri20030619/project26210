from pathlib import Path
import datetime,hashlib,json,shutil,csv
root=Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1")
d=root/"phase176_L2_MTA_myeloid_bridge_20260916_v1"
stamp=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
report="""# L2第一步：MTA—髓系—SPP1证据核查
## 结论
已有直接实验支持MTA改变单核/巨噬细胞功能状态；本项目已识别具有脂质处理和溶酶体表达特征的SPP1富集C4。两项事实可以构成研究问题的逻辑衔接，当前证据不能将其合并为已确认的“MTA诱导SPP1/C4”机制。

上游按用户决定暂时收束于Phase175：代谢背景和分段支持保留，SMS及多胺短路径未获稳定支持的结果保留，不再更换上游评分追求阳性。第一部分仍是生信发现起点，客户湿实验承担后续验证。

## 直接作用来源
| 来源 | 实验对象与设计 | 核查结果 | 本阶段用途 |
|---|---|---|---|
| Hansen 2022 | RAW264.7、BV2、PMA预处理THP-1及M-CSF培养的人原代单核来源细胞；外源MTA | ARG1/VEGFA等改变；含A2B和STAT3干预 | 支持MTA可作用于髓系，不能等同SPP1诱导 |
| Keyel 2014 | 小鼠BMDM和RAW报告细胞，MTA联合TLR刺激 | 炎症与共刺激响应改变，有刺激背景依赖 | 支持功能状态可塑性 |
| GSE249082 | 人成纤维细胞对照与MTA处理，8个文库 | 测序对象为成纤维细胞 | 不作为直接巨噬MTA转录组 |

Hansen实验多数RNA读出为12小时；原代细胞方法写M-CSF培养48小时后加MTA24小时，图2J图注却写12小时，保留该处差异。图2J列100/200 μM；其他实验剂量不可直接沿用。图注n≥4或n≥5为实验重复，原代独立供体数未明确，不能按供体样本量解读。作者未可靠测到所测时间点STAT3 Tyr705磷酸化变化，不能改写成已证明MTA升高pSTAT3。
[原始论文](https://www.nature.com/articles/s41598-022-07697-0)

Keyel主要功能实验200–500 μM、过夜；部分转录/NF-κB实验4小时，图示3–5次实验。TNF、CD86等响应降低；IL10依赖刺激条件，IL1β并非同步降低。没有据此确认SPP1。
[原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC4130577/)

GSE249082研究涉及成纤维细胞代谢后影响巨噬的旁分泌过程，但其公开RNA测序对比不能改称巨噬直接扰动。既有Phase77跨环境签名未通过预设资格，继续保留原判定。
[GEO官方记录](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249082)

## 项目实际表型及衔接
C4并非按SPP1阈值预选：冻结分群后多基因鉴定为GPNMB/CCL18相关、SPP1富集巨噬。26名有效配对患者中24名SPP1升高，中位log2(CPM+1)差1.448645，P=6.6778e-05，原基因家族FDR=0.000544714。五项研究的中位方向均正，属于内部一致性。

C4胆固醇处理和溶酶体相关程序升高，抗原呈递、炎症及干扰素相关程序较低。可以据此描述实际表达状态，不直接写成已测定免疫抑制或吞噬功能。来源为冻结的part2_rebuild_v1报告与统计。

Phase163受体面板以24名有效患者比较C4与其他巨噬：
- ADORA2B在C4中位检测比例为0，双方均无患者达到10%细胞检出；ARG1未检出。
- IL10在21/24患者较低，差值中位数-0.609921，P=4.42266e-05，原面板FDR=7.73966e-05。
- VEGFA在23/24患者较低，差值中位数-0.865160，P=5.96046e-07，FDR=1.39078e-06。
- STAT3与SOCS3均在23/24患者较低；完整7基因结果保存在随附原表。

因此，当前数据不支持C4复现该文献整套转录标志。RNA低检出不能证明受体蛋白不存在，STAT3 RNA不能代表其活性；上述差异也不能用来否定所有MTA髓系作用。M2、整体SPP1状态、新C4和旧SPP1-state保持区分。

## 当前可采用的叙事
“TYMP相关代谢结果将MTA确立为值得进一步研究的候选代谢物。已有直接扰动研究表明MTA能够改变髓系细胞功能状态，因此进一步考察本项目单核/巨噬细胞的状态组成。患者层分析显示，SPP1富集巨噬具有一致的脂质处理和溶酶体相关表达特征，由此提出MTA相关代谢环境与该状态形成之间的待检验联系。”

这段是候选机制的推理衔接，不把外部MTA实验与本项目C4观察写成同一实验的因果结果。尚缺直接MTA暴露后SPP1及对应程序的合格髓系扰动证据。

## 方法与本轮完成范围
对两篇原始作用论文及GEO记录进行有界核查，服务器保存原始HTML/XML/SOFT。补充检索词包括methylthioadenosine macrophage SPP1/osteopontin及macrophage RNA-seq/GEO。检索命中的大鼠颅骨分化研究不属于巨噬细胞，排除直接髓系证据用途。本轮未找到可确认合格的直接巨噬MTA转录组，不宣称全球不存在该数据。
复用冻结患者统计，不新作推断检验、不重新计算FDR、不重聚类、不绘图。此次完成L2来源及表型对照第一步，不代表整个L2与两条下游分支完成。

## 下一分析断点
1. 从冻结整体单核/巨噬结果整理SPP1连续状态、单核到巨噬的轨迹及程序变化；轨迹只作为状态连续性，不能代替谱系追踪。
2. 以整体巨噬为功能分析主体、C4为定位对象；复核已有患者配对及研究调整结果，先确认哪些关联独立于组成和研究差异。保留关键阴性，不强制三TF或A2B—STAT3作中介。
3. 只有找到对象正确、对照和生物重复合格的直接扰动数据才建立外部响应签名；当前不移植成纤维签名冒充巨噬证据。
4. 完成髓系状态衔接后分别进入免疫功能与空间组织两条分支，各自定义终点、统计家族和结论。
"""
(d/"L2_BRIDGE_RESULTS_CN.md").write_text(report,encoding="utf8")
src=root/"phase163_polyamine_MTA_serial_receiver_20260915_v1/MTA_RECEIVER_PANEL_SUMMARY.tsv"
shutil.copy2(src,d/"MTA_RECEIVER_PANEL_SUMMARY_UNCHANGED.tsv")
assert (d/"MTA_RECEIVER_PANEL_SUMMARY_UNCHANGED.tsv").read_bytes()==src.read_bytes()
checks=[]
for f in sorted((d/"sources").iterdir()):
 checks.append({"file":str(f.relative_to(d)),"bytes":f.stat().st_size,"sha256":hashlib.sha256(f.read_bytes()).hexdigest()})
(d/"SOURCE_CHECKSUMS.json").write_text(json.dumps(checks,indent=2),encoding="utf8")
status={"phase":176,"status":"L2_SOURCE_PHENOTYPE_REVIEW_COMPLETE","full_L2_complete":False,"new_models":0,"new_figures":0,"source_downloads":3,"receiver_table_byte_identical":True,"updated_at":stamp}
(d/"STATUS.json").write_text(json.dumps(status,indent=2),encoding="utf8")
backup=d/"root_record_backups";backup.mkdir(exist_ok=True)
entry=f"""
## Phase176 上游暂收束，L2来源及表型核查完成（{stamp}）
用户“先暂时这样，记录后下一步”：Phase175代谢背景/分段支持保留，不继续更换SMS或多胺模型；不确认中介。
已核对Hansen2022、Keyel2014与GSE249082并保存3份官方来源。MTA直接髓系功能作用有依据，SPP1/C4诱导未建立；Hansen含人源实验，不能仅写鼠实验。其原代供体数未明确，方法/图注时间差异保留。
冻结C4患者SPP1富集成立；IL10/VEGFA等与外部MTA标志不一致，低RNA检出不证明蛋白或活性缺失。旧成纤维MTA签名不冒充巨噬扰动。
报告：{d}/L2_BRIDGE_RESULTS_CN.md
本轮仅来源核查和已有表型对照，无新增模型、P/FDR、聚类或图。L2整体尚未完成。下一断点为整体单核/巨噬SPP1状态、轨迹与患者关联整理；免疫和空间后续仍分支。
"""
for name in ["PROJECT_STATUS.md","DECISIONS.md","CODEX_HANDOFF.md"]:
 p=root/name
 b=backup/name
 if not b.exists():shutil.copy2(p,b)
 old=p.read_text(encoding="utf8")
 if "## Phase176 上游暂收束" not in old:p.write_text(entry+"\n"+old,encoding="utf8")
p=root/"RESULTS_MANIFEST.tsv"
if not (backup/p.name).exists():shutil.copy2(p,backup/p.name)
if "phase176_L2_source_review" not in p.read_text():
 with p.open("a",encoding="utf8",newline="") as f:
  csv.writer(f,delimiter="\t").writerow(["phase176_L2_source_review","L2_source_phenotype_review","report",str(d/"L2_BRIDGE_RESULTS_CN.md"),stamp,"REVIEW_COMPLETE","3 archived sources; unchanged receiver table","Bounded source review; full L2 pending; no new models"])
delivery=d/"delivery";delivery.mkdir(exist_ok=True)
for n in ["DECISION_AND_SCOPE_CN.md","L2_BRIDGE_RESULTS_CN.md","MTA_RECEIVER_PANEL_SUMMARY_UNCHANGED.tsv","SOURCE_CHECKSUMS.json","STATUS.json"]:shutil.copy2(d/n,delivery/n)
receipt={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in delivery.iterdir()}
(d/"DELIVERY_SHA256.json").write_text(json.dumps(receipt,indent=2),encoding="utf8")
print(json.dumps({"status":status,"delivery":receipt},ensure_ascii=False))
