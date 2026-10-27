
from pathlib import Path
import pandas as pd,json,datetime,hashlib,shutil
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=P/'phase165_A_input_confirmation_20260915_v1'
d=pd.read_csv(O/'SECTION_INPUT_QC.tsv',sep='\t');c=pd.read_csv(O/'CANDIDATE_COVERAGE.tsv',sep='\t');t=pd.read_csv(O/'PROVISIONAL_TECHNICAL_PILOT.tsv',sep='\t')
m=json.loads((O/'metadata/SAMPLE_METADATA_EXTRACT.json').read_text())
pd.DataFrame([dict(section=x['gsm'],geo_title=x['title'][0],patient_id_confirmed='',mapping_status='UNRESOLVED',potential_duplicate='c_2 across preparations' if x['gsm'] in ['GSM5924030','GSM5924050'] else '') for x in m]).to_csv(O/'PATIENT_MAPPING_REVIEW.tsv',sep='\t',index=False)
checks=['coordinates_finite','coordinates_unique','barcodes_unique','abundance_exact_order','abundance_finite','counts_nonnegative_integer','total_counts_match','all_candidates_present']
ok=bool(d[checks].all().all() and len(d)==24)
report=f"""# A 方案确认记录
日期：2026-09-15。范围：A1 输入验收与 A2 设计确认；本轮未运行空间关联模型、未重做解卷积、未绘图。

结论：{'基础技术输入通过，可进入小规模技术试点；患者映射与细胞归属仍需补充，尚不具备宣称全部方案验收通过的条件。' if ok else '存在基础 QC 未通过项，见明细；暂不启动试点。'}

## 输入验收
24 张切片，FFPE 与 frozen 各 12 张；共 {int(d.spots.sum())} 个 spot。FFPE 17,943 基因，frozen 36,601 基因。
计数层非负整数、测序量核对、坐标有限且唯一、条码唯一及与细胞丰度表顺序一致：{str(ok)}。
TYMP、ODC1、SRM、SMS、MTAP、AMD1、MAT2A、SLC16A3、SPP1 及 COL1A1/COL1A2/DCN/LUM 在各切片的覆盖和检出数见 CANDIDATE_COVERAGE.tsv。基因可检测不等于关系成立。
现有空间坐标为图像像素；未经尺度核验不报告微米距离。

## 患者身份
已归档 GEO 官方 SOFT 与 24 个 GSM 的标题、特征；未提供明确患者 ID。
ffpe_c_2 与 frozen_c_2 可能来自同一患者，不能据命名直接确认为同一患者，也不能视为两个独立患者。
A/B/C 不可直接视为三个患者。官方分析脚本的旧 slide_list 与 GEO 名称也不完全一致，不能机械替换。
当前独立患者 n 未确认；正式跨切片推断暂停，技术试点可做描述性验证。
后续优先查原文补充临床表或作者样本映射；未解决时仅报告逐切片效应，不提供伪患者 P 值。

## 参考与终点修正
现有 13 状态参考不含新 C4，也不含单列成纤维细胞。
总体成熟巨噬丰度定义为 C1QC-state Macrophage + Other Macrophage + SPP1-state Macrophage，单核单列。
旧 SPP1-state 不改名为新 C4。C4 保留既有单细胞胆固醇处理状态定位，空间投射须另做验证。
基质评分只作为 RNA 表达代理，不写成成纤维细胞比例。
spot 中 TYMP/SPP1 不可直接归属恶性细胞/巨噬细胞。初步空间终点表述为“恶性细胞富集区域的候选表达与巨噬富集邻域 SPP1”，细胞来源由既有单细胞结果与后续参考敏感性检验辅助。
既有解卷积来自同一 RNA，需审核 TYMP/SPP1/候选基因参与参考的情况，必要时用剔除目标基因的参考作敏感性验证；不得将同源评分间相关当独立验证。

## 冻结的分阶段计划
1. 技术试点：每种制备选 spot 数距组内中位数最近的切片，平局按 GSM 排序，完全不按相关系数或 P 选片。预选 {', '.join(t.section)}。只验证邻接、边界、检测覆盖与模型可识别性；这两张切片不承担患者推断。
2. 主候选：SRM、SMS（分别检验，不拼成已成立中介）；SLC16A3 并行。ODC1/AMD1/MAT2A 作代谢背景，MTAP 单列清除背景，不预设缺失。
3. 上游问题为 TYMP—候选表达，接收端问题为候选邻域—SPP1。MTA 未实测，不生成“预测 MTA 浓度”充当验证。既有 MTA 串联阴性和分段正向结果均保留。
4. 邻域预设六个最近非自身 spot；先检查孔洞、断裂和异常长边，再冻结距离截断。主关系用未平滑表达；GraphPCA 仅作结构/敏感性选项，不为提高相关而选参数。区域阈值和边长阈值须在技术试点后、查看关联结果前写入最终预注册。
5. 调整层次预设为原始、测序深度及恶性/总体巨噬组成、再加基质/缺氧表达背景。明确各层效应变化；避免小区域过多协变量。检查共线性和有效空间信息后锁定模型。
6. 正式检验以独立患者为单位，FFPE/frozen 分层，同患者重复不得重复计 n。上游三候选和下游三候选分别形成同类 FDR 家族；报告效应、P、FDR及患者一致性。n<10 按探索性；无有效生物学重复只报告效应。
7. 与 Phase75/78 旧空间结果及 Phase163 串联模型逐项比较；新方法产生的差异须说明来自终点、区域、参考还是控制方式。先审统计文字，再决定扩展和绘图。

## 当前尚未执行
未进行技术试点邻接计算、空间关联、患者汇总、GraphPCA、SpaMTP 安装或新 C4 解卷积。
确认结论不等于 TYMP—MTA—SPP1 已获空间支持。

## 来源与明细
- GEO：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175540
- 官方脚本：https://raw.githubusercontent.com/maximemeylan/Meylan_et_al_2022/main/00_visium_preprocessing.R
- SECTION_INPUT_QC.tsv、CANDIDATE_COVERAGE.tsv、PATIENT_MAPPING_REVIEW.tsv、PROVISIONAL_TECHNICAL_PILOT.tsv。
"""
(O/'A_CONFIRMATION_CN.md').write_text(report)
status=dict(stage='Phase165',status='TECHNICAL_INPUT_PASS_PATIENT_MAPPING_PENDING' if ok else 'QC_FAILED',sections=len(d),spots=int(d.spots.sum()),technical_qc_pass=ok,patient_n=None,patient_mapping_complete=False,pilot_selected=list(t.section),pilot_executed=False,association_models_run=0,new_figures=0)
(O/'STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
stamp=datetime.datetime.now().astimezone().isoformat();back=O/'authority_backups';back.mkdir(exist_ok=True)
note=f"\n## Phase165 A方案确认（{stamp}）\n基础输入 QC：{ok}；24切片，患者独立性未确认，旧13状态不含新C4或成纤维细胞。试点预选但未执行，无新关联或图。详见 {O}/A_CONFIRMATION_CN.md。\n"
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 f=P/name
 if f.exists():
  shutil.copy2(f,back/name);f.write_text(note+'\n'+f.read_text())
f=P/'RESULTS_MANIFEST.tsv'
shutil.copy2(f,back/f.name)
with f.open('a') as h:
 for name in ['A_CONFIRMATION_CN.md','STATUS.json','SECTION_INPUT_QC.tsv','CANDIDATE_COVERAGE.tsv','PATIENT_MAPPING_REVIEW.tsv','PROVISIONAL_TECHNICAL_PILOT.tsv']:
  h.write('\t'.join(['phase165_'+name,'Phase165','confirmation',str(O/name),stamp,'TECHNICAL_PASS_IDENTITY_PENDING', '24 section counts coordinates barcodes and gene coverage audited','No association tests; patient mapping unresolved'])+'\n')
print(json.dumps(status,ensure_ascii=False));print('REPORT_SHA256',hashlib.sha256((O/'A_CONFIRMATION_CN.md').read_bytes()).hexdigest())

