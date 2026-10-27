from pathlib import Path
import json,re,hashlib,datetime,shutil,csv,io
import pandas as pd,numpy as np,h5py
O=Path(__file__).resolve().parent;ROOT=O.parent
now=datetime.datetime.now().astimezone().isoformat()
s=pd.read_csv(O/'SECTION_ENDPOINT_SUMMARY.tsv',sep='\t');d=pd.read_csv(O/'DESIGN_ENDPOINT_ELIGIBILITY.tsv',sep='\t')
v=pd.read_csv(O/'VALIDATION.tsv',sep='\t');w=pd.read_csv(O/'RELOAD_VALIDATION.tsv',sep='\t')
assert v.passed.all() and w.passed.all()
meta=json.loads((ROOT/'phase165_A_input_confirmation_20260915_v1/metadata/SAMPLE_METADATA_EXTRACT.json').read_text())
rows=[]
for x in meta:
 title=x['title'][0];m=re.search(r'(ffpe|frozen)_([abc]_\d+)',title)
 assert m is not None
 rows.append(dict(section=x['gsm'],geo_title=title,preparation_from_geo='FFPE' if m[1]=='ffpe' else 'frozen',candidate_label_not_patient_id=m[2],patient_id_confirmed='',mapping_status='UNRESOLVED',source='GEO sample metadata + original author repository review'))
pm=pd.DataFrame(rows)
pm['label_shared_across_preparations']=pm.candidate_label_not_patient_id.duplicated(keep=False)
assert len(pm)==24
assert pm.merge(s,on='section').eval('preparation_from_geo == preparation').all()
pm.to_csv(O/'PATIENT_MAPPING_REVIEW.tsv',sep='\t',index=False)
tree=json.loads((O/'metadata/author_repository_tree.json').read_text())
author=(O/'metadata_api/00_visium_preprocessing.R').read_text()
block=author.split('slide_list <- c(',1)[1].split(')',1)[0]
labels=re.findall(r'"([^"]+)"',block)
pd.DataFrame(dict(author_script_slide_label=labels,exact_label_present_in_geo=[x in set(pm.candidate_label_not_patient_id) for x in labels])).to_csv(O/'AUTHOR_GEO_LABEL_COMPARISON.tsv',sep='\t',index=False)
coverage=s.groupby('preparation').agg(sections=('section','size'),source_eligible=('source_all_designs_eligible','sum'),receiver_eligible=('receiver_all_designs_eligible','sum'))
with h5py.File(ROOT.parent/'project26210_external_joint_v1_phase54_gse175540_spatial_v1/run_v1/sections_v1/GSM5924034/processed_section_v1.h5ad') as f:dtype=str(f['layers/counts/data'].dtype)
audit=dict(author_repository='https://github.com/maximemeylan/Meylan_et_al_2022',commit=tree['sha'],files_reviewed=8,issues_reviewed=5,geo_samples=24,confirmed_patient_ids=0,repeated_candidate_label='c_2',repeated_sections=['GSM5924030','GSM5924050'],author_slide_labels=labels,exact_labels_shared_with_geo=sorted(set(labels)&set(pm.candidate_label_not_patient_id)),conclusion='No explicit GSM-to-independent-patient crosswalk established from retrieved GEO metadata and author repository; no inferred patient IDs assigned',publisher_page_access='Web retrieval returned Internal Error; supplementary clinical table not obtained',original_count_dtype=dtype)
(O/'PATIENT_SOURCE_REVIEW.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2))
issues=json.loads((O/'metadata_api/issue_2_comments.json').read_text())
assert any(x['user']['login']=='maximemeylan' and 'low tumor content' in x['body'] for x in issues)
report=f"""# 24张切片检查完成：24张上游端点、20张邻域接收端技术可用

## 结论
按Phase166严格规则逐片检查全部24张切片。上游TYMP—SRM/SMS/SLC16A3的设计在24/24张均可识别；源区候选—邻域SPP1接收端在20/24张满足冻结技术门槛。线性及逐列秩变换、基础及扩展调整的资格判定一致。尚未计算这些关系的效应量、关联检验、P/FDR或中介模型，未绘图。

| 制备 | 检查切片 | 上游端点可用 | 邻域接收端可用 |
|---|---:|---:|---:|
| FFPE | 12 | 12 | 9 |
| frozen | 12 | 12 | 11 |
| 合计 | 24 | 24 | 20 |

## 不合格接收端完整保留
| 切片 | 制备 | 合格接收邻域对应源spot数 | 主要原因 |
|---|---|---:|---|
| GSM5924034 | FFPE | 2 | 少于50，设计不满秩 |
| GSM5924038 | FFPE | 3 | 少于50，设计不满秩 |
| GSM5924039 | FFPE | 11 | 少于50；扩展设计不满秩 |
| GSM5924051 | frozen | 17 | 少于50；SRM线性扩展设计条件数32.48 |

这4张保留上游端点，不降低门槛挽救接收端。原作者对ffpe_c_10（GSM5924034）的公开回复说明该片肿瘤含量低、无TLS；这与技术覆盖不足相容，但不是本项目关联证据。原作者回复：https://github.com/maximemeylan/Meylan_et_al_2022/issues/2 。

## 设计与空间结构
24张共576个技术设计组合（24×2端点×2调整×3候选×2变换），528个可用、48个不适用；这里的组合数不是独立统计样本量。所有设计的结局均有非零方差，可用设计均满秩。可用上游条件数最大：线性3.9131、秩变换4.0674；接收端：线性13.2651、秩变换9.6905。秩变换未新增排除切片。

全部切片的距离筛选与Visium第一圈拓扑一致，邻接对称、无自身连接、最多6邻居。GSM5924048保留5个连通区：两个主要区域894和646 spot，接收端源spot分别238和113；另4个spot分布在3个小区。不强行跨区连接。后续应保留连通区标识，并补充分区内效应/分区调整敏感性；本阶段未将分区加入新的主模型。

接收邻域存在共享spot，复用次数最高6。邻域重叠详见NEIGHBOR_OVERLAP.tsv；不把源spot、接收spot或重复邻域当独立患者，也不将设计满秩解释为空间独立性。没有病理图像边界验收。

## 患者映射仍未解决
已复核GEO原始SOFT、24条样本信息、原作者仓库全部8个代码/README文件以及5个既有issue和回复。GEO标题可以确认12 FFPE/12 frozen，与输入清单完全一致，但未建立GSM至独立患者的明确对照表。ffpe_c_2与frozen_c_2仍是潜在重复个体线索，不据此宣称已经确认配对或23位患者。

原作者预处理脚本列出的12个slide标签与当前GEO标签仅部分一致；没有可审计重命名对照，不能依靠名称相似性补齐患者ID。原文出版社页面访问失败，未取得补充临床表，不能声称已穷尽原文所有患者信息。
- GEO：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175540
- 原作者仓库：https://github.com/maximemeylan/Meylan_et_al_2022 （锁定提交 {tree['sha']}）
- 原文：https://doi.org/10.1016/j.immuni.2022.02.001

患者映射表中确认ID留空。当前允许按切片描述效应、FFPE/frozen分层展示；不能将24张或20张切片作为24位或20位患者计算合并P/FDR。

## Methods
使用Phase165输入清单内全部24张processed_section_v1.h5ad、原始counts、spatial坐标及Phase78冻结13状态q05丰度。每张切片独立读取，不修改输入。18个基因均存在，表达为log1p(counts/spot总counts×10000)，本轮显式使用float64。合成候选SRM/SMS与SLC16A3分别检查；MTAP保留清除背景。

源区为归一化Malignant q05≥0.2且不低于每个其他类别与合并成熟巨噬；接收区合并C1QC-state/Other/SPP1-state成熟巨噬q05≥0.1且不低于本片中位数，单核单列。至多六邻居，距离≤本片最近邻中位距离×1.3，且符合Visium阵列第一圈；不含自身，每个源spot至少3个合格接收邻居。

上游设计以TYMP为预测变量，分别以SRM/SMS/SLC16A3为结局；接收设计以候选为预测变量、邻居SPP1均值为结局，调整源区TYMP。基础调整源区深度、恶性/合并巨噬丰度；接收端追加邻域对应变量。扩展增加基质COL1A1/COL1A2/DCN/LUM及缺氧CA9/EGLN3/VEGFA/BNIP3/PDK1的全片标准化均值，并对接收端追加邻域评分。分别检查原尺度标准化和端点区域内逐列平均秩转换后的矩阵。冻结门槛：≥50源spot、满列秩、条件数<30；额外核查有限值与结局非零方差。输出共享接收spot复用和连通区清单。

## 验证与执行记录
初始输入/图结构240/240通过，输出重载及输入完整性232/232通过，共472/472。52项输入与来源文件执行前后SHA256一致。两片历史试点邻接与全部spot成员完全一致，24个旧线性设计的秩与资格判定完全一致。
历史试点counts为{dtype}，当前显式float64，环境numpy版本亦与旧记录不同；条件数最大相对差9.0648e-7。严格1e-7数值容差曾失败，核对后用float32相容rtol=1e-5复核通过，未改变区域或资格阈值。保留PILOT_DESIGN_COMPARISON.tsv及旧验证脚本。首版验证脚本曾将transform列误写为同名方法属性，已用显式列索引修复；均未影响分析输出。
raw.githubusercontent.com下载过慢，已停止本轮慢下载助手，改用GitHub官方blob API获得同一提交的8文件；无作者联系、无新issue。软件版本与执行路径见COMPUTE_STATUS.json。

## 后续可执行范围
技术可用清单已冻结。下一分析可对24片上游及20片接收端计算逐片方向和效应量，并保留制备分层、区域结构及全部反向结果；患者独立性未确认前不做患者层汇总检验。旧参考包含目标基因，正式生物学解释仍需目标基因剔除敏感性。旧13状态没有新C4、未单列成纤维细胞：不将旧SPP1-state改名C4，也不将基质评分等同成纤维比例。MTA多胺轴、MTAP背景、SLC16A3并行和整体巨噬/C4分工保持。此次检查未新增支持或反对TYMP—X—SPP1的生物学关联结果。
"""
(O/'ENDPOINT_CHECK_CN.md').write_text(report,encoding='utf-8')
status=json.loads((O/'COMPUTE_STATUS.json').read_text())
status.update(status='PASS_TECHNICAL_CHECK_PATIENT_MAPPING_UNRESOLVED',source_sections_eligible=24,receiver_sections_eligible=20,FFPE_receiver_eligible=9,frozen_receiver_eligible=11,design_combinations=576,eligible_design_combinations=int(d.eligible.sum()),patient_mapping='UNRESOLVED',confirmed_patient_ids=0,validation_checks=len(v)+len(w),validation_passed=int(v.passed.sum()+w.passed.sum()),completed_at=now)
(O/'STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
(O/'EXECUTION_NOTES.json').write_text(json.dumps(dict(command='OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 '+status['executable']+' -u '+str(O/'run_check.py'),validation_command='OPENBLAS_NUM_THREADS=1 '+status['executable']+' '+str(O/'validate_check_v3.py'),first_validator_failure='pandas transform attribute collides with method; corrected explicit column indexing',second_validator_failure='condition-number relative tolerance 1e-7 too strict for float32 historical vs float64 current computation; max relative difference 9.064785485596893e-7; use rtol1e-5, retain exact graph/membership/rank/eligibility requirements',ssh='Intermittent connection close; read/write attempts retried; no computation failure'),ensure_ascii=False,indent=2))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
hashes=[]
for p in sorted(O.rglob('*')):
 if p.is_file() and 'authority_backups' not in p.parts and p.name not in ['ARTIFACT_HASHES.tsv','FINAL_STATE_VALIDATION.json']:
  hashes.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)))
pd.DataFrame(hashes).to_csv(O/'ARTIFACT_HASHES.tsv',sep='\t',index=False)
backup=O/'authority_backups';backup.mkdir(exist_ok=False)
statefiles=['PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','CODEX_HANDOFF.md']
for name in statefiles:shutil.copy2(ROOT/name,backup/name)
notice=f"""
## Phase167 24片空间端点检查完成（{now}）
冻结严格规则：上游24/24可用，接收端20/24可用（FFPE9/12、frozen11/12），576设计528可用；线性/秩、基础/扩展一致。接收端排除GSM5924034/38/39及GSM5924051，仅2/3/11/17源spot。472/472核查PASS，52输入/来源SHA不变，旧试点成员/图/资格一致；float32与float64条件数差异记录保留。
患者映射仍UNRESOLVED；GEO及原作者8文件/5issue无明确GSM→独立患者对照，c_2潜在跨制备重复。GSM5924048有两个大连通区保留。无关联效应/P/FDR/中介/新图。
报告：{O/'ENDPOINT_CHECK_CN.md'}
下一步可逐片效应描述：24上游、20接收，制备分层；患者未确认不汇总患者P/FDR。目标基因剔除敏感性未实施，旧状态不改名C4。原MTA/MTAP/SLC16A3与整体巨噬+C4定位保留。

"""
for name in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']:
 p=ROOT/name;tmp=p.with_name(p.name+'.phase167.tmp');tmp.write_text(notice+p.read_text(),encoding='utf-8');tmp.replace(p)
p=ROOT/'RESULTS_MANIFEST.tsv'
with p.open(newline='') as f: reader=csv.reader(f,delimiter='\t');allrows=list(reader)
header=allrows[0];assert len(header)==8,header
assert not any(r and str(r[0]).startswith('phase167_') for r in allrows[1:])
reg=['ENDPOINT_CHECK_CN.md','STATUS.json','SECTION_ENDPOINT_SUMMARY.tsv','DESIGN_ENDPOINT_ELIGIBILITY.tsv','GRAPH_QC.tsv','NEIGHBOR_OVERLAP.tsv','CONNECTED_COMPONENT_QC.tsv','REGION_GENE_COVERAGE.tsv','PATIENT_MAPPING_REVIEW.tsv','PATIENT_SOURCE_REVIEW.json','AUTHOR_GEO_LABEL_COMPARISON.tsv','INPUT_MANIFEST.tsv','VALIDATION.tsv','RELOAD_VALIDATION.tsv','PILOT_DESIGN_COMPARISON.tsv','FROZEN_CONFIG.json','COMPUTE_STATUS.json','EXECUTION_NOTES.json','run_check.py','validate_check_v3.py','finalize_check.py','ARTIFACT_HASHES.tsv']
for name in reg:
 path=O/name;assert path.exists()
 allrows.append(['phase167_'+name,'Phase167','technical_endpoint_check',str(path),now,'PASS','472/472 checks; 24 source, 20 receiver sections','Patient mapping unresolved; zero association tests and figures; sha256='+sha(path)])
tmp=p.with_name(p.name+'.phase167.tmp')
with tmp.open('w',newline='') as f:csv.writer(f,delimiter='\t',lineterminator='\n').writerows(allrows)
tmp.replace(p)
gate=dict(stage='Phase167',four_backups_exist=all((backup/n).exists() for n in statefiles),three_latest_notices=all((ROOT/n).read_text().startswith(notice) for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']),manifest_rows=len(reg),manifest_paths_exist=all((O/n).exists() for n in reg),status=status['status'],authority_hashes={n:sha(ROOT/n) for n in statefiles})
assert gate['four_backups_exist'] and gate['three_latest_notices'] and gate['manifest_paths_exist']
(O/'FINAL_STATE_VALIDATION.json').write_text(json.dumps(gate,ensure_ascii=False,indent=2))
print(json.dumps(status,ensure_ascii=False,indent=2));print(json.dumps(gate,ensure_ascii=False,indent=2))
