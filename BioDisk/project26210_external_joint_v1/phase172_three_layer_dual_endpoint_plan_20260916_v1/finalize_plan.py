from pathlib import Path
import hashlib,json,csv,shutil,datetime
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=R/'phase172_three_layer_dual_endpoint_plan_20260916_v1'
p=O/'DETAILED_ANALYSIS_PLAN_CN.md'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
expected='36ac081892f14330085561210a1fda20dcd14b36579be893d57a2e4768c1b119'
assert sha(p)==expected,'Transfer hash mismatch'
s=p.read_text(encoding='utf-8')
required=['## 三、第一层','## 四、第二层','## 五、第三层A','## 六、第三层B','## 七、两条分支','## 八、统一统计','## 九、湿实验','## 十、阶段顺序','## 十一、最终交付','## 十二、服务器来源']
assert all(x in s for x in required)
assert '用户确认目前提供的湿实验依据只有既有代谢组图' in s
assert '第一部分报告已完成' in s
assert '免疫功能和空间结构是不同终点' in s
assert '新增分析未启动' in s
sources=[
'phase162_MTA_polyamine_part1_reassessment_20260915_v1/MTA_POLYAMINE_PART1_REASSESSMENT_CN.md',
'phase163_polyamine_MTA_serial_receiver_20260915_v1/SERIAL_MTA_AND_RECEIVER_REVIEW_CN.md',
'phase149_part2_scope_initial_comparison_20260915_v1/PART2_COMPLETENESS_AND_INITIAL_COMPARISON_CN.md',
'phase151_TYMP_SPP1_threeTF_role_20260915_v1/TYMP_SPP1_THREE_TF_ROLE_REVIEW_CN.md',
'phase144_parts2_3_4_audit_remote_v1/PARTS2_3_4_AUDIT_CN.md',
'phase168_reference_source_audit_20260916_v1/SPATIAL_26_SECTION_PATIENT_MAP.tsv',
'phase169_GSE175540_treatment_audit_20260916_v1/TREATMENT_AUDIT_CN.md',
'phase170_MTA_spatial_three_part_analysis_20260916_v1/THREE_PART_SPATIAL_RESULTS_CN.md']
assert all((R/n).exists() for n in sources)
now=datetime.datetime.now().astimezone().isoformat()
local='C:/Users/31937/Downloads/26210/phase172_detailed_plan_20260916_v1/TYMP_MTA_SPP1_后续三层与免疫空间双分支详细计划_v1.md'
manifest=dict(stage='Phase172',created_at=now,server_plan=str(p),local_delivery=local,sha256=sha(p),bytes=p.stat().st_size,local_server_identical=True,source_report=dict(path='D:/电脑管家迁移文件/xwechat_files/wxid_913cjq0lvmfd32_3e57/msg/file/2026-09/26210_TYMP_report_selected_Fig1_to_Fig7_draft_v6.pdf',sha256='204d1cf400e360746755c8473f3cf131db6de37bf41f155a628612fbe8f77f85'),sources=[dict(path=str(R/n),sha256=sha(R/n)) for n in sources])
(O/'DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
status=dict(stage='Phase172',status='PLAN_COMPLETE_DELIVERED',new_analyses_started=False,new_figures=0,validation='required sections present; eight source paths verified; local-server SHA256 identical',completed_at=now)
(O/'STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
B=O/'authority_backups';B.mkdir(exist_ok=False)
names=['PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','CODEX_HANDOFF.md','AGENTS.md']
for n in names:shutil.copy2(R/n,B/n)
notice=f"""
## Phase172 详细计划已交付：承接第一部分，免疫/空间双分支（{now}）
用户明确：第一部分报告已出，按26210_TYMP_report_selected_Fig1_to_Fig7_draft_v6.pdf的逻辑往后衔接；湿实验承担后续验证，不作为生信叙事起点。用户确认当前仅提供既有代谢图，不再等待其他湿实验说明。此条纠正Phase171中的起点和待确认状态。
计划三层：TYMP代谢背景引出MTA→MTA衔接整体单核/巨噬及SPP1状态→第三层分A免疫功能和B空间组织，分别端点、统计和交付，彼此不替代。固定C4，不强制SRM/SMS、三TF或SLC16A3中介；复用已有证据，保留关键阴性。
服务器：{p}
本地：{local}
SHA256：{expected}
本次只生成详细计划及双份交付，未启动新增分析或绘图。下次按计划从L1已有证据衔接与L2直接作用来源对照开始，不重做第一部分。
"""
for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']:
 path=R/n;tmp=path.with_name(n+'.phase172.tmp');tmp.write_text(notice+path.read_text(),encoding='utf-8');tmp.replace(path)
agents=R/'AGENTS.md'
rules="""

## 2026-09-16 用户确认后续主线与两类独立终点
以已完成第一部分26210_TYMP_report_selected_Fig1_to_Fig7_draft_v6.pdf为固定起点，不重新搭建第一部分。研究叙事由TYMP相关代谢发现引出MTA，再衔接SPP1富集巨噬及下游；客户湿实验承担后续验证，不作为生信发现起点。
用户确认当前只有既有代谢组图，不继续等待额外湿实验说明。免疫功能与空间组织是两个独立终点，必须分支分析、分别固定检验家族和报告，不能互相代替。
详细计划：phase172_three_layer_dual_endpoint_plan_20260916_v1/DETAILED_ANALYSIS_PLAN_CN.md。本次仅计划与交付，新增分析未启动。先审核文字统计，暂不绘图；不得以完整强因果中介作为所有生信结果可用的硬门槛，也不得凭叙事填补缺失箭头。
"""
tmp=agents.with_name('AGENTS.md.phase172.tmp');tmp.write_text(agents.read_text()+rules,encoding='utf-8');tmp.replace(agents)
path=R/'RESULTS_MANIFEST.tsv'
with path.open(newline='') as f:rows=list(csv.reader(f,delimiter='\t'))
assert len(rows[0])==8 and not any(r and r[0].startswith('phase172_') for r in rows[1:])
for n in ['DETAILED_ANALYSIS_PLAN_CN.md','DELIVERY_MANIFEST.json','STATUS.json','finalize_plan.py']:
 rows.append(['phase172_'+n,'Phase172','detailed_analysis_plan',str(O/n),now,'PASS','Plan only; dual endpoints; local-server delivery verified','sha256='+sha(O/n)])
tmp=path.with_name(path.name+'.phase172.tmp')
with tmp.open('w',newline='') as f:csv.writer(f,delimiter='\t',lineterminator='\n').writerows(rows)
tmp.replace(path)
gate=dict(stage='Phase172',plan_sha256=sha(p),local_server_match=True,sources_verified=8,plan_sections_verified=True,four_state_files_updated=all((R/n).read_text().startswith(notice) for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']),operating_rules_updated=True,backups=len(names))
(O/'FINAL_STATE_VALIDATION.json').write_text(json.dumps(gate,indent=2))
print(json.dumps(manifest,ensure_ascii=False,indent=2));print(json.dumps(gate))
