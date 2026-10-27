from pathlib import Path
import json,hashlib,datetime,shutil,re
import pandas as pd,numpy as np
from scipy.stats import false_discovery_control
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');R=P/'part2_rebuild_v1';O=P/'phase149_part2_scope_initial_comparison_20260915_v1';E=O/'evidence'
checks=[]
def ck(name,v):checks.append(dict(check=name,passed=bool(v)))
manifest=pd.read_csv(O/'SOURCE_MANIFEST.tsv',sep='\t')
for _,r in manifest.iterrows():
 p=Path(r.path);ck('source_hash_'+str(p),p.stat().st_size==r.bytes and hashlib.sha256(p.read_bytes()).hexdigest()==r.sha256)
d=pd.read_csv(E/'MATCHED_UNIVERSE_2X2_SUMMARY.tsv',sep='\t');p=pd.read_csv(E/'MATCHED_UNIVERSE_2X2_PATIENT_EFFECTS.tsv',sep='\t')
ck('12_matched_summary_rows',len(d)==12);ck('same_24_patients',set(d.n_patients)=={24});ck('patient_keys_unique',not p.duplicated(['patient','network','label','regulon']).any())
sets=[set(g.patient) for _,g in p.groupby(['network','label','regulon'])];ck('same_patients_all_12_endpoints',all(s==sets[0] for s in sets))
for key,g in d.groupby(['network','label']):ck(str(key)+'_BH',np.allclose(false_discovery_control(g.P),g.FDR_posthoc_threeTF))
cov=pd.read_csv(E/'TRAJECTORY_COVERAGE.tsv',sep='\t')
ck('primary_C4_442',int(cov.loc[cov['mode']=='primary','C4_with_pseudotime'].iloc[0])==442)
ck('sensitivity_C4_1743',int(cov.loc[cov['mode']=='without_stress_cycle','C4_with_pseudotime'].iloc[0])==1743)
a=pd.read_csv(E/'NEW_C4_VS_JULY_SEVEN_SUBTYPES.tsv',sep='\t');ck('all_C4_accounted',a.new_C4_cells.sum()==2674)
ck('July_overlap_36',a.loc[a.legacy_subtype=='SPP1+Macro','new_C4_cells'].iloc[0]==36)
for fn in ['PART2_COMPLETENESS_AND_INITIAL_COMPARISON_CN.md','REVIEW_SUMMARY_CN.md']:
 t=(O/fn).read_text();ck(fn+'_Chinese_intact','第二部分' in t and '???' not in t);ck(fn+'_no_templates',not any(x in t for x in ['{matched}','{coverage}','{major}','{traj}']))
ck('no_figures',not any(p.suffix.lower() in ['.png','.pdf','.svg','.jpg'] for p in O.rglob('*')))
assert all(c['passed'] for c in checks),[c for c in checks if not c['passed']]
(O/'VALIDATION.json').write_text(json.dumps(dict(status='PASS',checks=len(checks),passed=len(checks),details=checks),ensure_ascii=False,indent=2))
B=O/'state_before_update';B.mkdir(exist_ok=False)
for q in [P/'CODEX_HANDOFF.md',P/'PROJECT_STATUS.md',P/'DECISIONS.md',P/'RESULTS_MANIFEST.tsv',R/'STATUS.json']:shutil.copy2(q,B/q.name)
note="""## 2026-09-15 Phase149 第二部分完成度与初版对照审阅完成

当前权威审阅：phase149_part2_scope_initial_comparison_20260915_v1/PART2_COMPLETENESS_AND_INITIAL_COMPARISON_CN.md；简明版REVIEW_SUMMARY_CN.md。原方案两份DOCX及提取文本哈希一致。第二部分核心方法/患者计算完成；轨迹流线图、regulon瀑布热图、最终客户稿待内容审核后制作。三TF一致级联上扬和普遍TYMP驱动终末分化不获支持，不能把此当作算法未完成。
新增核对：C4与7月SPP1+Macro交集36/2674（1.35%），与8月规则SPP1-state交集1681/2674（62.86%）；三套定义不能等同。主Monocle3有根C4仅442/2674（16.53%），去应激/增殖1743/2674（65.18%）；主模型8/8、敏感性19/19患者C4晚于单核，但不普遍晚于其他成熟巨噬。
新版JUN/JUNB/CEBPD在C4相对其他巨噬均降低；JUN主轨迹正向，JUNB方法依赖，CEBPD主轨迹无支持。PPARG状态及多方法轨迹方向更一致，TFEB为状态差异，THRA仅5靶基因。固定相同24患者的2标签×2网络事后对照确认CEBPD结论随群定义和网络改变，不能仅归因样本量。
9月15日已完成SCENIC（05:35）、后续患者统计/交接/独立验收（05:40，1364/1364 PASS），235 regulons/67949预测边/33098细胞；上游1500/1500 PASS。本次审阅未重训模型、未绘图。细胞/患者对照、全关键统计和源哈希保存在evidence及SOURCE_MANIFEST。
当前阶段：审阅完成，待用户讨论内容与图件；不自动新建图，不重跑第三四部分。后续如使用新C4通信或空间丰度须按新定义衔接。
"""
repair=[]
for fn in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 path=P/fn;text=path.read_text();parts=re.split(r'(?=^## )',text,flags=re.M);kept=[];removed=[]
 for part in parts:
  if part.startswith('## 2026-09-15') and part.count('?')>20:
   removed.append(part)
  else:kept.append(part)
 if removed:
  (B/(fn+'.garbled_automation_prefix.txt')).write_text(''.join(removed))
  repair.append(dict(file=fn,removed_garbled_sections=len(removed),backup=str(B/(fn+'.garbled_automation_prefix.txt'))))
 path.write_text(note+'\n'+''.join(kept),encoding='utf-8')
(O/'RECORD_ENCODING_REPAIR.json').write_text(json.dumps(dict(reason='Earlier PowerShell stdin used a non-UTF8 encoding for new automation prose; original scientific sections preserved; restored readable milestones from machine-readable statuses.',repairs=repair),ensure_ascii=False,indent=2))
path=R/'STATUS.json';d=json.loads(path.read_text());d['remaining_analysis'].update(status='COMPUTATION_AND_CONTENT_AUDIT_COMPLETE_PRESENTATION_PENDING',current_step='user_content_review',supervisor_pid=None,review=str(O/'PART2_COMPLETENESS_AND_INITIAL_COMPARISON_CN.md'));path.write_text(json.dumps(d,ensure_ascii=False,indent=2))
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
with (P/'RESULTS_MANIFEST.tsv').open('a') as f:
 for name in ['PART2_COMPLETENESS_AND_INITIAL_COMPARISON_CN.md','ORIGINAL_PLAN_COVERAGE.tsv','evidence/MATCHED_UNIVERSE_2X2_SUMMARY.tsv','VALIDATION.json']:
  f.write('\t'.join(['phase149_'+Path(name).stem,'phase149','audit',str(O/name),now,'REVIEW_COMPLETE','Original plan + July/August initial comparisons','No figures; no model rerun'])+'\n')
status=dict(status='COMPLETE_SCOPE_AND_INITIAL_COMPARISON_REVIEW',validation_checks=len(checks),validation_passed=len(checks),time_utc=now,core_part2_computation_complete=True,original_full_deliverable_complete=False,figures_created=False,remaining=['trajectory_figure','regulon_waterfall_figure','final_client_narrative_after_review'])
(O/'STATUS.json').write_text(json.dumps(status,indent=2))
files=['PART2_COMPLETENESS_AND_INITIAL_COMPARISON_CN.md','REVIEW_SUMMARY_CN.md']
receipt=[]
for f in files:
 path=O/f;receipt.append(dict(file=f,bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
(O/'DELIVERY_MANIFEST.json').write_text(json.dumps(receipt,indent=2))
print(json.dumps(status));print(json.dumps(receipt));print('encoding_repair',repair)

