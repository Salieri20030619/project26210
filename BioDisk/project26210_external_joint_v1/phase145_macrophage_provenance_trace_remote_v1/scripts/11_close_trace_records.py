from pathlib import Path
import json,csv,hashlib,datetime,shutil
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=P/'phase145_macrophage_provenance_trace_remote_v1'
v=json.loads((O/'TRACE_VALIDATION.json').read_text());assert v['checks']==v['passed']==23
z=json.loads((O/'identity/CURRENT_RULE_STATE_PROPORTION_STATUS.json').read_text());assert z['n_patients']==40 and z['n_macrophages']==42546 and z['n_spp1_state']==8892
f=O/'MACROPHAGE_SPP1_PROVENANCE_REPORT_CN.md';s=f.read_text().replace('按barcode回填宏细胞标签','按barcode回填巨噬细胞标签')
extra='''\n### 患者内阈值对丰度解释的约束\n\n现行reference-v2共有42,546个巨噬细胞、40名患者，其中8,892个SPP1-state细胞。各患者SPP1-state比例最小18.61%、中位21.37%、最大24.32%。这个范围直接受患者内上四分位及互斥规则影响，不是完全由自然聚类频率决定。它适用于比较高程序状态的条件特征，但不能独立证明每个患者都自然存在相近比例的SPP1亚群；解释患者间丰度关联时应明确这一约束。该事实不改动现有统计值。详见identity/CURRENT_RULE_STATE_PATIENT_PROPORTIONS.tsv。\n\n'''
s=s.replace('## 六 与当前重建主线逐细胞核对',extra+'## 六 与当前重建主线逐细胞核对');f.write_text(s)
v['proportion_constraint_checked']=True;v['proportion_summary']=z;(O/'TRACE_VALIDATION.json').write_text(json.dumps(v,ensure_ascii=False,indent=2))
# Preserve compact timestamp/call evidence without duplicating thousands of feature names.
ls=(O/'legacy_object/COMMANDS.txt').read_text().splitlines();out=[]
for i,l in enumerate(ls):
 if l.endswith('$time') or l.endswith('$call'):
  j=i+1
  while j<len(ls) and ls[j].strip():j+=1
  out.extend(ls[i:j]+[''])
(O/'legacy_object/COMMAND_TIMELINE_EXCERPT.txt').write_text('\n'.join(out))
block='''\n## 2026-09-14 Phase145 巨噬细胞与SPP1详细追溯完成\n\n- 用户说明此前这部分不是其本人完成，要求先记录当前情况再详细追溯。启动快照CURRENT_STATE_BEFORE_TRACE_CN.md已在执行前保存；追溯目录phase145_macrophage_provenance_trace_remote_v1。正式报告MACROPHAGE_SPP1_PROVENANCE_REPORT_CN.md。\n- 更正：旧对象确实做过聚类。Seurat commands记录2026-07-18 RNA pca前30维邻居图及0.05/0.1/0.15/0.2聚类；最终41,105细胞、七亚型，全部与res.0.1一致。41,105/41,105匹配旧275,150细胞注释对象中的Macrophage；旧45,309巨噬中4,204未进入最终对象，删除理由/完整marker表/生成脚本尚未恢复，不能推定作者。\n- 旧SPP1群11,290细胞，10,604(93.92%)来自GSE242299；多数其他群也高度研究集中。旧七群不是已经确认的跨队列稳定细分。\n- 当前75文库重建保留旧对象38,692细胞；其余2,413均有当前QC排除记录。旧SPP1群保留10,276，其中9,232现为Macrophage、1,044为其他注释；现行规则下仅1,889对应SPP1-state，其余1,629 C1QC、5,714 Other。不是旧标签准确率，而是定义不等同。\n- 8月5日历史会话确认因空间参考敏感性而改用新版QC对象及患者内SPP1/C1QC程序，不继承旧硬标签。该规则后续用于Phase21/33/34/109。当前40患者SPP1比例18.61%-24.32%受阈值约束，不能当自然聚类丰度发现。\n- 扫描1,984脚本/文本、17个旧项目会话和3个8月5日会话；关键来源SHA256已冻结，23/23逐细胞/数量独立核对通过，原RDS原位hash完成。没有重新聚类或改写旧结果，没有联系他人，本地没有新增进度记录。\n- 当前结论：有历史分群基础，但当前数据的无监督巨噬精细分群→完整marker发现SPP1→下游验证链仍未完备。旧对象可作历史/敏感性候选，不能直接拼接当前程序状态下游。当前断点：用户审阅详细追溯；后续精细分群须按新明确任务启动。\n'''
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md']:
 f=P/name;t=f.read_text();i=t.find('\n');tmp=f.with_suffix('.md.phase145.tmp');tmp.write_text(t[:i+1]+block+t[i+1:]);tmp.replace(f)
with (P/'DECISIONS.md').open('a') as h:h.write(block)
files=['CURRENT_STATE_BEFORE_TRACE_CN.md','MACROPHAGE_SPP1_PROVENANCE_REPORT_CN.md','TRACE_VALIDATION.json','SOURCE_FILE_MANIFEST.tsv','legacy_object/OBJECT_SUMMARY.txt','legacy_object/COMMAND_TIMELINE_EXCERPT.txt','identity/IDENTITY_STATUS.json','identity/LEGACY_PARENT_STATUS.json','identity/LEGACY_REUSE_SUMMARY.tsv','identity/legacy_subtype_by_study.tsv','identity/legacy_subtype_by_current_reference_state.tsv','identity/CURRENT_RULE_STATE_PATIENT_PROPORTIONS.tsv','SOURCE_SEARCH_FILES.json','project_session_trace/MATCH_INDEX.json','reference_transition_trace/MATCH_INDEX.json']
with (O/'DELIVERABLE_SHA256.tsv').open('w') as h:
 w=csv.writer(h,delimiter='\t');w.writerow(['file','sha256'])
 for name in files:w.writerow([name,hashlib.sha256((O/name).read_bytes()).hexdigest()])
f=P/'RESULTS_MANIFEST.tsv';header=next(csv.reader(f.open(),delimiter='\t'))
with f.open('a') as h:
 w=csv.DictWriter(h,fieldnames=header,delimiter='\t',lineterminator='\n')
 for i,name in enumerate(files,1):
  w.writerow(dict(result_id=f'phase145_trace_{i:02d}',stage='Phase145',artifact_type=Path(name).suffix.lstrip('.'),absolute_path=str(O/name),generated_at=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),validation_status='PASS_TRACE_WITH_EXPLICIT_PROVENANCE_GAPS',validation_evidence=str(O/'TRACE_VALIDATION.json'),notes='Read-only historical provenance; no new clustering; SHA256='+hashlib.sha256((O/name).read_bytes()).hexdigest()))
assert len([r for r in csv.DictReader(f.open(),delimiter='\t') if r['stage']=='Phase145'])==len(files)
for r in csv.DictReader((O/'DELIVERABLE_SHA256.tsv').open(),delimiter='\t'):assert hashlib.sha256((O/r['file']).read_bytes()).hexdigest()==r['sha256']
status={'status':'TRACE_COMPLETE_WITH_DOCUMENTED_MISSING_ORIGINAL_MATERIALS','four_state_files_updated':True,'artifacts_registered':len(files),'source_files_hashed':len(list(csv.DictReader((O/'SOURCE_FILE_MANIFEST.tsv').open(),delimiter='\t'))),'checks_passed':23,'new_clustering_run':False,'background_jobs_remaining':False}
(O/'FINAL_COMPLETION.json').write_text(json.dumps(status,ensure_ascii=False,indent=2));print(json.dumps(status,ensure_ascii=False))
