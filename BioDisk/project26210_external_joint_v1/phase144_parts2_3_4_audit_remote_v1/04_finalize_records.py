from pathlib import Path
import json,csv,hashlib,datetime,shutil,io
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=P/'phase144_parts2_3_4_audit_remote_v1'
v=json.loads((O/'REPORT_VALIDATION.json').read_text());assert v['passed']==v['checks']==1387
backup=O/'state_before';backup.mkdir(exist_ok=False)
for name in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md','RESULTS_MANIFEST.tsv']:shutil.copy2(P/name,backup/name)
block='''## 2026-09-14 Phase144 第二、三、四部分已有结果审计完成

- 用户明确要求审计目前已有的第二、第三、第四部分结果。按Y01-Y03髓系、C01-C05通信/临床、S01-S06空间核对，加入Phase109与Phase70/97关键补充，共16项状态。
- 权威审计目录：`phase144_parts2_3_4_audit_remote_v1`。报告`PARTS2_3_4_AUDIT_CN.md`；完整统计附件`ALL_PRIMARY_STATISTICS.md`；模块状态`MODULE_COVERAGE.tsv`。
- 盘点2167文件、冻结并复核46份来源SHA256；历史检验口径及round_trip读取下1387/1387字段级数值检查通过。初轮读取精度与SciPy并列秩auto差异完整保留，现代auto敏感性不改变Y03显著性定位。不是全流程重跑或大型对象重新验算。
- 第二部分：状态排序和regulon结果可用；SPP1/C1QC相对程序下降，不能称统一替代；JUNB保留研究反向。第三部分：Phase109已填补SPP1细分通信旧缺口，SPP1与内皮双向兼容性成立，上游TYMP-SPP1丰度/特异性仍不支持。第四部分：SPP1状态程序-血管重塑获FFPE及frozen分别支持，TYMP-SPP1一般化链不支持；正式S05病理仍未完成。
- 审计处理：小样本statistical_support不等于FDR显著；Cox主log-rank不能转移显著性给HR/交互。新增事后审计BH单独标注，ICAM5主要轴调整OS HR的q=0.05999。Phase97底图q05丰度与气泡组成调整表达程序是不同端点，复用必须明确图注。S04 n=0/2占位P/FDR不得作为有效阴性推断。
- 未修改原分析结果、客户湿实验、原报告或图件；不在本地写进度。无新模型/后台计算。当前断点：用户审阅三部分审计，随后按用户指定修订/出图；此前Phase143第一部分27页报告及选图状态保持。

'''
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md']:
 f=P/name;s=f.read_text();pos=s.find('\n');out=s[:pos+1]+'\n'+block+s[pos+1:];tmp=f.with_suffix('.md.phase144.tmp');tmp.write_text(out);tmp.replace(f)
f=P/'DECISIONS.md';s=f.read_text();tmp=f.with_suffix('.md.phase144.tmp');tmp.write_text(s+'\n\n'+block);tmp.replace(f)
files=['PARTS2_3_4_AUDIT_CN.md','ALL_PRIMARY_STATISTICS.md','MODULE_COVERAGE.tsv','SOURCE_SHA256_ALL.tsv','REPORT_VALIDATION.json','numerical_v2/NUMERICAL_CHECKS.tsv','numerical_v2/AUDIT_COMPUTATION_STATUS.json','numerical_v2/Y03_SCIPY_AUTO_SENSITIVITY.tsv','numerical_v2/C05_cox_AUDIT_BH.tsv','numerical_v2/C05_interaction_AUDIT_BH.tsv']
with (O/'DELIVERABLE_SHA256.tsv').open('w') as h:
 w=csv.writer(h,delimiter='\t');w.writerow(['path','sha256'])
 for rel in files:w.writerow([str(O/rel),hashlib.sha256((O/rel).read_bytes()).hexdigest()])
f=P/'RESULTS_MANIFEST.tsv';orig=f.read_text();header=next(csv.reader(io.StringIO(orig),delimiter='\t'));buf=io.StringIO();w=csv.DictWriter(buf,fieldnames=header,delimiter='\t',lineterminator='\n')
for i,rel in enumerate(files,1):
 w.writerow(dict(result_id=f'phase144_audit_{i:02d}',stage='Phase144',artifact_type=Path(rel).suffix.lstrip('.'),absolute_path=str(O/rel),generated_at=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),validation_status='PASS_RESULT_TABLE_AUDIT',validation_evidence=str(O/'REPORT_VALIDATION.json'),notes='Existing-results audit only; SHA256='+hashlib.sha256((O/rel).read_bytes()).hexdigest()))
tmp=f.with_suffix('.tsv.phase144.tmp');tmp.write_text(orig+('' if orig.endswith('\n') else '\n')+buf.getvalue());tmp.replace(f)
assert 'Phase144' in (P/'CODEX_HANDOFF.md').read_text()[:200]
assert len([r for r in csv.DictReader(f.open(),delimiter='\t') if r['stage']=='Phase144'])==len(files)
(O/'FINAL_COMPLETION.json').write_text(json.dumps({'status':'COMPLETE','four_state_files_updated':True,'manifest_artifacts':len(files),'report':str(O/'PARTS2_3_4_AUDIT_CN.md')},ensure_ascii=False,indent=2))
print((O/'FINAL_COMPLETION.json').read_text())
