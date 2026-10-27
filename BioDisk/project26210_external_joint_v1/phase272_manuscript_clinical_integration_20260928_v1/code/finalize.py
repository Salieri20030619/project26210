from pathlib import Path
import json,hashlib,shutil,datetime,csv,re
p=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase272_manuscript_clinical_integration_20260928_v1');r=p.parent
b=json.loads((p/'report/report_materialized.json').read_text());build=json.loads((p/'qc/BUILD.json').read_text());render=json.loads((p/'qc/RENDER_VALIDATION.json').read_text())
assert build['pages']==43 and build['tables']==34
for x in json.loads((p/'qc/INPUTS.json').read_text()):assert hashlib.sha256(Path(x['path']).read_bytes()).hexdigest()==x['sha256']
i=next(i for i,x in enumerate(b) if x.get('text','').startswith('2.13 '))
j=next(j for j in range(i+1,len(b)) if b[j]['type']=='h2')
section='\n\n'.join(('## ' if x['type']=='h2' else '')+x['text'] for x in b[i:j])
(p/'delivery/临床小节正文.md').write_text(section,encoding='utf-8')
pdf=p/'report/26210_manuscript_with_clinical_section_v1.pdf';shutil.copy2(pdf,p/'delivery'/pdf.name)
msg="""# 本次正文更新
当前版本：Phase272，2026-09-28，43页。

- 新增2.13临床小节（第19页），接在C4状态排序及图4说明之后，空间分析之前；空间结果顺延为2.14-2.21。
- 正文精选文献SPP1五基因的CheckMate OS复现、连续及临床调整模型，以及CM025 OS治疗交互；之后介绍本项目C4状态评分及时间依赖性PFS关联。
- 使用Phase271统一性别标签后的最终值：文献评分调整后OS HR1.270/P0.01611/FDR0.06444；C4时变增量P0.0007173/FDR0.001435。旧调整值未混入正文。
- 同步补充摘要、1.11临床方法、第1.19节统计口径、第3.4节讨论及结论。关键阴性和JAVELIN复现差异集中第3.8节；适用范围集中第3.9节。
- 新增附表S19/S20（第40页），参考文献补至22条。原32张表的数据及5张主图、补充图S1保持原内容。未新增模型或图件。
- 配套clinical_statistics保留完整C4模型、交互、时间模型、敏感性及文献再分析表。文献的未调整/中位数组、JAVELIN调整组和统一性别后的CheckMate调整组分别存放，避免把旧四水平性别模型当最终值。
- 全43页已渲染检查；最后版式修订仅改变1、5、39-41页并重新检查，临床附表与表注同页。正文文字、编号、页码、图件跳转和源文件校验通过。

旧Phase267完整保留，当前正式稿以本目录PDF为准。
"""
(p/'delivery/更新说明.md').write_text(msg,encoding='utf-8')
(p/'README_REPRODUCE.md').write_text("Authoritative final text: report/report_materialized.json. Set PYTHONPATH to Phase198/deps and Phase206/deps; run code/render_report.py with batch_bench Python, then code/validate_render.py. No statistical rerun required. Final renderer keeps clinical tables and notes together. Initial edit script and input hashes preserve source provenance; final JSON includes the documented wording/layout polish.\n")
val=dict(status='PASS',completed_at=datetime.datetime.now().astimezone().isoformat(),pages=43,old_tables_preserved=32,new_tables=2,main_figures=5,supplementary_figures=1,new_section='2.13',clinical_result_page=19,clinical_supplement_page=40,methods_sections=19,results_sections=21,reference_count=22,new_analyses=0,all_paragraphs_verified=render['paragraphs_verified'],links_checked=render['links'],source_hashes_unchanged=True,visual_review_all_pages=True,final_changed_pages_rechecked=[1,5,39,40,41],pdf_sha256=build['sha256'])
(p/'qc/FINAL_VALIDATION.json').write_text(json.dumps(val,ensure_ascii=False,indent=2))
shutil.copy2(p/'qc/FINAL_VALIDATION.json',p/'delivery/FINAL_VALIDATION.json')
manifest=[]
for f in sorted((p/'delivery').rglob('*')):
 if f.is_file():manifest.append(dict(file=str(f.relative_to(p/'delivery')),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
(p/'delivery/DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
note=f"""## Phase272（2026-09-28）：临床小节整合正式稿完成
用户要求在分群后、空间前补充临床小节，并精选此前回顾性复现。已以Phase267为底稿生成43页正式PDF，新增结果2.13（第19页），位于原2.12状态排序/图4说明之后、原空间2.13之前；空间顺延2.14-2.21。精选SPP1文献五基因CheckMate OS中位数/连续/临床调整结果与CM025 OS交互作为参照，再展开本项目C4十九基因状态评分和时变PFS关联，Sunitinib条件关联保留。
同步新增方法1.11，原空间方法顺延1.12-1.18，统计1.19；摘要纳入311/726临床队列和关键关联；讨论3.4加入临床异质性，3.8集中OS/交互阴性及JAVELIN精确复现差异，3.9说明回顾性/时间风险集/空间临床不同患者等边界，结论同步。附表S19/S20第40页，参考文献19-22补全。未新增分析或新图；原32表数值完全不变，5主图+S1原矢量源SHA相同。
统计使用Phase271统一性别标签的最终值（Jiang5 adjusted OS HR1.269884/P0.01611/q0.06444；C4时变P0.0007173/q0.001435），没有沿用旧四水平性别结果。完整文献及C4统计按版本拆分随delivery/clinical_statistics交付。
权威目录phase272_manuscript_clinical_integration_20260928_v1；report/report_materialized.json为最终正文源；PDF 26210_manuscript_with_clinical_section_v1.pdf。43页视觉检查及最后变更页1/5/39-41复查通过，正文/表格文字、编号、页码、{render['links']}图件链接通过，PDF SHA256 {build['sha256']}。旧Phase267留档，当前稿以Phase272为准。本地交付目录Downloads/26210/Phase272_临床小节整合_v1。无后台计算。

"""
bk=p/'qc/root_record_backup';bk.mkdir(exist_ok=True)
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 f=r/name;shutil.copy2(f,bk/name);f.write_text(note+f.read_text(),encoding='utf-8')
f=r/'RESULTS_MANIFEST.tsv';shutil.copy2(f,bk/f.name)
with f.open('a') as out:
 wr=csv.writer(out,delimiter='\t',lineterminator='\n')
 for item in manifest:wr.writerow(['phase272_'+item['file'],'phase272','clinical_manuscript_integration',str(p/'delivery'/item['file']),val['completed_at'],'PASS',str(p/'qc/FINAL_VALIDATION.json'),'New section2.13 before spatial; selected retrospective references plus C4 state; original figures/statistics preserved'])
print(json.dumps(val,ensure_ascii=False));print('DELIVERY_FILES',len(manifest))
