from pathlib import Path
import json,hashlib,shutil,datetime,csv
import pymupdf as fitz
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');p=r/'phase273_clinical_suppfig_wording_20260928_v1'
b=json.loads((p/'report/report_materialized.json').read_text());plot=json.loads((p/'qc/INDEPENDENT_PLOT_CHECK.json').read_text());render=json.loads((p/'qc/RENDER_VALIDATION.json').read_text());build=json.loads((p/'qc/BUILD.json').read_text())
assert render['pages']==45 and render['tables']==34 and len(render['figures'])==7
for x in json.loads((p/'qc/INPUTS.json').read_text()):assert hashlib.sha256(Path(x['path']).read_bytes()).hexdigest()==x['sha256']
assert not any('再分析' in x.get('text','') or '再分析' in x.get('caption','') for x in b)
for source,dest in [
(p/'report/26210_manuscript_clinical_suppfig_v1.pdf',p/'delivery/26210_manuscript_clinical_suppfig_v1.pdf'),
(p/'figures/Supplementary_Figure_S2.pdf',p/'delivery/Supplementary_Figure_S2.pdf'),
(p/'figures/Supplementary_Figure_S2.png',p/'delivery/Supplementary_Figure_S2.png'),
(p/'report/Supplementary_Figure_S2_caption_CN.md',p/'delivery/Supplementary_Figure_S2_caption_CN.md')]:shutil.copy2(source,dest)
shutil.copytree(p/'data',p/'delivery/figure_data',dirs_exist_ok=True)
i=next(i for i,x in enumerate(b) if x.get('text','').startswith('2.13 '));j=next(j for j in range(i+1,len(b)) if b[j]['type']=='h2')
(p/'delivery/临床小节正文.md').write_text('\n\n'.join(('## ' if x['type']=='h2' else '')+x['text'] for x in b[i:j]),encoding='utf-8')
note="""# 文字与临床补充图更新
当前版本：Phase273，2026-09-28，45页。

正文、摘要、讨论及附表标题中的“再分析”已按语境改为“临床队列分析”“临床转录组分析”或“临床关联”。方法与讨论保留已发表公开队列和回顾性资料说明。

新增补充图S2，单页A4，包含四个面板：
- A：CheckMate nivolumab组SPP1巨噬评分高低组OS曲线，含95%CI、删失标记和在险人数。
- B：文献五基因评分的治疗交互；同时显示CM025 OS/PFS和JAVELIN PFS，原始P/FDR均保留。
- C：C4完整调整模型的时间HR曲线，显示0-24个月的点态95%CI和在险人数。
- D：固定12个月HR的完整模型、试验分层及CM025子集比较，右侧P/FDR为完整时间效应的联合检验。

图已并入报告第44页，中文图注第45页；正文2.13仍在第19页，已增加对应面板跳转。独立交付矢量PDF和300 dpi PNG，便于后续排版。

34张表的数值、既有五张主图和补充图S1均保持不变。本轮无新假设检验或截点搜索；绘图所需曲线按原模型重建，独立核对KM估计、置信区间、风险人数及C4冻结时点效应。最终8个变更页已视觉检查，其余37页与上一版逐像素一致。

完整既有统计位于clinical_statistics；新增绘图数据位于figure_data。旧Phase272保留，当前报告以本目录PDF为准。
"""
(p/'delivery/更新说明.md').write_text(note,encoding='utf-8')
val={'status':'PASS','completed_at':datetime.datetime.now().astimezone().isoformat(),'pages':45,'tables':34,'main_figures':5,'supplementary_figures':2,'new_figure':'S2','new_panels':4,'figure_page':44,'legend_page':45,'clinical_section_page':19,'reanalysis_word_absent_in_manuscript':True,'retrospective_public_source_disclosed':True,'original_table_values_unchanged':True,'source_hashes_unchanged':True,'new_hypothesis_tests':0,'KM_independent_error':plot['max_KM_error'],'C4_curve_verified':True,'visual_review_changed_pages':plot['changed_pages'],'unchanged_pages_pixel_identical':37,'internal_links_verified':render['links'],'pdf_sha256':build['sha256']}
(p/'qc/FINAL_VALIDATION.json').write_text(json.dumps(val,ensure_ascii=False,indent=2));shutil.copy2(p/'qc/FINAL_VALIDATION.json',p/'delivery/FINAL_VALIDATION.json')
manifest=[]
for f in sorted((p/'delivery').rglob('*')):
 if f.is_file():manifest.append({'file':str(f.relative_to(p/'delivery')),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
(p/'delivery/DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
(p/'README_REPRODUCE.md').write_text('Use frozen source tables. code/export_plot_data.R reconstructs and verifies plot inputs; code/draw_s2.py renders the A4 figure; code/edit_report.py applies wording and appends S2; code/render_report.py, code/validate_render.py and code/check_plot.py rebuild and verify the report. Use batch_bench Python with Phase198/deps and Phase206/deps in PYTHONPATH. Source and previous outputs are preserved.\n')
record=f"""## Phase273（2026-09-28）：临床措辞修订与补充图S2完成
用户要求不再使用“再分析”表述，并将数张临床小图放入补充图。已在Phase272基础上将8处相关文字按语境改成临床队列分析/临床转录组分析/临床关联，正文PDF不再出现该词；方法与讨论保留公开已发表资料和回顾性来源，不改为新增独立临床队列。
新增A4补充图S2四面板：A文献五基因nivolumab OS KM（n181，高84/低97，123事件，原311队列中位数；CI/删失/风险人数，logrankP0.01484/q0.05934）；B文献五基因治疗交互CM025 OS/PFS及JAVELIN PFS（含阴性比较及原P/FDR）；C固定19基因C4完整调整的0-24月时变HR/95%CI及风险人数；D12月HR完整模型/试验分层/CM025子集，标注完整2df效应P/FDR，不作为独立验证。原模型/基因/阈值不变，无新假设检验。
图并入第44页、中文图注第45页，2.13临床小节仍第19页并补面板链接。报告45页、34表、5主图+S1/S2；旧34表数值及原图源不变。独立KM乘积限/Greenwood-log CI核对误差{plot['max_KM_error']:.3g}；C4曲线公式及冻结时点HR/CI、风险人数通过。8个变更页[1,2,19,24,27,40,44,45]视觉检查完成，其余37页逐像素与Phase272相同；{render['links']}个内部链接及所有文字/页码核对通过。
权威目录phase273_clinical_suppfig_wording_20260928_v1；正式PDF report/26210_manuscript_clinical_suppfig_v1.pdf，SHA256 {build['sha256']}。delivery含报告、S2矢量PDF/300dpi PNG、图注、正文摘录、全部既有统计和绘图数据。旧Phase272留档，当前稿以Phase273为准。本地交付目录Downloads/26210/Phase273_临床补充图与正文修订_v1。无后台任务。

"""
bk=p/'qc/root_record_backup';bk.mkdir(exist_ok=True)
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 f=r/name;shutil.copy2(f,bk/name);f.write_text(record+f.read_text(),encoding='utf-8')
f=r/'RESULTS_MANIFEST.tsv';shutil.copy2(f,bk/f.name)
with f.open('a') as out:
 wr=csv.writer(out,delimiter='\t',lineterminator='\n')
 for item in manifest:wr.writerow(['phase273_'+item['file'],'phase273','clinical_supplementary_figure_and_wording',str(p/'delivery'/item['file']),val['completed_at'],'PASS',str(p/'qc/FINAL_VALIDATION.json'),'Four-panelS2; wording revised; public retrospective origin retained; source statistics unchanged'])
print(json.dumps(val,ensure_ascii=False));print('FILES',len(manifest))
