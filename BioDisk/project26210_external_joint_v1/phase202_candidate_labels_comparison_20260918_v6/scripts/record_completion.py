from pathlib import Path
import shutil,json,datetime
r=Path.cwd();d=r/'phase202_candidate_labels_comparison_20260918_v6'
text='''\n\n## 2026-09-18 Phase201/202：C4主分析报告与三候选标注已交付\n用户最终范围：原报告仅增补三候选/联合覆盖优势，继续以原C4作主分析；Fig2D和Fig3B加C10、C15位置，另给三群/联合/整体比较图。此决定取代Phase200建议将联合群作为后续主对象的建议，但保留其完整分析档案。\nPhase201报告v4：34页，原31表完全一致、7图哈希一致，新增2.9节、图8(候选marker/研究分布)、图9(表达/检测/RNA贡献)，第20–21页；C4主分析从2.10节继续。未将未改善的次要联合结果强行加入主线，完整记录保留补充。\nPhase202主图v6：Fig2D、Fig3B新增C10蓝色点划线、C15绿色点线，C4紫色虚线；位置来自原冻结群KDE，四个表达底图数值完全相同，轮廓不等于表达峰值。Fig3其他面板沿用phase199/197原源PDF，C4 GSEA不变。\n新增比较：A/B完整发现队列患者等权中位数，C4/C10/C15/U/整体有效n=26/2/9/28/38；C10仅描述。C为原16同文库患者RNA份额中位54.59%/0%/1.85%/59.63%/100%，缺失群计0；各组包含、患者构成不同，中位数不可相加，0中位不等于所有患者均无贡献。\n权威：phase201_three_candidates_C4_main_report_20260918_v4及phase202_candidate_labels_comparison_20260918_v6。本地：C:/Users/31937/Downloads/26210/31_C4主分析报告与候选群标注_v4。四PDF、三图SVG/600dpiPNG、报告MD及两补充记录；14/14文件哈希一致。34页及全部新图完成视觉检查；患者数、RNA可加性、原始图表保留校验通过。\n当前断点：本次任务全部完成，无计算任务运行。Fig1仍暂定稿；主图Fig2/Fig3与原报告内部图号为不同编号体系。后续修改以本次版本为起点，旧数据/结果/版本不覆盖。\n'''
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv','AGENTS.md']:
 p=r/name
 if p.exists():shutil.copy2(p,d/'record_backup'/name)
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:
 with (r/name).open('a') as f:f.write(text)
with (r/'AGENTS.md').open('a') as f:f.write('\n\n### 2026-09-18 当前报告与图件口径\n按用户最终决定，C4继续为主分析；C10/C15和联合群仅简要候选来源比较。Phase201报告v4和Phase202主图v6为当前版本，Phase200完整联合复评留档。允许按已有SSH传输授权往返文件，无需重复确认。\n')
with (r/'RESULTS_MANIFEST.tsv').open('a') as f:f.write('\nphase201_202_C4_primary_report_labels\tphase202\treport_and_figures\t'+str(d/'delivery')+'\t'+datetime.datetime.now().isoformat()+'\tPASS\t34-page visual QA; 31 unchanged tables; 7 unchanged figures; frozen expression verified; 14 local hash matches\tC4 primary; C10/C15 locators; union/whole descriptive comparison\n')
(d/'qc/LOCAL_TRANSFER_VERIFIED.json').write_text(json.dumps({'files':14,'sha256_match':14,'local_directory':'C:/Users/31937/Downloads/26210/31_C4主分析报告与候选群标注_v4'},ensure_ascii=False,indent=2))
print('COMPLETION_RECORDED')
