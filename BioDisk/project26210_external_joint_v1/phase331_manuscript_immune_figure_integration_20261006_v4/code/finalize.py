from pathlib import Path
import json,hashlib,shutil,re
P=Path(__file__).resolve().parents[1];R=P.parent;O=R/'phase327_abstract_clarity_revision_20261006_v1';F=R/'phase330_figure_F_tissue_immune_20261006_v3';sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
b=json.loads((P/'report/report_materialized.json').read_text());old=json.loads((O/'report/report_materialized.json').read_text())
for nb in b:
 if nb['type']=='table':
  ob=old[nb['_phase327_index']];assert all(nb[k]==ob[k] for k in ['headers','rows','caption'])
assert [x['text'] for x in b if x['type'] in ['h1','h2']]==[x['text'] for x in old if x['type'] in ['h1','h2']]
assert b[9]['text']==old[9]['text']
assert sha(P/'figures/Figure5.pdf')==sha(F/'delivery/26210_immune_spatial_clinical_A4_v3.pdf')
for f in (O/'figures').iterdir():assert sha(f)==sha(P/'figures'/('Figure6.pdf' if f.name=='Figure5.pdf' else f.name))
# Preserve final renderer choices in the preparation script as well.
ip=P/'code/integrate.py';s=ip.read_text();needle="(P/'code/render_report.py').write_text(s)"
add='''s=s.replace("else 125))", "else 80 if x['text'].startswith('2.26 ') else 125))")
s=s.replace("  tt=Table(data,", "  assert len(weights)==n and all(v>0 for v in weights), (x['caption'],n,weights)\\n  tt=Table(data,")
'''
if add not in s:ip.write_text(s.replace(needle,add+needle))
compile(ip.read_text(),str(ip),'exec')
review=json.loads((P/'qc/PAGE_REVIEW_MAP.json').read_text());build=json.loads((P/'qc/BUILD.json').read_text());assert build['pages']==54 and build['figures']==8 and build['tables']==41
q={'status':'PASS','pages':54,'tables_unchanged':41,'all_section_headings_unchanged':True,'abstract_revision_preserved':True,'figure5_page':25,'figure5_caption_page':26,'figure6_page':29,'new_figure_scale':1.0,'all_old_figures_hash_identical':True,'new_figure_hash_identical':True,'internal_links':78,'changed_pages_visually_reviewed':review['changed_pages'],'other_pages_pixel_identical_excluding_footer':len(review['same_content_and_pixels_excluding_footer']),'final_review_images':[f'review_{i:02d}.png' for i in range(1,8)],'resolved_draft_layout_issues':['Figure-number branch replacement was limited to figure flowable; table widths restored and asserted.','Short final Results subsection kept on preceding page; unnecessary near-empty page removed.'],'logical_checks':['New figure follows clinical Results 2.17 and introduces spatial Results 2.18 onward.','Figure 5 panel references match immune, bulk, spatial and clinical measurements.','Former Figure 5 is Figure 6 throughout methods/results/caption and bookmarks.','Clinical treatment-interaction and experimental evidence retained in text/supplement despite absent composite panels.','No new numerical analyses or revised evidence claims.']}
name='26210_manuscript_immune_figure_integrated_v4.pdf';shutil.copy2(P/'report'/name,P/'delivery'/name);q['pdf_sha256']=sha(P/'delivery'/name);q['pdf_bytes']=(P/'delivery'/name).stat().st_size
(P/'qc/COMPLETE.json').write_text(json.dumps(q,ensure_ascii=False,indent=2));(P/'delivery/整合说明.txt').write_text('综合图v3已作为图5加入完整文章第25页，图注第26页。原空间图顺延为图6（第29页），对应正文引用、书签和链接同步更新。\n全文54页，章节编号和41张统计表保持；旧稿保留。已检查新增与重排页面，另41页除页码外与原稿逐像素一致。\nPDF SHA256: '+q['pdf_sha256']+'\n',encoding='utf-8-sig')
manifest=[{'name':f.name,'bytes':f.stat().st_size,'sha256':sha(f)} for f in sorted((P/'delivery').iterdir())];(P/'delivery/文件校验.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
e=f'''\n## Phase331（2026-10-06）：用户接受综合图v3，已正式插入完整文章v4\n用户“暂时先这样吧，加入到完整文章中”，已将Phase330最终8面板综合图作为图5置于2.17临床结果之后、2.18空间结果之前。图5完整A4原尺寸嵌入第25页，图注第26页；原图5顺延为图6（第29页），所有正文/方法/图注/附录图件概览及图内跳转书签同步。章节编号不变，41张表及全部原统计值、旧7图文件哈希均保持。\n当前正式全文：{P/'delivery'/name}，54页，{q['pdf_bytes']}字节，SHA256 {q['pdf_sha256']}。继承Phase327摘要修订；Phase327和所有历史图稿保留。新图源SHA与Phase330一致，scale=1，8图共78个内部链接核对通过。13个改变页面逐页视觉检查，另外41页除页码外与已审Phase327逐像素一致；图表宽度及短节分页问题在本轮排版检查中修复。验收qc/COMPLETE.json，原始映射/输入/渲染脚本均在本阶段目录。\n交付Downloads/26210/Phase331_综合图并入全文_v4。当前断点：综合图已经采用并入全文，后续以本v4为正式入口；不再使用“独立图未插入”的旧状态。无新分析或统计计算。\n'''
for fn in ['CODEX_HANDOFF.md','PROJECT_STATUS.md']:
 f=R/fn;f.write_text(e+'\n'+f.read_text(),encoding='utf-8')
with (R/'DECISIONS.md').open('a') as f:f.write(e)
with (R/'RESULTS_MANIFEST.tsv').open('a') as f:f.write('\t'.join(['phase331_full_manuscript_v4','Phase331','full_manuscript',str(P/'delivery'/name),'2026-10-06','PASS',str(P/'qc/COMPLETE.json'),'54 pages; Figure5=accepted immune composite; Figure6=former spatial figure; SHA256='+q['pdf_sha256']])+'\n')
print('FINAL',json.dumps({'pages':54,'sha256':q['pdf_sha256'],'bytes':q['pdf_bytes']}))
