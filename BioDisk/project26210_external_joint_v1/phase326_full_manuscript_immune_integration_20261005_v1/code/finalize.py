from pathlib import Path
import json,csv,hashlib,shutil
from datetime import datetime,timezone,timedelta
P=Path(__file__).resolve().parents[1];R=P.parent
now=datetime.now(timezone(timedelta(hours=8))).isoformat(timespec='seconds')
v=json.loads((P/'qc/RENDER_VALIDATION.json').read_text())
assert v['status']=='STRUCTURAL_PASS_VISUAL_PENDING' and v['pages']==52
changed=json.loads((P/'qc/CHANGED_FINAL_PAGES.json').read_text())
assert changed==[2,3,8,9,22,32,46],changed
b=json.loads((P/'report/report_materialized.json').read_text())
base=json.loads((R/'phase273_clinical_suppfig_wording_20260928_v1/report/report_materialized.json').read_text())
imm=json.loads((R/'phase325_immune_results_mechanism_revision_20261005_v1/inputs/DOCUMENT_CONTENT.json').read_text())
tb=[x for x in b if x['type']=='table']
assert len(tb)==41
for x in [x for x in base+imm if x['type']=='table']:
 assert any(x['headers']==y['headers'] and x['rows']==y['rows'] for y in tb)
pdf=P/'report/26210_manuscript_immune_integrated_v2.pdf';sha=hashlib.sha256(pdf.read_bytes()).hexdigest()
review={'status':'PASS','completed_at':now,'full_manuscript_read':True,'pages_reviewed':list(range(1,53)),
'final_changed_pages_rechecked':changed,'unchanged_pages_byte_identical_to_reviewed_render':45,
'logic_checks':[
'Immune methods remain detailed and are distributed by analytic type in Methods; existing clinical method retained and duplicate removed.',
'Results2.13–2.16 follow C4 and macrophage state results; existing clinical result is2.17; spatial results continue2.18–2.26.',
'C4 between-state contrast is explicitly distinguished from within-patient TYMP continuous associations and the independent11-patient immune cohort.',
'Tumour TYMP–macrophage SPP1 cross-compartment null is preserved separately from myeloid intracellular positive association.',
'SPP1 expression-program/ECM association is distinguished from old abundance-neighborhood null and lack of stable SPP1-linked CD8 interface effect.',
'Clinical effects merged once into2.17, with different endpoints/scores and negative interactions preserved; supplementary S2 retained.',
'Abstract, question, working model and conclusion reflect immune findings without upgrading association to causality.',
'Original S1–S20 statistical appendix and all source figures retained;35 references and13 main tables use consistent numbering.'
],
'visual_checks':'All52 pages reviewed; final7 changed pages reviewed again; no clipped text, missing glyphs, orphan heading, numeric wrapping or broken table headers. Full-page native figures retain original dimensions.',
'new_biological_tests':0,'all41_table_values_unchanged':True,'original_pdf_unchanged':True}
(P/'qc/FULL_MANUSCRIPT_REVIEW.json').write_text(json.dumps(review,ensure_ascii=False,indent=2))
complete={**v,'status':'PASS','visual_status':'PASS','logic_review':'PASS','completed_at':now,'sha256':sha,'bytes':pdf.stat().st_size}
(P/'qc/COMPLETE.json').write_text(json.dumps(complete,ensure_ascii=False,indent=2))
delivery=P/'delivery'/pdf.name;shutil.copy2(pdf,delivery)
record=f"""## Phase326（{now[:10]}）：免疫章节已正式并入完整文章并完成全文浏览
用户指定Windows微信目录的26210_manuscript_clinical_suppfig_v1.pdf，其SHA256与Phase273可编辑底稿完全一致（25e580c0b0f2e56567d6f431012fc425c50188c832d06fefb5e7004048311112），授权插入合适位置、沿用章节编号并仔细检查全文逻辑。
已将Phase325免疫机制解释稿并入正式全文。详细方法按分析类型分布并完整保留；新增结果2.13–2.16位于原C4状态之后，原临床合并为2.17，空间结果续为2.18–2.26，其中2.22为SPP1/基质及外部功能证据，2.23合并T/CD8界面。临床与空间重复内容去重；摘要、研究问题、讨论、工作模型与结论同步。
全文52页，方法1.1–1.25、结果2.1–2.26、讨论3.1–3.12、附表说明4.1–4.16。原34表与新增7表共41表数值不变；主表1–13，既有附表S1–S20保持；35文献统一映射；原5主图+S1/S2共7矢量图原尺寸、源哈希不变；68个图件链接及PDF书签有效。
52页全部视觉浏览，末次修改7页[2,3,8,9,22,32,46]复查，其余45页与已审图像逐字节一致。全文逻辑核查重点完成C4群间差异与TYMP细胞内关联、肿瘤来源与髓系来源、表达程序阳性与旧丰度邻近阴性、空间分区与特异屏障的区分。完整证据：{P}/qc/FULL_MANUSCRIPT_REVIEW.json及COMPLETE.json。
当前完整文章采用{delivery}，SHA256 {sha}。原Phase273、Phase325和所有历史稿保留，无新生物学分析。当前断点为整合全文审阅；旧“尚未插入完整稿”断点已完成。
"""
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md']:
 f=R/name;f.write_text(record+'\n'+f.read_text())
with (R/'DECISIONS.md').open('a') as f:f.write('\n'+record)
manifest=R/'RESULTS_MANIFEST.tsv'
with manifest.open() as f:header=next(csv.reader(f,delimiter='\t'))
assert header==['result_id','stage','artifact_type','absolute_path','generated_at','validation_status','validation_evidence','notes']
with manifest.open('a',newline='') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n')
 for rid,kind,path in [('phase326_integrated_full_pdf','full_manuscript',delivery),('phase326_content_source','editable_report_source',P/'report/report_materialized.json'),('phase326_full_review','logic_and_visual_review',P/'qc/FULL_MANUSCRIPT_REVIEW.json'),('phase326_complete','final_validation',P/'qc/COMPLETE.json')]:
  w.writerow([rid,'phase326_full_manuscript_immune_integration',kind,str(path),now,'PASS',str(P/'qc/COMPLETE.json'),'52 pages; original34+immune7 tables; 7 native figures; no new statistical tests'])
(P/'qc/DELIVERY_MANIFEST.json').write_text(json.dumps({'filename':delivery.name,'sha256':sha,'bytes':delivery.stat().st_size},indent=2))
print(json.dumps({'status':'PASS','filename':delivery.name,'sha256':sha,'bytes':delivery.stat().st_size,'pages':52}))
