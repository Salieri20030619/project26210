from pathlib import Path
import json,hashlib,datetime,csv,io,shutil,os
P=Path(__file__).resolve().parents[1];R=P.parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
now=datetime.datetime.now().astimezone().isoformat()
report=P/'delivery/文献方法对照与下一步取舍.md'
sources=[
('HCC_TIB_2023','https://doi.org/10.1016/j.jhep.2023.01.011','publisher_abstract_and_accessible_excerpts;exact_supplement_thresholds_not_retrieved'),
('CRC_NicheNet_2022','https://www.nature.com/articles/s41467-022-29366-6','Fig7_results_and_method_excerpts'),
('RCC_spatial_topography_2021','https://www.nature.com/articles/s41379-021-00864-0','methods_and_results_excerpts'),
('TNBC_spatial_classifier_2021','https://www.nature.com/articles/s41467-021-25962-0','methods_and_results_excerpts;coefficients_not_retrieved'),
('DDR1_collagen_2021','https://pmc.ncbi.nlm.nih.gov/articles/PMC8839149/','full_text_methods_SHG_CT_FIRE;27_gene_reference_not_extracted'),
('TabulaTIME_2025','https://pmc.ncbi.nlm.nih.gov/articles/PMC12643937/','results_methods_excerpts;reference_not_downloaded'),
('TYMP_pancreas_2026','https://link.springer.com/article/10.1186/s13046-026-03712-8','full_text;Results_Discussion_CD3_direction_conflict_confirmed')]
with (P/'resources/SOURCE_AUDIT.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['study','url','read_scope']);w.writerows(sources)
entry="""## Phase294（2026-09-30）：回到原始文献重审方法
用户要求从相关文献找方法。本轮只做有针对性方法对照，未新分析/绘图/改标签/文章。记录：phase294_literature_method_reassessment_20260930_v1/delivery/文献方法对照与下一步取舍.md。
重点取舍：优先接收端实际靶基因驱动的NicheNet候选输出—间质程序联系，区别SPP1分子与SPP1相关巨噬状态；TGFB1/IL1输出文献依据不替代项目结果。免疫表型应区分冷/浸润/排斥样，不能低T=屏障；借RCC组织分型及TNBC发表分类器框架，但原权重/归一化未取得前不当现成评分。独立TabulaTIME参考用于身份/程序对齐，允许unknown、不将壁细胞强改CAF。
DDR1原研究SHG/CT-FIRE胶原几何与ECM转录量不同，当前无对应图像不可直接复现；引用的27基因签名尚未追溯完整。TYMP胰腺文献Results/Discussion CD3方向冲突再次确认，不能当T进入改善依据。HCC2023精确补充阈值本轮未取得；原文方法细节未核清的项目只列候选。
所有Phase280–293正负证据与来源定义保持。下一实施需按已发表方法取原始资源并冻结新问题，不重复同批阈值搜索。审阅完成，无后台。

"""
status=dict(status='LITERATURE_METHOD_REVIEW_COMPLETE',at=now,studies=len(sources),new_analysis=False,full_systematic_review=False,report_sha256=sha(report),method_resource_gaps=['HCC supplement exact thresholds','classifier original coefficients/preprocessing','27gene signature original list','TabulaTIME reference inputs'])
(P/'qc/COMPLETION.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
updates={}
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv']:
 old=(R/name).read_bytes()
 with (P/'root_record_backups'/name).open('xb') as f:f.write(old)
 if name!='RESULTS_MANIFEST.tsv':updates[name]=entry.encode()+old
 else:
  fields=old.decode('utf-8-sig').splitlines()[0].split('\t');out=io.StringIO();w=csv.DictWriter(out,fieldnames=fields,delimiter='\t',lineterminator='\n')
  w.writerow(dict(result_id='phase294_literature_methods',stage='Phase294',artifact_type='literature_method_audit',absolute_path=str(report),generated_at=now,validation_status='REVIEWED_WITH_ACCESS_LIMITS',validation_evidence=str(P/'resources/SOURCE_AUDIT.tsv'),notes='7 primary studies;no new analysis;resource gaps explicit;SHA256='+sha(report)))
  updates[name]=old+(b'' if old.endswith(b'\n') else b'\n')+out.getvalue().encode()
for name,content in updates.items():
 target=R/name;assert target.read_bytes()==(P/'root_record_backups'/name).read_bytes()
 tmp=R/(name+'.phase294.tmp')
 with tmp.open('xb') as f:f.write(content)
 shutil.copymode(target,tmp);os.replace(tmp,target)
(P/'qc/ROOT_UPDATES.json').write_text(json.dumps({n:sha(R/n) for n in updates},indent=2))
print(json.dumps(status,ensure_ascii=False))
