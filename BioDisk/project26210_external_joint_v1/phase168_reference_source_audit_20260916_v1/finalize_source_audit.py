from pathlib import Path
import json,csv,hashlib,datetime,shutil
import pandas as pd
O=Path(__file__).resolve().parent;R=O.parent
status=json.loads((O/'MAPPING_VALIDATION.json').read_text())
m=pd.read_csv(O/'SPATIAL_26_SECTION_PATIENT_MAP.tsv',sep='\t',keep_default_na=False)
assert len(m)==26 and m.section.is_unique
assert len(m[m.dataset=='GSE175540'].patient_id.unique())==23
assert len(m[m.dataset=='GSE210041'].patient_id.unique())==2
tree=json.loads((O/'author_repo_tree.json').read_text());log=json.loads((O/'AUTHOR_API_FETCH_LOG.json').read_text())
assert len(log)==33 and all('error' not in r for r in log)
notes={
 'scienceadv_reference':'10.1126/sciadv.aed0733; user supplied published PDF, Methods page14 and Data availability page19',
 'combined_pdf':'/mnt/BioDisk/wangmingyuan/aed0733_CombinedPDF_v4 (1).pdf',
 'combined_version':'72-page prepublication combined manuscript; final supplement pages69-70 visually reviewed; not asserted identical to final publisher supplement',
 'original_source':'10.1016/j.immuni.2022.02.001; official Elsevier supplementary Table S1 (mmc2) and S2 (mmc3), Figure S1/S4 (mmc1)',
 'author_code_commit':tree['sha'],
 'zenodo_archive':'Record obtained; full ZIP retrieval stopped after slow transfer; 33 relevant code/text/table files fetched successfully through GitHub API at v1 tree',
 'formal_sciadv_supplement':'publisher HTTP403; preserved response is not a valid PDF',
 'mapping_v1_failure':'Earlier GSM hits on pages40/59 are not duplicate clinical tables; corrected to verified table at69-70, no change to source values',
 'images':'Only read-only source-page renderings; no new analysis plots',
 'unknowns':['GSE175540 prior treatment/clinical stage not provided in mapped per-sample table','GSM5924036 pT missing','TLS Table S2 a_4 not a_3: no silent relabel','Cross-dataset patient overlap not tested; do not combine dataset patient counts blindly']
}
(O/'AUDIT_NOTES.json').write_text(json.dumps(notes,ensure_ascii=False,indent=2))
report="""# 原文资料核查：24张GSE175540切片对应23位患者，原始补充表已确认

## 核心结果
用户参考的Science Advances论文为Ge等，2026，10.1126/sciadv.aed0733。其正式正文第14页列出GSE175540（24张）和GSE210041（2张）。本轮沿参考论文、Zenodo/作者代码、原始GEO、原始Immunity研究补充材料追溯，已经解决Phase165–167留下的GSE175540患者独立性问题。

原始Immunity研究官方补充表S1（original_mmc2.xlsx，Sheet1）第74行：ID c_2（B74），队列BIONIKK（A74），pT1a（C74），空间制备Frozen / FFPE（K74）。结合GEO精确标题，GSM5924030=ffpe_c_2与GSM5924050=frozen_c_2是同一研究对象的两种制备。其余22个空间ID各对应一条原始对象记录；全部24切片的ID、制备都与原始表一对一/一对二匹配。因此GSE175540计23位独立研究对象，不能按24位患者计算。

| 队列 | 切片 | 独立对象 | 制备 |
|---|---:|---:|---|
| ExhauCRF | 4 | 4 | frozen |
| IMM | 4 | 4 | frozen |
| BIONIKK | 16 | 15 | 12 FFPE + 4 frozen |
| GSE175540合计 | 24 | 23 | 12 FFPE + 12 frozen |

## 对Phase167结果的更新
Phase167冻结规则与技术输出保持原样，患者映射以本轮为准：
- 上游可用24张切片，去重后23位患者。
- 邻域接收端可用20张切片，去重后19位患者。
- c_2的两张切片均通过上游及接收端，后续不能把它们当两个独立重复。
- FFPE和frozen分别分析时各包含12位对象，但两组共享c_2，不能称完全独立的两个人群。
- 可在冻结分析方案后开展患者层方向/效应量与推断性统计；必须先明确同患者双制备的汇总或配对处理及队列/制备差异。不要为增强关联而任选其中一张。
- 本轮仅补全来源、临床标签与独立样本数；未运行关联、P/FDR、中介或生物学新图。原Phase167的24/20技术可用结论有效，患者未解决的边界已被本轮证据更新。

## 新增临床信息
GSE175540原始补充表提供23/24切片的pT标签（22/23独立对象）；GSM5924036=c_21为NA。这些标签与72页合并稿补充表2的26行临床信息全部一致。T分期与总体临床分期分开保留；pT3a不能自动填成临床Stage III。该表中GSE175540的Stage、Treatment均为NA，不能把缺失治疗填作未治疗，也不能从BIONIKK队列名称推断具体样本用药时点/方案。

GSE210041从GEO原始SOFT获得：
| GSM | 原始患者ID | 临床分期 | ISUP | 既往治疗 |
|---|---|---|---|---|
| GSM6415705 | Patient 1 | Stage IV | ISUP3 | Nivolumab+Ipilimumab |
| GSM6415706 | Patient 7 | Stage III | ISUP2 | Naive |

两例都为FFPE。临床分期/治疗与合并稿补充表一致。它们尚未进入Phase167检查，不自动扩充为26个可用上游或22个接收端；两个数据集的ID单独命名，跨数据集是否重叠未核实，不简单宣称总计25位独立患者。

## 原始实验设计与TLS信息
官方Figure S1显示原始130例手术切除原发ccRCC肿瘤来自ExhauCRF46、IMM25、BIONIKK59；经TLS相关选择，12 frozen+12 FFPE进入空间测序。当前24张是选择后的空间子集，不能当无选择的普通ccRCC人群。
官方Table S2（original_mmc3.xlsx）含12个Visium BCR分析ID及TLS状态。与当前12张frozen精确匹配11个ID；表中的a_4没有对应当前frozen条目，当前a_3未在该表中出现。保留这一不一致，不将a_4静默替换成a_3。官方Figure S4图注明确将a_3列作TLS阳性，可作为单独来源证据；本轮TLS_SOURCE_REVIEW保留精确表匹配与图注线索分开。
TLS maturation及其他实验列的x表示进行了该分析，不直接改写成TLS阳性标签。

## 参考论文和代码审查
72页合并稿第69–70页的26样本临床表已视觉核对。它是已有较早版本，正文一处把另一空间数据写作GSE210042，补充表及正式出版PDF写GSE210041；旧稿GSE256632与正式正文GSE156632也有差异。数据编号采用正式正文与GEO核对结果，不以旧稿错字重新下载错误队列。
Zenodo记录19403993确认对应作者仓库qgggt61-commits/001 v1。本轮成功读取33个代码/文本/表文件。Figure2脚本使用sample_list、samples_dir、sample_names，但文件未定义完整列表；当前合并代码只合并sample_objects[[1]]与[[2]]。Figure3脚本输入路径明确指向GSM5924033，内部sample_name为UKF275T，不能把这类对象名当患者ID。公开代码可提供方法线索，不能替代本轮原始补充表的样本映射。
仓库BayesPrism结果以EA等bulk样本编号为行，图像距离表以TIB/Normal等分组为主；本轮检索的文本资料没有提供额外GSM至患者对照。未分析或重算作者湿实验。

## 来源与审计
1. Science Advances参考文章： https://doi.org/10.1126/sciadv.aed0733
2. 参考代码记录： https://zenodo.org/records/19403993
3. 原始研究官方补充表S1： https://ars.els-cdn.com/content/image/1-s2.0-S1074761322000814-mmc2.xlsx
4. 原始研究官方补充表S2： https://ars.els-cdn.com/content/image/1-s2.0-S1074761322000814-mmc3.xlsx
5. 原始研究官方补充图： https://ars.els-cdn.com/content/image/1-s2.0-S1074761322000814-mmc1.pdf
6. GEO： https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175540
7. GEO： https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210041

原始XLSX只读提取，源文件保留，患者映射逐行带Sheet1行号和单元格位置。映射26/26唯一GSM、GSE175540全部24条制备匹配、23对象/唯一c_2双制备、26行临床表与原始信息一致。源哈希、软件与失败记录另存。
Science Advances正式补充PDF入口仍HTTP403，未取得；完整Zenodo ZIP慢下载已停止，代码使用官方GitHub API锁定版本成功获取。原始Immunity补充表/补充图已完整取得并核验格式，所以关键患者结论不依赖访问失败的资源。所有分析/状态仅记录服务器，既有输入和Phase167输出未改写。

## 建议的后续分析
以本轮映射固定患者层设计：GSE175540上游23例、接收端19例；保留制备及三来源队列标识，预先规定c_2双制备处理，再开展效应量与统计。不把GSE210041两例直接并入发现队列，可单独技术核查后作为描述性补充。MTA/多胺轴、MTAP背景、SLC16A3并行及整体巨噬/C4定位规则继续适用。
"""
(O/'REFERENCE_SOURCE_AUDIT_CN.md').write_text(report,encoding='utf-8')
now=datetime.datetime.now().astimezone().isoformat();status.update(status='PASS_PATIENT_MAPPING_RESOLVED',completed_at=now,author_code_files_reviewed=33)
(O/'STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
src=[R.parent/'aed0733_CombinedPDF_v4 (1).pdf',R/'phase165_A_input_confirmation_20260915_v1/metadata/SAMPLE_METADATA_EXTRACT.json',R/'phase167_A_24section_endpoint_check_20260916_v1/SECTION_ENDPOINT_SUMMARY.tsv']
pd.DataFrame([dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in src]).to_csv(O/'INHERITED_SOURCE_HASHES.tsv',sep='\t',index=False)
valid=[dict(check='26_unique_GSM',passed=len(m)==26 and m.section.is_unique),dict(check='23_GSE175540_patients',passed=status['GSE175540_patients']==23),dict(check='c2_only_duplicate',passed=set(m[m.patient_id.duplicated(keep=False)].source_id)=={'c_2'}),dict(check='24_source_sections_23_patients',passed=status['phase167_source_patients']==23),dict(check='20_receiver_sections_19_patients',passed=status['phase167_receiver_patients']==19),dict(check='26_clinical_rows_reconciled',passed=status['clinical_table_26_rows_verified']),dict(check='33_code_text_files_downloaded',passed=len(log)==33 and all('error' not in r for r in log))]
pd.DataFrame(valid).to_csv(O/'VALIDATION.tsv',sep='\t',index=False);assert all(x['passed'] for x in valid)
hashes=[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(O.rglob('*')) if p.is_file() and 'authority_backups' not in p.parts and p.name not in ['ARTIFACT_HASHES.tsv','FINAL_STATE_VALIDATION.json']]
pd.DataFrame(hashes).to_csv(O/'ARTIFACT_HASHES.tsv',sep='\t',index=False)
B=O/'authority_backups';B.mkdir(exist_ok=False)
names=['PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','CODEX_HANDOFF.md']
for n in names:shutil.copy2(R/n,B/n)
notice=f"""
## Phase168 原始补充表解决空间患者映射（{now}）
用户参考Ge et al. Sci Adv 2026 eaed0733。已追溯官方Immunity2022补充表S1：GSE175540共24切片=23对象，c_2（Sheet1 B74/K74）明确Frozen / FFPE，对应GSM5924030/5924050。ExhauCRF4、IMM4、BIONIKK15对象（16片）；本轮取代Phase165–167患者未解决的断点。
Phase167技术资格不变，现为上游24片/23患者、接收20片/19患者。双制备按同患者处理，先冻结患者汇总方案再检验；不按24/20独立患者。GSE210041两片分别Patient1/Patient7，GEO有Stage IV/III、ISUP3/2、Nivo+Ipi/Naive，尚未进入技术检查，不盲目跨数据集合计患者。
补齐GSE175540 pT23/24片；c_21缺失。原始TLS表a_4与当前a_3不一致保留，FigureS4另有a_3阳性说明。原始图S1确认TLS选择子集；不据队列名推断每例用药。无新关联/图。
报告：{O/'REFERENCE_SOURCE_AUDIT_CN.md'}
映射：{O/'SPATIAL_26_SECTION_PATIENT_MAP.tsv'}
"""
for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']:
 p=R/n;t=p.with_name(n+'.phase168.tmp');t.write_text(notice+p.read_text(),encoding='utf-8');t.replace(p)
p=R/'RESULTS_MANIFEST.tsv'
with p.open(newline='') as f:rows=list(csv.reader(f,delimiter='\t'))
assert len(rows[0])==8 and not any(x and x[0].startswith('phase168_') for x in rows[1:])
reg=['REFERENCE_SOURCE_AUDIT_CN.md','SPATIAL_26_SECTION_PATIENT_MAP.tsv','PHASE167_ELIGIBILITY_WITH_PATIENTS.tsv','SCIADV_COMBINED_TABLE2_26_ROWS.tsv','FROZEN_TLS_SOURCE_REVIEW.tsv','original_mmc1.pdf','original_mmc2.xlsx','original_mmc3.xlsx','ORIGINAL_TABLES_EXTRACT.json','GSE210041_family.soft','MAPPING_VALIDATION.json','VALIDATION.tsv','AUDIT_NOTES.json','STATUS.json','ARTIFACT_HASHES.tsv','build_patient_map_v2.py','finalize_source_audit.py']
for n in reg:rows.append(['phase168_'+n,'Phase168','reference_source_audit',str(O/n),now,'PASS','Original-source mapping: 24 sections/23 patients; receiver20/19','sha256='+sha(O/n)])
t=p.with_name(p.name+'.phase168.tmp')
with t.open('w',newline='') as f:csv.writer(f,delimiter='\t',lineterminator='\n').writerows(rows)
t.replace(p)
gate=dict(stage='Phase168',manifest_entries=len(reg),all_registered_paths_exist=all((O/n).exists() for n in reg),four_backups=True,three_state_headers_updated=all((R/n).read_text().startswith(notice) for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']),authority_sha256={n:sha(R/n) for n in names})
assert gate['all_registered_paths_exist'] and gate['three_state_headers_updated']
(O/'FINAL_STATE_VALIDATION.json').write_text(json.dumps(gate,ensure_ascii=False,indent=2))
print(json.dumps(status,ensure_ascii=False,indent=2));print(json.dumps(gate,ensure_ascii=False,indent=2))
