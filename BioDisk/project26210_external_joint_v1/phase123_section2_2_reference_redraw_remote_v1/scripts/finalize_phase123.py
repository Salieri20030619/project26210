from pathlib import Path
import json,hashlib,datetime,shutil,subprocess
w=Path(__file__).resolve().parents[1];root=w.parent;qc=w/'qc'
val=json.loads((qc/'VALIDATION.json').read_text());ind=json.loads((qc/'INDEPENDENT_VALIDATION.json').read_text());ora=json.loads((qc/'ORA_INDEPENDENT_RECHECK.json').read_text());delivery=json.loads((qc/'DELIVERY_VALIDATION.json').read_text());rows=json.loads((qc/'DELIVERY_MANIFEST.json').read_text())
assert ind['status']=='PASS' and ora['status']=='PASS' and delivery['sha256_matches']==21
assert val['CNV_patient_n']==19 and val['clone_DE_ORA_up_FDR_significant']==17
checks=[]
for row in rows:
 f=w/'figures'/row['name'];assert hashlib.sha256(f.read_bytes()).hexdigest()==row['sha256']
 if f.suffix=='.pdf':
  p=subprocess.run(['pdfinfo',str(f)],capture_output=True,text=True);assert p.returncode==0 and any(l.startswith('Pages:') and l.split(':')[1].strip()=='1' for l in p.stdout.splitlines());checks.append({'file':f.name,'pages':1,'status':'PASS'})
(qc/'PDF_VALIDATION.json').write_text(json.dumps({'status':'PASS','tool':'pdfinfo','checks':checks},indent=2))
visual={'status':'PASS_VISUAL_REVIEW','all_ten_families':'final contact sheet reviewed for layout, labels, legends and clipping','full_resolution_review':['03_E_program_gene_expression_heatmap_v1.png','05_I_gene_contribution_reference_blue_v2.png','10_F1_high_clone_up_pathway_enrichment_v1.png'],'I_final':'v2 corrected annotation bracket outside grid; prior v1 archived','user_review':'PENDING'}
(qc/'VISUAL_REVIEW.json').write_text(json.dumps(visual,indent=2))
manifest=json.loads((qc/'SOURCE_MANIFEST.json').read_text());known={r['path'] for r in manifest};base=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1')
more=[base/'02_integration_annotation/cell_annotation_v1/part1_scvi_integrated_cell_annotated_v1.h5ad',base/'05_cnmf/full_ranks4_8_v1/full_models_v1.rds',base/'05_cnmf/full_ranks4_8_v1/run_qc_v1.tsv',root/'local_records_migration_20260908_150455/references/image9.png',root/'local_records_migration_20260908_150455/references/image10.png']
for f in more:
 if str(f) in known:continue
 h=hashlib.sha256()
 with f.open('rb') as rd:
  for b in iter(lambda:rd.read(1048576),b''):h.update(b)
 manifest.append({'path':str(f),'bytes':f.stat().st_size,'sha256':h.hexdigest()})
(qc/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
backup=w/'state_backup_before_completion';backup.mkdir(exist_ok=False)
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv','AGENTS.md']:shutil.copy2(root/name,backup/name)
section='## Phase123 section2.2 cNMF reference redraw completed (2026-09-08)\n\n- 用户明确要求进入原恶性上皮章节2.2，按原方案参考image9及image10重画。当前交付为第一版图件复核批次，下一断点：用户审阅2.2，不自动进入2.3。\n- 服务器唯一权威目录：`/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase123_section2_2_reference_redraw_remote_v1`。所有进度、脚本、输入身份、统计表及审计只保存在服务器。\n- 交付：`C:/Users/31937/Downloads/26210/03_2.2_cNMF模块与TYMP_v1`。10 PNG + 10 PDF + 中文阅读说明，正式文件21/21 SHA256匹配；I图最终为v2，两张I v1保存在历史草稿子目录。\n- 图件：rank QC；冻结UMAP；E模块RNA热图；H真实140分量共识相关性及聚类树；I指定蓝色基因贡献矩阵；模块队列构成；患者层F1–TYMP和F1–CNV；模块top100载荷功能富集；F1高克隆配对DE；F1高克隆上调基因集富集。\n- 复用全部11,950个正式恶性细胞、2,500个cNMF输入基因、rank7及20次原始运行。未重新拟合cNMF/UMAP/克隆。140分量复建共识载荷误差4.51e-16；20位患者效应复算最大误差4.72e-16；PDF 10/10单页结构通过，最终视觉检查完成。\n- F1–TYMP：20位患者17位正向，中位rho=0.1149743，后续18终点双侧审核P=0.0153122，BH FDR=0.0196871。F1主导细胞4,079个，99.6568%来自GSE207493。F1–CNV仅19位可估计患者，8位正向，中位rho=-0.0332641，原rank筛选单侧P=0.6782408、FDR=0.9234400。两面板保留不同来源检验，不合并校正。\n- **重要历史更正：旧Phase116/117和后续交接中“F1-high clone pathway ORA全部阴性/全部FDR=1”的概括错误。原始表与Phase116副本SHA256一致，上调17/201条FDR显著（min q=0.0003101934），下调0/201且q均为1。原脚本取各方向按edgeR P排序前250基因，非只取FDR显著基因；背景15,308基因，两个方向各201条通路独立BH。按原定义独立复算误差<=1.89e-15。该更正取代旧阴性摘要；未改变原DE/ORA结果。**\n- E是RNA表达而非VIPER活性；I是相对基因载荷而非DE统计。TYMP不在原cNMF 2,500 HVGs内，不伪造TYMP loading；其表达单独取自原RNA矩阵。F1筛选使用TYMP证据，本批不构成新独立生物学验证，也不建立TYMP驱动更恶性状态的因果性。\n- 完整审核见 `qc/VALIDATION.json`、`INDEPENDENT_VALIDATION.json`、`ORA_INDEPENDENT_RECHECK.json`、`HISTORICAL_SUMMARY_CORRECTION.md`、`SOURCE_MANIFEST.json`、`DELIVERY_MANIFEST.json`、`VISUAL_REVIEW.json`。\n'
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:
 f=root/name;f.write_text(section+'\n\n'+f.read_text(encoding='utf8'),encoding='utf8')
f=root/'AGENTS.md';text=f.read_text(encoding='utf8');text='\n'.join('当前任务与断点以服务器 CODEX_HANDOFF.md 最新阶段为准。' if l.startswith('当前任务：') else l for l in text.splitlines());f.write_text(text+'\n',encoding='utf8')
with (root/'RESULTS_MANIFEST.tsv').open('a') as mf:
 for i,row in enumerate(rows):mf.write('phase123_delivery_'+str(i+1)+'\tphase123_section22_redraw\treview_artifact\t'+str(w/'figures'/row['name'])+'\t'+datetime.datetime.now().isoformat()+'\tPASS_REVIEW_CHECKS\tSHA256 '+row['sha256']+'\tFrozen source redraw; historical ORA summary corrected\n')
 mf.write('phase123_ORA_correction\tphase123_section22_redraw\tevidence_correction\t'+str(qc/'HISTORICAL_SUMMARY_CORRECTION.md')+'\t'+datetime.datetime.now().isoformat()+'\tPASS_INDEPENDENT_RECHECK\tUp 17/201; down 0/201; original file hashes identical\tSupersedes historical all-negative ORA summary\n')
completion={'status':'PASS_COMPLETE_PENDING_USER_REVIEW','stage':123,'files_delivered':21,'source_hashes':len(manifest),'PDF_checks':len(checks),'all_primary_state_files_updated':True,'active_plot_jobs':False,'next':'User reviews section 2.2 first batch'}
(qc/'COMPLETION.json').write_text(json.dumps(completion,indent=2));(qc/'CURRENT_TASK.json').write_text(json.dumps(completion,indent=2));print(json.dumps(completion))
