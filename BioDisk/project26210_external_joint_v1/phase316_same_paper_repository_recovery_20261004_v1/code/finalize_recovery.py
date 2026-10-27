from pathlib import Path
import json,hashlib,datetime,csv,zipfile,html,shutil,subprocess
import pandas as pd
P=Path(__file__).resolve().parents[1];R=P.parent;B=R/'phase315_CRC_SPP1_barrier_paper_reproduction_20261004_v1';now=datetime.datetime.now().astimezone().isoformat()
subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True)
ident=json.loads((P/'qc/P7_IDENTITY_VALIDATION.json').read_text());assert ident['status']=='PASS'
rcheck=pd.read_csv(P/'qc/R_INDEPENDENT_CHECKS.tsv',sep='\t');assert rcheck['pass'].all() and len(rcheck)==100
base=json.loads((P/'qc/BASELINE_PHASE315_SHA256.json').read_text());changed=[n for n,v in base.items() if hashlib.sha256((B/n).read_bytes()).hexdigest()!=v];assert not changed
result=json.loads((P/'qc/RECOVERY_DECISION.json').read_text());assert result['P7_primary_direction_pass'] and result['R3_four_patient_direction_pass']
validation={'status':'PASS','time':now,'identity_checks':len(ident['checks']),'independent_R_checks':100,'baseline_hashes_unchanged':len(base),'baseline_changes':changed,'candidate_official_MD5':'PASS','coordinates_labels_UMI_features_exact_match':4248,'raw_count_normalization_max_error':ident['normalization_max_error'],'scientific_scope':'Same-paper P7 missing-input recovery and frozen descriptive immune readout; not new cohort validation; original Pearson/Spearman reporting mismatch retained'}
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(validation,indent=2))
search=[{'source':'Zenodo12634581 / UCASpatial','url':'https://zenodo.org/records/12634581','status':'P7_RECOVERED_VERIFIED','detail':'Official CRC spatial example, MD5 verified; exact 4248 coordinates, labels, UMI and feature counts match original P7; internal orig_ident P2 is mirror name.'},{'source':'Author GitHub','url':'https://github.com/youqiongye/CRC_scRNAseq','status':'NO_NEW_SPATIAL_FILES_OR_CORRECTION','detail':'Current main same commit tree as V1.0.0:33fe60d146e5493bdc1485f1ae4ebda26d4d6d5e; issues API returned empty.'},{'source':'SpatialTME author portal','url':'https://www.spatialtme.yelab.site/','status':'NO_DATA_USED','detail':'Homepage/JS/GitHub accessible; HRA-filtered sample response empty; all-sample request TLS timeout. Did not infer direct access to server-side data directories.'},{'source':'Zenodo14602110','url':'https://zenodo.org/records/14602110','status':'NOT_SELECTED','detail':'Crescendo dataset; metadata did not establish same-paper P7 identity; no expression data downloaded.'},{'source':'UCASpatial Figure34.zip','url':'https://zenodo.org/records/12634581','status':'DIRECTORY_AND_CODE_ONLY','detail':'25,015,791,543 byte public archive; 20,953 bytes range read for directory, 109,599 byte code region retrieved; 180 member directory and selected code CRC verified; full archive not downloaded or whole-file MD5 validated.'},{'source':'Zenodo18482631','url':'https://zenodo.org/records/18482631','status':'METADATA_ONLY','detail':'10.26GB additional archive identified, not downloaded once P7 recovered.'}]
(P/'resources/REPOSITORY_SEARCH_FINDINGS.json').write_text(json.dumps(search,indent=2))
rows=pd.read_csv(P/'tables/P7_RECOVERED_IMMUNE_RESULTS.tsv',sep='\t')
text='\n'.join(['项目26210 Phase316：原文P7数据已从公开仓库补齐','',f'完成时间：{now}',
'一、结果结论',
'本次补全成功。UCASpatial官方Zenodo仓库中的CRC示例对象，经逐spot身份核验确认就是Qi等2022年SPP1/FAP屏障论文的P7切片。已用真实表达矩阵补齐Phase315缺失的免疫空间读出；未使用其他患者、插补值或新边界替代。',
'补齐P7后，原文四位患者均呈现SPP1/FAP联合状态评分正相关、ECM相关程序富集，以及肿瘤区较低T细胞转录信号的方向。前两层继承Phase315原结果，本轮只新增P7免疫描述。',
'需要区分：数据补全与生物学方向重复已成功；原文P7图值0.49与Pearson图注的差异仍存在。Pearson=0.439746、Spearman=0.485071，后者按两位小数对应0.49。没有改动原冻结统计规则，不能宣称整篇论文所有方法和图值逐项一致。','',
'二、同一样本的证据',
'来源：https://zenodo.org/records/12634581',
'文件：Example_Data_Human_CRC_ST_data.rds；129,839,528字节。',
'官方MD5：5e9b2b225ac5c6b30bce99b43e073209，实测一致。SHA256：044e8988c55c74141d32d35f777b26c36cdab9c9fafb15b106c1db3e44062d98。',
'对象维度：36,601基因×4,248 spots，12,042,958非零计数。镜像内部orig_ident=P2，未按名字直接认定样本。',
'使用对象自带lowres缩放系数0.05905512，无拟合坐标变换。4,248个坐标与原文Figure6e逐项完全一致，坐标误差0且一对一；4,248个原DefineTypes标签全部一致；原Figure8d-f的UMI总量及检出基因数也均4,248/4,248完全一致。',
'原始矩阵列和与对象元数据一致；所有目标计数非负整数。对象保存的normalized data与log1p(count/library_total×10000)逐项一致，最大误差5.33e-15。',
'P7原QC工作表实际叫Figure 8d-f，且坐标列名是imageRow/imageCol。首次验证器仅搜索Supplementary Figure8，导致未找到QC表。保留失败记录后只修复表名/大小写映射，未放宽门槛或改科学结果。最初系统R缺少SeuratObject，改用已有r_seurat环境成功读取，未安装或修改环境。','',
'三、P7补齐结果',
'沿用Phase315冻结比较：作者Malignant epithelial cells区域与其余spot比较，以CD3D均值和检出比例为主要描述；CD8A等7个原图标记完整保留。1,657个肿瘤区spot、2,591个其他spot；它们来自1位患者，不作为独立生物学重复。',
'标记\t肿瘤区均值\t其他区域均值\t均值差\t肿瘤区检出比例\t其他区域检出比例'])+'\n'
for a in rows[rows.comparison.eq('all_nonmalignant')].itertuples():text+=f'{a.gene}\t{a.mean_malignant:.9f}\t{a.mean_other:.9f}\t{a.delta_mean:.9f}\t{a.detection_malignant:.9f}\t{a.detection_other:.9f}\n'
text+='\nCD3D均值0.139256对0.251325，检出17.44%对24.31%；CD8A均值0.120908对0.182061，检出15.39%对18.14%。全部8个原图标记的均值和检出比例均较低。固定排除Unknown的敏感性保留，主要CD3D/CD8A方向不变。\n'
text+='\n四、四患者证据更新与边界\nP6/P7/P8/P9联合区Pearson r分别0.580833/0.439746/0.559819/0.553774；原40项GO原始P精确复算结果保持。新增P7后CD3D/CD8A肿瘤区平均表达与检出比例4/4均低于其他区域，仍按4个独立患者的探索性方向描述，无新患者总体P/FDR。\n这是同一论文数据的补全与复算，不是新的独立外部队列。逐spot低T信号与屏障解释方向一致，不等于直接测量物理阻挡、迁移或SPP1处理效应。原文AddModuleScore和全GO检验家族未由原始矩阵重建，Phase315的证据层级保持。P7 Pearson/Spearman图注差异继续显式保留。未进入外部肾癌验证或TYMP定位，也未重写正式文章。\n'
text+='\n五、保存、运行与验收\n服务器新目录：'+str(P)+'\n原Phase315的98个文件SHA256全部保持不变；本轮新增结果在Phase316，未覆盖冻结结果。\n'
text+='主要输入：downloads/Example_Data_Human_CRC_ST_data.rds、Phase315原source.xlsx与作者坐标标签表。新建脚本、协议、镜像清单、身份对照、16行P7结果和64行四患者结果均保存。\n实际命令：r_seurat/bin/Rscript --vanilla code/extract_candidate_identity.R；batch_bench/bin/python code/verify_P7_identity.py；同Python运行code/complete_P7_immune_readout.py；r_seurat/bin/Rscript --vanilla code/validate_P7_recovery.R。以上环境根为/mnt/BioDisk/wangmingyuan/miniconda3/envs。每个分析阶段先执行项目50GiB内存门禁。\n'
text+='环境：r_seurat（R4.1.3，具体包版本见resources/VALIDATION_R_SESSION.txt）；batch_bench（Python、pandas、SciPy沿Phase315环境）。\n17项身份核验通过；独立R直接重读原RDS与原XLSX，以原始counts重新归一化、汇总和核对，100项检查通过。正式提取、身份核验、补算及独立R验收均exit0。内存门禁通过，实际可用约376.6–379.2GiB；Swap仅余240kB已记日志。\n'
text+='根目录PROJECT_STATUS.md、CODEX_HANDOFF.md、DECISIONS.md、RESULTS_MANIFEST.tsv同步更新；Phase315保留原始部分复现状态作为历史，Phase316登记缺失数据恢复成功及尚存统计标注差异。\n'
text+='\n六、轻量结果包\n包内提供报告、16/64行结果、完整复现脚本、目标基因计数和标准化表达、坐标/原标签、官方元数据与校验记录。129.8MB完整RDS留在服务器，可用包内来源链接和官方MD5获取；不把25GB归档列为已下载。引用来源见resources/REPOSITORY_SEARCH_FINDINGS.json。\n'
name='26210_Phase316_P7数据补全与复现更新'
(P/'delivery'/f'{name}.txt').write_text(text,encoding='utf-8');(P/'delivery'/f'{name}.html').write_text('<!doctype html><html lang="zh"><meta charset="utf-8"><title>Phase316 P7数据补全</title><style>body{max-width:1080px;margin:36px auto;font:16px/1.65 system-ui;color:#17212b;padding:0 20px}pre{white-space:pre-wrap;overflow-wrap:anywhere;font:inherit}</style><body><pre>'+html.escape(text)+'</pre></body></html>',encoding='utf-8')
state={'status':'SAME_PAPER_P7_RECOVERED_VALIDATED','time':now,'recovery_success':True,'R3_four_patient_direction_pass':True,'exact_Pearson_plot_match_unresolved':True,'external_validation_started':False,'TYMP_positioning_started':False,'phase315_baseline_preserved':True,'analysis_jobs':'none','validation':validation}
(P/'STATUS.json').write_text(json.dumps(state,indent=2))
for name0 in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv','AGENTS.md']:shutil.copy2(R/name0,P/'root_record_backups'/name0)
note='## Phase316终态（2026-10-05）：同论文P7数据已从公开镜像补齐，四患者方向一致\n用户授权仅搜索其他仓库补全，失败则恢复Phase315。已在UCASpatial官方Zenodo12634581找到Example_Data_Human_CRC_ST_data.rds（129839528字节，官方MD5匹配）。内部名字P2但4248坐标/原区域标签/UMI总量/检出基因数与原文P7全部逐项一致，确认同一切片。补齐真实计数与归一化表达，不是其他患者或插补。\nP7按原冻结区域比较，CD3D肿瘤0.139256/其他0.251325、检出17.44%/24.31%；CD8A0.120908/0.182061、检出15.39%/18.14%。排Unknown同向。合并既有三例后CD3D/CD8A均值及检出4/4肿瘤区较低；原SPP1/FAP4/4正向及40项GO复算继承，不重算原结果。无新增生物学P/FDR。\n数据补全成功，17项身份/100项独立R验收PASS，Phase315全部98文件哈希未变。原P7图0.49对应Spearman0.485071、Pearson0.439746与图注差异仍保留；不能改称严格方法逐项完全复现。无外部验证、无TYMP定位、无正式稿改写。目录phase316_same_paper_repository_recovery_20261004_v1；交付delivery/26210_Phase316_P7数据补全与复现更新.txt及轻量ZIP。系统R缺包和源QC表名差异造成的早期失败均留档，正式阶段均exit0；无计算后台。\n\n'
for n in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:
 f=R/n;f.write_text(note+f.read_text())
f=R/'AGENTS.md';f.write_text(f.read_text()+'\n\n## 2026-10-05 同论文数据补全范围\n用户授权搜索其他仓库补齐原文数据，失败则回到现有复现。Phase316已成功恢复并严格验证P7，Phase315原文件全部保留。当前仅补全同论文数据与原固定读出；不因补全自动启动外部队列或TYMP定位。P7原论文Pearson/Spearman标注差异继续保留，不改旧统计判据。\n')
with (R/'RESULTS_MANIFEST.tsv').open('a',newline='') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n')
 for rid,typ,path in [('phase316_P7_identity','identity_validation','qc/P7_IDENTITY_VALIDATION.json'),('phase316_P7_immune','same_paper_recovery','tables/P7_RECOVERED_IMMUNE_RESULTS.tsv'),('phase316_four_patient_immune','descriptive','tables/FOUR_PATIENT_IMMUNE_RESULTS.tsv'),('phase316_report','report',f'delivery/{name}.txt')]:w.writerow([rid,'Phase316',typ,str(P/path),now,'RECOVERY_PASS;EXACT_METHOD_DISCREPANCY_RETAINED',str(P/'qc/FINAL_VALIDATION.json'),'17 identity and 100 independent R checks;Phase315 98 hashes unchanged;no external validation'])
files=[f for f in P.rglob('*') if f.is_file() and not any(k in f.parts for k in ['root_record_backups','downloads']) and f.suffix not in ['.zip','.bin'] and f.name not in ['DELIVERY_VALIDATION.json','PACKAGE_FILES.json']]
manifest={str(f.relative_to(P)):{'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in files};(P/'delivery/PACKAGE_FILES.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False));files.append(P/'delivery/PACKAGE_FILES.json')
archive=P/'delivery/26210_Phase316_P7_recovery_results.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in files:z.write(f,str(f.relative_to(P)))
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for n,v in manifest.items():assert hashlib.sha256(z.read(n)).hexdigest()==v['sha256']
out={'status':'PASS','zip':str(archive),'bytes':archive.stat().st_size,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'report_sha256':hashlib.sha256((P/'delivery'/f'{name}.txt').read_bytes()).hexdigest(),'files':len(files),'all_member_CRC_SHA':'PASS','full_RDS_stored_server_only':True}
(P/'qc/DELIVERY_VALIDATION.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
