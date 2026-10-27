from pathlib import Path
import json,hashlib,datetime,csv,zipfile,html,subprocess,platform,shutil
import pandas as pd,numpy,scipy,openpyxl
P=Path(__file__).resolve().parents[1];R=P.parent
subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True)
now=datetime.datetime.now().astimezone().isoformat()
v=json.loads((P/'qc/FINAL_VALIDATION.json').read_text());assert v['status']=='PASS'
d=json.loads((P/'qc/REPRODUCTION_DECISION.json').read_text());assert d['success'] is False
c=pd.read_csv(P/'tables/REPRODUCED_CORRELATIONS.tsv',sep='\t');g=pd.read_csv(P/'tables/GO_STATISTICS_REPRODUCTION.tsv',sep='\t');im=pd.read_csv(P/'tables/IMMUNE_MAP_QUANTIFICATION.tsv',sep='\t')
lines=['项目26210 Phase315 文献原数据复现结果','空间共现与ECM方向得到复现；完整复现未通过，已停止','',f'日期：{now}','原文：Qi et al. Nature Communications 2022;13:1742','https://www.nature.com/articles/s41467-022-29366-6','',
'1. 本轮判定',
'当前判定为 PARTIAL_REPRODUCTION_STOPPED。主要生物学方向在可取得的原文源数据中重复，但未满足事前冻结的完整复现标准。停止原因是原文相关方法标注与图值不一致，以及一张切片免疫源数据缺失/重复，不是检出了与文献相反的主要生物学方向。',
'用户明确要求先复现，成功后才考虑外部验证定位，失败则停止并报告。本轮未启动外部肾癌验证、TYMP定位或另一篇文章替代复现。',
'本轮使用作者公开的逐spot评分、坐标、区域标签、免疫基因表达和GO计数信息进行源数据层复算。没有从FASTQ、原始空间count矩阵重新生成细胞定义或AddModuleScore。公开Supplementary Data3/4为54,103个单细胞的矩阵与注释，不是这四张空间切片；原空间原始数据HRA000979受控，未申请访问。','',
'2. SPP1相关巨噬与FAP成纤维状态在四位患者中均正相关',
'使用原文指定联合区域，完整保留四位患者。下列为Pearson r及源表spot级原文模型P；这些P不作为独立患者推断。',
'患者\t联合区spot数\tPearson_r\t原文spot模型P\tSpearman_r（固定诊断）']
for a in c[c.scope.eq('author_joint_cluster')].itertuples():lines.append(f'{a.patient}\t{a.n_spots}\t{a.pearson_r:.9f}\t{a.paper_spot_model_P:.8g}\t{a.spearman_r_diagnostic:.9f}')
lines+=['',
'P6复算Pearson r=0.580833，与图6j标注0.58一致。P7按图注所述Pearson复算r=0.439746，与图6k的0.49不一致；事前保留的Spearman诊断为0.485071，按两位小数恰为0.49。图内标记Rs，而图注写Pearson，说明存在统计方法标注不一致。它不推翻P7正相关方向，但不能把换用Spearman当作原Pearson复现通过；冻结判据未改。',
'P6主图源表含全片4457行，按图注必须取作者联合区772行；P7主图表已经是联合区1187行。全片与联合区的全部结果分别保留，未挑取更显著的范围。',
'主图旧MARCO列名与补充源表SPP1评分的逐行差异不超过3.58466e-5，P6主图评分与带坐标补充表完全一致；作者联合区标签、标记表和发表图支持将其识别为发表的SPP1相关状态评分。保留旧列名及映射依据，不将该核对声称为从原始基因重新生成SPP1程序。','',
'3. ECM相关富集的40项原始统计复算一致',
'Figure6l四位患者共40项GO源统计，使用原GeneRatio/BgRatio重建超几何右尾与Fisher精确检验。全部原始P复算一致，SciPy与源表最大log10(P)绝对差1.42e-14；独立log-gamma组合数求和最大误差2.85e-11。',
'正文预定的ECM组织、胶原纤维组织和TGFβ反应三项，在四位患者中均为富集且原始P<0.05。',
'患者\tGO\t富集倍数\t复算原始P\t作者原FDR']
for a in g[g.selected_primary].itertuples():lines.append(f'{a.patient}\t{a.GO}\t{a.fold_enrichment:.6f}\t{a.hypergeom_P:.8g}\t{a.author_FDR:.8g}')
lines+=['',
'作者FDR原值完整保留；没有对仅展示的10项重新BH。缺少原完整差异基因背景和全部GO检验家族，本轮确认计数到原始P的数值复算，不声称完成原始表达到差异分析、GO及全家族FDR的端到端重建。GO的基因集合统计P也不是4名患者的处理效应P。','',
'4. 可评估三张切片中，肿瘤区T细胞信号较低；P7未能完成',
'原文免疫地图属于定性展示。本轮事前登记一个描述性量化：按作者Malignant epithelial cells标签，与全部其余spot比较CD3D平均表达和检出比例；CD8A等另外7个原图标记完整列出。未修改区域、距离或阈值；不计算生物学P值。',
'患者\t标记\t肿瘤区均值\t其他区域均值\t均值差\t检出比例差']
for a in im[im.comparison.eq('all_nonmalignant')&im.gene.isin(['CD3D','CD8A'])].itertuples():lines.append(f'{a.patient}\t{a.gene}\t{a.mean_malignant:.9f}\t{a.mean_other:.9f}\t{a.delta_mean:.9f}\t{a.delta_detection:.9f}')
lines+=['',
'P6、P8、P9的CD3D和CD8A平均表达与检出比例均在肿瘤区较低；固定的排除Unknown比较也保持方向。它们是3名患者的描述性组织对应，不是新迁移实验，也不证明SPP1导致T细胞排斥。',
'P7缺口：Supplementary Figure9s有4457行，与Figure6m（P6）的10列全部逐项相同。P7标签表为4248行，二者坐标交集为0。源工作簿中全部带坐标CD3D表只有Figure6m、Supplementary Figure9s/9t/9u，前三个独立空间输入实际对应P6、P8、P9，没有可核验的P7免疫表。未将P6行按顺序拼到P7，也未用图像估计或其他患者补齐。完整四患者R3因此未通过。','',
'5. 可以保留的结论与本轮停止边界',
'作者公开数据重复支持SPP1相关巨噬/FAP成纤维状态共现及联合区域ECM富集，可评估的三张切片也呈现与T细胞肿瘤内分布不足一致的方向。不能报告完整四患者、精确方法一致的屏障复现成功。',
'这是原论文源数据复算，不能计作新的独立生物学验证；同spot评分共现不等于细胞直接接触，混合区域低T信号不等于已测物理阻挡。原文机制背景及项目TYMP结果保持各自证据来源。',
'按用户停止条件，本阶段已收束。没有扩大统计条件、变更程序、替换研究或启动外部验证。若以后补齐作者P7源数据并澄清图6统计方法，才有重新开启本阶段的明确输入条件；本轮未联系作者。','',
'6. 数据、代码与验收',
'服务器目录：'+str(P),
'实际源文件：downloads/source.xlsx（34,709,129字节）；downloads/41467_2022_29366_MOESM7_ESM.xlsx（3,399,734字节）；source_file_descriptions.pdf；resources/github_tree.json；resources/author_Fig6.png。',
'作者源表SHA256：e32d1708174e0a46bc0ae2aad3c86b6e77b3189248338abf5d23b85c473a5981。',
'完整协议：resources/REPRODUCTION_PROTOCOL_FROZEN.json；结果：tables/REPRODUCED_CORRELATIONS.tsv、GO_STATISTICS_REPRODUCTION.tsv、IMMUNE_MAP_QUANTIFICATION.tsv。',
'命令：/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python code/reproduce_spatial_source.py；同环境python code/validate_reproduction.py。各脚本内部每阶段调用/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh。',
'环境：batch_bench，Python '+platform.python_version()+'；pandas '+pd.__version__+'；NumPy '+numpy.__version__+'；SciPy '+scipy.__version__+'；openpyxl '+openpyxl.__version__+'。',
'独立核验266项全部通过，包含原XLSX重新读取、手工相关/秩相关及t尾概率、均值/检出比例、GO组合数求和、坐标交集、重复表与原件/协议SHA。该PASS代表结果及停止原因核查通过，不代表完整科学复现成功。',
'正式分析与独立验证退出码0。网络首轮SSL错误与重试日志完整保留；必要源数据取得成功，作者README下载未完成，代码树元数据已取得。全片相关极小P在浮点计算中有下溢到0的情况，仅属原spot模型诊断，不能解释为数学上的P=0。',
'内存门禁均通过，MemAvailable约377.7–377.9 GiB；Swap总量约2 GiB、仅余240 kB。未重训、未修改旧分析或正式稿。',
'已更新PROJECT_STATUS.md、CODEX_HANDOFF.md、DECISIONS.md、RESULTS_MANIFEST.tsv，原记录已备份。结果停止状态与数值验收状态分开登记。']
text='\n'.join(lines)+'\n';name='26210_Phase315_文献复现与停止报告'
(P/'delivery'/f'{name}.txt').write_text(text,encoding='utf-8')
(P/'delivery'/f'{name}.html').write_text('<!doctype html><html lang="zh"><meta charset="utf-8"><title>Phase315 文献复现结果</title><style>body{max-width:1050px;margin:36px auto;font:16px/1.65 system-ui;color:#17212b;padding:0 20px}pre{white-space:pre-wrap;overflow-wrap:anywhere;font:inherit}a{color:#1265a8}</style><body><pre>'+html.escape(text)+'</pre></body></html>',encoding='utf-8')
(P/'resources/SOFTWARE.json').write_text(json.dumps({'python':platform.python_version(),'pandas':pd.__version__,'numpy':numpy.__version__,'scipy':scipy.__version__,'openpyxl':openpyxl.__version__,'environment':'batch_bench'},indent=2))
(P/'qc/AUTHOR_FIGURE_VISUAL_REVIEW.json').write_text(json.dumps({'status':'REVIEWED','source':'resources/author_Fig6.png','not_new_generated_figure':True,'P6_display':.58,'P7_display':.49,'display_statistic':'Rs','caption_statistic':'Pearson','review':'Original figure inspected; joint regions and published coefficient labels checked'},indent=2))
state={'status':'PARTIAL_REPRODUCTION_STOPPED','time':now,'scientific_full_reproduction_success':False,'source_directional_findings':'4/4 FAP-SPP1 positive; 40/40 GO raw P replicated; 3 available patients lower malignant CD3D/CD8A','stop_reasons':['P7 Pearson 0.439746 vs published 0.49; Spearman 0.485071 explains plotted rounding, caption conflicts','P7 immune source absent: candidate is exact 4457-row P6 duplicate, 0 coordinates match P7'],'numerical_validation':v,'external_validation_started':False,'TYMP_positioning_started':False,'analysis_jobs':'none','report':str(P/'delivery'/f'{name}.txt')}
(P/'STATUS.json').write_text(json.dumps(state,indent=2,ensure_ascii=False))
for file in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','AGENTS.md','RESULTS_MANIFEST.tsv']:
 if not (P/'root_record_backups'/file).exists():shutil.copy2(R/file,P/'root_record_backups'/file)
note='## Phase315终态（2026-10-04）：原文关系部分复现，未过完整复现门槛，按用户要求停止\n用户要求先复现，成功后才考虑外部验证定位，失败则停止报告。首篇固定Qi et al. Nat Commun 2022 13:1742，未改为其他论文。作者逐spot源数据四患者FAP/SPP1联合区均正：Pearson r=0.580833/0.439746/0.559819/0.553774；40项GO原始P准确复算，三主程序4/4富集；可用P6/P8/P9肿瘤区CD3D/CD8A均值及检出均更低。\n完整复现未通过：P7图值0.49对应Spearman0.485071而图注Pearson0.439746；P7候选免疫表实际与P6全4457行10列完全重复，与P7的4248坐标零匹配。状态PARTIAL_REPRODUCTION_STOPPED，属于方法标注/输入完整性问题，不是主要方向反向。266项独立数值核验PASS，分析/终验exit0；不将该PASS当科学复现成功。\n目录phase315_CRC_SPP1_barrier_paper_reproduction_20261004_v1，协议/全部结果/原文源表/原始失败与哈希保留。无外部验证、无TYMP定位、未联系作者、未改旧对象/正式稿。原始ST受控未申请，本轮仅源数据层复算，非端到端重建。交付delivery/26210_Phase315_文献复现与停止报告.txt与结果ZIP。当前任务按停止条件收束，不自动换研究或阈值继续。\n\n'
for file in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:
 f=R/file;f.write_text(note+f.read_text())
f=R/'AGENTS.md';f.write_text(f.read_text()+'\n\n## 2026-10-04 文献复现先行与失败停止\n用户明确：先复现SPP1下游文献关系，成功后再考虑外部验证定位；失败则停止并报告。Phase315为原论文源数据部分复现，因方法标注和P7数据缺口未过完整门槛，已停止。不得自动扩大外部验证、TYMP定位、换论文或调条件继续；后续由用户指令或明确补齐输入再开启。\n')
with (R/'RESULTS_MANIFEST.tsv').open('a',newline='') as f:
 wr=csv.writer(f,delimiter='\t',lineterminator='\n')
 for rid,kind,file in [('phase315_protocol','protocol','resources/REPRODUCTION_PROTOCOL_FROZEN.json'),('phase315_reproduced_correlations','partial_reproduction','tables/REPRODUCED_CORRELATIONS.tsv'),('phase315_GO_reproduction','reproduction','tables/GO_STATISTICS_REPRODUCTION.tsv'),('phase315_immune_maps','partial_reproduction','tables/IMMUNE_MAP_QUANTIFICATION.tsv'),('phase315_decision','stop_decision','qc/REPRODUCTION_DECISION.json'),('phase315_report','report',f'delivery/{name}.txt')]:
  wr.writerow([rid,'Phase315',kind,str(P/file),now,'PARTIAL_REPRODUCTION_STOPPED;NUMERICAL_QC_PASS',str(P/'qc/FINAL_VALIDATION.json'),'266 checks; no external validation/TYMP positioning; source-level only'])
files=[f for f in P.rglob('*') if f.is_file() and 'root_record_backups' not in f.parts and f.suffix!='.zip' and f.name not in ['DELIVERY_VALIDATION.json','PACKAGE_FILES.json']]
manifest={str(f.relative_to(P)):{'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in files}
(P/'delivery/PACKAGE_FILES.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False));files.append(P/'delivery/PACKAGE_FILES.json')
archive=P/'delivery/26210_Phase315_reproduction_results.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in files:z.write(f,str(f.relative_to(P)))
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for n,vv in manifest.items():assert hashlib.sha256(z.read(n)).hexdigest()==vv['sha256']
out={'status':'PASS','zip':str(archive),'bytes':archive.stat().st_size,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'members':len(files),'all_member_SHA_CRC':'PASS','scientific_status':'PARTIAL_REPRODUCTION_STOPPED','report_sha256':hashlib.sha256((P/'delivery'/f'{name}.txt').read_bytes()).hexdigest()}
(P/'qc/DELIVERY_VALIDATION.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2,ensure_ascii=False))
