from pathlib import Path
import json,hashlib,subprocess
import pandas as pd
w=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase127_section2_4_pyscenic_redraw_remote_v1');q=w/'qc'
note="2.4 恶性上皮 cNMF–pySCENIC 图件（v1）\n\n本包包含7张图，每张提供360 dpi PNG和PDF，共14个图件及本说明。建议依次查看03热图、04连续关联、05小样本分组效应、06环形网络；01/02为两种算法的完整F1–F7特异性排名，07为候选调控子的UMAP展示。\n\n结果\nILF2、HOXD3和HNF4A调控子在两种算法中均与F1呈正向患者层关联（20位患者，各算法全调控子检验家族BH FDR<0.05）。F1主导组比较包括9位配对患者，作为探索性结果报告方向和效应量。历史4终点规则下，无候选达到4/4；ILF2和HOXD3的额外分组证据仅为名义显著，不能称为分组FDR显著或已经确认的因果驱动。\nHIF1A调控子在两套正式AUCell矩阵中均未检出；此前8项预设靶向敏感性运行也均未恢复HIF1A，因此本包不生成HIF1A调控活性图。\n\n图件说明\n01/02 A：分别展示GENIE3与GRNBoost2的F1–F7调控子特异性排名。为确保两种算法展示同一候选空间，均使用197个共有调控子。标注每组RSS前三位；这与根据患者层连续关联确定的候选TF是不同的排序问题。\n03 B：24个调控子的患者平衡AUCell热图。每模块依两算法RSS平均名次选前三位，再加入预先确定的ILF2、HOXD3、HNF4A。每患者×模块至少10个细胞，先取该患者模块内AUCell中位数，再跨患者取中位数；各算法分别在F1–F7间按行z标准化。色阶[-2,2]，超出范围仅颜色饱和；完整值在服务器表中。左右相邻列是同一模块的两种算法，不是治疗或响应分组。树为显示矩阵的欧氏距离average linkage，不定义新的生物亚型。\n04 C：20位患者的F1–AUCell Spearman相关。每点为一个患者，黑菱形为中位数，标注效应量、双侧Wilcoxon P和原检验家族BH FDR。\n05 D：F1主导细胞与其他主导模块细胞的患者内AUCell中位数差。两组每患者均至少10个细胞，共9位配对患者。每点为一个患者，黑菱形为中位数；正文优先展示方向一致性和效应量，不添加显著性星号。两种算法使用各自的AUCell数值尺度。\n06 E：三个候选TF的环形网络。图顶端的GENIE3/GRNBoost2节点是同一TF在两种算法中的网络表示。每条曲线都来自正式TF–靶基因预测。节点颜色区分共同支持、仅GENIE3及仅GRNBoost2；每类按方法内权重百分位的平均值最多显示10个靶点。排除TF自身节点；不按TYMP是否出现来挑选靶基因。完整未截取网络保存在服务器。\n07 F：沿用全部11,950个恶性细胞的冻结UMAP，每种算法分别展示ILF2、HOXD3、HNF4A AUCell。每个面板用其实际完整范围，数值不平滑、不截尾。\n\n方法\n复用原11,950细胞、26患者、5队列的正式cNMF标签及AUCell矩阵，GENIE3为299个调控子，GRNBoost2为314个，两者共有197个。RSS采用已安装官方pySCENIC源码的1−Jensen–Shannon距离定义，并用KL-to-mixture等式独立复算。\n连续终点按患者内相关汇总，FDR分别在各算法全调控子家族内计算；分组终点按“算法×cNMF模块”分别计算FDR。原P/FDR保留，不把连续、分组、RSS或网络混入一个校正池。本次只新增描述性绘图汇总，不重跑GRN、motif剪枝、AUCell、cNMF或UMAP。\n热图各模块实际患者数：F1=9，F2=12，F3=2，F4=4，F5=4，F6=2，F7=4；共37份患者×模块谱。同患者可贡献多个模块谱，图中各列n不能相加作为独立患者数。\n\n讨论与解读范围\nAUCell衡量表达推断调控子基因集的富集，不是直接蛋白测量或TF结合证据。两种算法共享数据来源，不等于独立生物学验证。RSS及模块热图受队列、模块细胞比例与可用患者数影响；展示峰值不能替代患者层检验，也不证明调控F1的因果方向。\nILF2的TYMP预测边仅在GRNBoost2中出现，GENIE3没有该边；HOXD3和HNF4A的两套网络均没有TYMP边。当前数据不支持将三个候选均称为已经验证的TYMP直接上游调控因子。\nHIF1A未检出指本次网络无法评估其调控子，不等于HIF1A在生物学上没有作用。\n"
cs=pd.read_csv(q/'CANDIDATE_CONTINUOUS_DISPLAY.tsv',sep='\t')
gs=pd.read_csv(q/'CANDIDATE_GROUP_COMPLETE_STATS.tsv',sep='\t')
note+='\n候选连续终点完整统计\n'
for r in cs.itertuples():note+=f'{r.tf} | {r.method} | n={r.n_eligible_patients} | median rho={r.median_patient_effect_internal:.7g} | P={r.wilcoxon_p_two_sided:.7g} | FDR={r.BH_FDR:.7g} | 正/负/零={r.n_positive}/{r.n_negative}/{r.n_zero}\n'
note+='\n候选分组终点完整统计（n=9，探索性）\n'
for r in gs.itertuples():note+=f'{r.tf} | {r.method} | n={r.n_eligible_patients} | median delta={r.median_patient_effect_internal:.7g} | P={r.wilcoxon_p_two_sided:.7g} | FDR={r.BH_FDR:.7g} | 正/负/零={r.n_positive}/{r.n_negative}/{r.n_zero}\n'
pf=pd.read_csv(w/'data/PATIENT_FACTOR_AUC_MEDIANS.tsv',sep='\t')
note+=f'\n37份热图谱来自{pf.patient.nunique()}位不同患者。\n'
note+='\n复核与记录\n原始输入SHA256完全一致；156条候选患者记录的效应及资格复算通过，原家族BH复算通过；RSS两种计算最大误差2.23e-16。七张PDF完成结构与最终渲染视觉检查。\n唯一记录目录：'+str(w)+'\n所有脚本、逐细胞AUCell、全患者统计、完整网络、源哈希、软件版本及审核材料只保存在服务器；本地仅交付图件和本说明。\n'
(w/'figures/README_2.4_pySCENIC.txt').write_text(note,encoding='utf8')
# Final edge identity validation after layout refinement.
edges=pd.read_csv(q/'CIRCULAR_DISPLAY_EDGES.tsv',sep='\t')
targets=pd.read_csv(w/'data/CANDIDATE_REGULON_TARGETS.tsv',sep='\t')
lookup=targets.set_index(['tf','target','method']).weight
for r in edges.itertuples():assert abs(lookup.loc[(r.tf,r.target,r.method)]-r.weight)<1e-12
(q/'VISUAL_REVIEW.json').write_text(json.dumps({'status':'PASS','reviewed_pdf_renderings':7,'render_scale':1800,'final_refinements':['D decimal ticks and two-line xlabel','E expanded lower label space'],'all_labels_legends_and_networks':'readable; no clipping or overlapping labels in final reviewed renderings','user_review':'PENDING'},indent=2))
checks=[]
for f in sorted((w/'figures').glob('*.pdf')):
 info=subprocess.run(['pdfinfo',str(f)],capture_output=True,text=True,check=True).stdout
 assert any(l.startswith('Pages:') and l.split(':')[1].strip()=='1' for l in info.splitlines())
 checks.append({'name':f.name,'pages':1,'final_structure':'PASS'})
assert len(checks)==7
(q/'PDF_VALIDATION.json').write_text(json.dumps(checks,indent=2))
manifest=json.loads((q/'SOURCE_MANIFEST.json').read_text());known={x['path'] for x in manifest}
for f in list((w/'scripts').glob('*.py'))+[w.parent/'local_records_migration_20260908_150455/references/image12.png']:
 if str(f) not in known:manifest.append({'path':str(f),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
(q/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
rows=[{'name':f.name,'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in sorted((w/'figures').iterdir()) if f.is_file()]
assert len(rows)==15
(q/'DELIVERY_MANIFEST.json').write_text(json.dumps(rows,indent=2))
print(json.dumps({'files':len(rows),'bytes':sum(r['bytes'] for r in rows),'heatmap_unique_patients':pf.patient.nunique(),'network_edges_verified':len(edges)}))
