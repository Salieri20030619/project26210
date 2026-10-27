from pathlib import Path
import json,csv,hashlib,datetime,shutil
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');R=P/'part2_rebuild_v1';C=R/(R/'provenance/CURRENT_CLUSTERING_RUN.txt').read_text().strip();M=R/'04_markers_annotation/run_v1';V=R/'05_patient_study_validation/run_v1';O=R/'11_reports/run_v1';O.mkdir()
validation=json.loads((R/'provenance/three_step_validation_v1/STATUS.json').read_text());assert validation['status']=='PASS'
visual=json.loads((R/'10_figures/run_v1/VISUAL_REVIEW.json').read_text());assert visual['status']=='PASS'
cl=json.loads((C/'STATUS.json').read_text());ann=pd.read_csv(M/'CLUSTER_ANNOTATIONS.tsv',sep='\t',dtype={'cluster':str});comp=pd.read_csv(V/'CLUSTER_COMPOSITION_QC.tsv',sep='\t',dtype={'cluster':str});merged=ann.merge(comp,on='cluster',validate='one_to_one')
rows=['# 第二部分重构 前三步结果报告','', '日期：2026-09-14。独立分支part2_rebuild_v1。旧第二部分对象、图表和报告均保留。','', '## 主要发现','', (M/'ANNOTATION_DECISION_CN.md').read_text(), '', '## 输入核查确认主分析范围','', '输入为当前75文库、五研究、40患者重建对象的42,546个高置信巨噬细胞。8,215个单核细胞作为独立边界参照，不参与巨噬精细分群。全部主分析细胞在既有QC中保留，无当前预测doublet或ambiguous注释。245个混合谱系复核标记保留，并完成去除这些细胞后重建图的敏感性检查。该标记用于复核，不作为已证实污染或双细胞的标签。', '', '## 无监督聚类及选择依据','', f"从原始counts重新选择研究分层的3,000高变基因，使用全转录组总量归一化、log1p CP10K、缩放PCA、按study的Harmony校正及30维kNN图；旧SPP1/C1QC或七亚型标签不参与。Leiden选择分辨率{cl['selected_resolution']}，得到{cl['clusters']}个簇。分辨率先按三种种子与80%分层抽样重建图的稳定性选择，再读取marker，不以出现SPP1簇为目标。", '', '初轮k30、0.2–1.0分辨率均未通过稳定性门槛，保留run_v1。随后在查看marker前固定k50与0.2、0.15、0.1、0.05的结构恢复计划，按相同门槛选择。正式分辨率与图敏感性表见provenance/CURRENT_CLUSTERING_RUN.txt所指目录内的RESOLUTION_STABILITY.tsv和GRAPH_SENSITIVITY.tsv。Harmony校正和UMAP展示不等于完全消除研究效应，簇的研究构成单独报告。', '', '## 分群后鉴定概览','', '|簇|注释|细胞数|患者级marker有效n|最大研究及比例|跨研究构成门槛|', '|---|---|---:|---:|---|---|']
for _,r in merged.iterrows():rows.append(f"|C{r.cluster}|{r.short_label}|{r.n_cells}|{r.patients_paired_marker_test}|{r.dominant_study} {100*r.dominant_study_fraction:.1f}%|{'通过' if r.broad_study_distribution_gate else '未通过'}|")
rows+=['','构成门槛是预定义的报告标记：至少3个研究各有≥20细胞及≥2患者、总体≥5患者、最大研究占比≤80%。它不是独立验证，也没有用来反向调整聚类。','', '## Marker统计与完整性','', '保留每簇全部37,795基因的表达、检出率、患者效应及可检验状态，不只展示SPP1。推断单位为患者；同一患者的簇内与其余巨噬细胞分别聚合原始counts，均≥20细胞才构成有效配对。比较log2(CPM+1)的配对差；预设表达门槛后采用双侧Wilcoxon（显式渐近法、连续性校正），在各簇对其余巨噬细胞的基因假设家族内做BH。n<3不作推断，n=3–9按探索性呈现。','', '聚类和marker发现使用同一表达数据，患者级FDR属于条件性marker发现证据，不是外部独立验证或因果证据。完整患者差值和伪bulk矩阵保存于04_markers_annotation/run_v1，研究内方向表位于05_patient_study_validation/run_v1。','', '## 验证与可复现性','', f"独立数值核查{validation['passed']}/{validation['checks']}项通过，涵盖输入身份、counts、图对称性、分辨率选择、完整基因家族BH、原始counts重新聚合的配对效应等；图件逐张视觉核查通过。程序、参数、来源及版本均保留。", '', '方法接口参考：[Scanpy高变基因](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html)、[Harmony接口](https://scanpy.readthedocs.io/en/stable/generated/scanpy.external.pp.harmony_integrate.html)、[Leiden接口](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.leiden.html)。实际运行软件版本见01_input_audit/run_v1/SOFTWARE_VERSIONS.json。','', '## 当前完成范围与后续衔接','', '本轮完成用户授权的三个步骤：输入核查、无监督精细聚类、分群后鉴定。轨迹、regulon、通信及空间模型没有在本轮重跑。后续应使用本分支实际确定的簇/状态身份逐项评估依赖，不把旧规则状态的结果直接改名成新簇的结果。','']
(O/'THREE_STEP_REPORT_CN.md').write_text('\n'.join(rows))
summary={'status':'THREE_STEPS_COMPLETE','completed_at':datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),'step1':'PASS_INPUT_SCOPE_AND_LINEAGE_REVIEW','step2':'PASS_UNSUPERVISED_CLUSTERING_REVIEW','step3':'PASS_POSTCLUSTER_IDENTIFICATION_WITH_REPORTED_BOUNDARIES','macrophages':42546,'patients':40,'libraries':75,'clusters':cl['clusters'],'selected_resolution':cl['selected_resolution'],'numerical_checks':validation['checks'],'visual_review':'PASS','new_trajectory_or_communication_models':False,'report':str(O/'THREE_STEP_REPORT_CN.md')}
(R/'STATUS.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2));(O/'COMPLETION.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
# Keep per-cell/patient files on server. Register primary review results only here.
artifacts=[O/'THREE_STEP_REPORT_CN.md',R/'STATUS.json',C/'STATUS.json',C/'RESOLUTION_STABILITY.tsv',C/'GRAPH_SENSITIVITY.tsv',M/'CLUSTER_ANNOTATIONS.tsv',M/'ANNOTATION_DECISION_CN.md',M/'TOP40_POSITIVE_PATIENT_MARKERS.tsv',V/'CLUSTER_COMPOSITION_QC.tsv',V/'TOP_AND_CANONICAL_MARKER_STUDY_EFFECTS.tsv',R/'provenance/three_step_validation_v1/STATUS.json',R/'10_figures/run_v1/VISUAL_REVIEW.json']
manifest=[]
for f in artifacts:manifest.append({'path':str(f),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
pd.DataFrame(manifest).to_csv(O/'DELIVERABLE_SHA256.tsv',sep='\t',index=False)
block=f'''## 2026-09-14 Phase147 第二部分独立重构前三步完成

按用户指定顺序完成输入核查、无监督精细聚类和分群后鉴定。新结果仅在part2_rebuild_v1。42,546巨噬细胞、40患者、75文库；另8,215单核作边界参照。选定Leiden分辨率{cl['selected_resolution']}，共{cl['clusters']}簇；参数选择先于marker鉴定，未使用旧SPP1标签定义聚类。完整注释与证据定位见04_markers_annotation/run_v1/ANNOTATION_DECISION_CN.md；报告11_reports/run_v1/THREE_STEP_REPORT_CN.md。数值{validation['passed']}/{validation['checks']}核查及图件视觉核查通过。当前结果是条件性发现，保留研究限制/谱系复核标记；未自动重跑轨迹、regulon或通信。旧第二部分完全保留。当前断点：用户审阅新分群和SPP1证据，再决定下游衔接。

'''
for n in ['CODEX_HANDOFF.md','PROJECT_STATUS.md']:
 f=P/n;f.write_text(block+f.read_text())
with (P/'DECISIONS.md').open('a') as h:h.write('\n'+block)
f=P/'RESULTS_MANIFEST.tsv';header=next(csv.reader(f.open(),delimiter='\t'))
with f.open('a') as h:
 w=csv.DictWriter(h,fieldnames=header,delimiter='\t',lineterminator='\n')
 for i,r in enumerate(manifest,1):w.writerow(dict(result_id=f'phase147_complete_{i:02d}',stage='Phase147',artifact_type=Path(r['path']).suffix.lstrip('.'),absolute_path=r['path'],generated_at=summary['completed_at'],validation_status='PASS_THREE_STEP_DISCOVERY_WITH_REPORTED_BOUNDARIES',validation_evidence=str(R/'provenance/three_step_validation_v1/STATUS.json'),notes='Independent new Part2. SHA256='+r['sha256']))
print(json.dumps(summary,ensure_ascii=False))
