from pathlib import Path
import json,hashlib,datetime,shutil
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');R=P/'part2_rebuild_v1'
C=R/'03_macrophage_clustering/run_v4_exploratory';M=R/'04_markers_annotation/run_v4_macrophage_only';V=R/'05_patient_study_validation/run_v4_macrophage_only';F=R/'10_figures/run_v3'
O=R/'11_reports/run_v1';O.mkdir(exist_ok=False)
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def writej(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')
a=pd.read_csv(C/'CELL_CLUSTER_ASSIGNMENTS.tsv.gz',sep='\t',dtype={'cluster':str},low_memory=False)
c4=a[a.cluster=='4'];assert len(c4)==2674 and c4.curated_patient_id.nunique()==37
cols=['cell_id','cluster','study','curated_patient_id','sample_id','source_tissue']
D=R/'09_downstream_crosswalk';c4[cols].to_csv(D/'C4_CANDIDATE_CELL_IDS.tsv.gz',sep='\t',index=False)
st=json.loads((C/'STATUS.json').read_text());ev=pd.read_csv(M/'C4_KEY_MARKER_EVIDENCE.tsv',sep='\t');spp=ev[ev.gene=='SPP1'].iloc[0]
ann=pd.read_csv(R/'04_markers_annotation/run_v2/CLUSTER_ANNOTATIONS.tsv',sep='\t',dtype={'cluster':str})
core=pd.read_csv(C/'CLUSTER_CORE_STABILITY.tsv',sep='\t',dtype={'cluster':str})
summary=ann.merge(core,on='cluster');summary['cells']=summary.cluster.map(a.cluster.value_counts());summary.to_csv(O/'ALL_CANDIDATES_AUDIT.tsv',sep='\t',index=False)
def mdtable(df):
 return '| '+' | '.join(map(str,df.columns))+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(str(v).replace('|','/').replace('\n',' ') for v in row)+' |' for row in df.itertuples(index=False,name=None))
state={'state_id':'P2R1_C04','name':'GPNMB/CCL18-high, SPP1-enriched macrophage candidate','status':'EXPLORATORY_DISCOVERY_WITH_SUPPORTED_ALGORITHMIC_CORE','definition':'Exact frozen unsupervised cluster 4 membership; no SPP1 threshold used for selection.','n_cells':2674,'patients_present':37,'paired_patients':26,'studies':5,'global_partition_stable':False,'minimum_seed_Jaccard':0.8573877271169009,'subsample_Jaccard':0.7589121646453509,'source_graph':str(C/'CLUSTERING_GRAPH.h5ad'),'source_graph_sha256':sha(C/'CLUSTERING_GRAPH.h5ad'),'assignments_sha256':sha(C/'CELL_CLUSTER_ASSIGNMENTS.tsv.gz'),'cell_ids_sha256':sha(D/'C4_CANDIDATE_CELL_IDS.tsv.gz'),'primary_markers':ev.to_dict(orient='records'),'scope':'Discovery dataset conditional marker contrasts; not external validation, lineage transition, mechanism, or clinical utility.','generated_at':now}
writej(D/'C4_STATE_DEFINITION.json',state)
report='''# 第二部分无监督发现三步执行与审计报告

## 无监督分群发现SPP1富集巨噬候选C4，多基因及患者方向支持其身份

本次完成输入与谱系QC、无监督精细聚类、分群后完整marker鉴定及患者/研究一致性核查。全部结果位于独立part2_rebuild_v1，未覆盖既有第二部分。C4包含2,674个细胞，出现在37名患者和五项研究；在26名符合巨噬细胞内配对比较条件的患者中，SPP1有24名正向、2名反向，患者配对log2(CPM+1)差中位数1.448645，P=6.677814821616047e-05，簇内基因家族BH FDR=0.000544714373533。GPNMB、CCL18、CYP27A1、LPL、APOE共同支持其多基因身份，可表述为“GPNMB/CCL18高表达、SPP1富集的巨噬细胞候选状态”。

完整17簇划分未达到预设全局稳定性门槛；当前数据不支持把全部17簇称为已确认的稳定生物亚型。C4逐簇核心稳定性达标，支持对该候选开展有边界的发现性解释。技术检查通过不等于独立生物学验证。

## 第一阶段确认75文库来源，并逐细胞保留纯化依据

来源对象：/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad。
SHA256：e4779c46af13a8643e8dec584eeff4322258181f5edab3e55dcd8c20bc28a9e7。
原对象258,120细胞、75文库、40名整理后患者、五研究，37,795个原始计数基因。五研究为GSE156632、GSE159115、GSE178481、GSE207493、GSE242299。

抽取广义Macrophage 42,546细胞，以及Monocyte 8,215细胞作边界背景，共50,761髓系细胞。检查细胞ID、来源、计数行和、有限非负整数和患者映射。已有918个Epithelial/Myeloid混合细胞未进入此候选集。初始245个混合谱系复核标记保留进入首轮审查，不按SPP1删除或保留。

首轮聚类提示低复杂度、血管周混合信号。随后设定精细分析至少500个检出基因、1000个总原始counts；排除首轮血管/周细胞混合C2；非髓系面板检出至少3基因且平均logCP10K超过髓系骨架者排除。非髓系面板包括上皮、T、NK、B、内皮、周细胞、成纤维及肥大细胞；髓系骨架为LST1/TYROBP/FCER1G/AIF1/CTSS/LYZ/CSF1R/CD68/CTSB。该步骤不以SPP1表达或旧SPP1/C1QC标签选择细胞。

排除并集12,872细胞：低复杂度标记11,819、首轮血管混合866、非髓系主导453，三类存在重叠。保留29,674候选细胞，仍覆盖40患者、75文库和五研究。逐细胞原因在02_myeloid_qc/run_v2/CELL_PURIFICATION_AUDIT.tsv.gz，保留集在PURIFIED_CELL_METADATA.tsv.gz。总counts不得跨不同测序协议一律写成UMI；早期记录中“UMI”列名/文字在本报告统一按原始counts解释。

## 第二阶段建立探索性候选目录，同时保留全局稳定性未通过的结果

首轮run_v1尝试的分辨率未通过稳定性门槛。run_v2低分辨率得到4簇，算法一致性很高，但事后完整标记表明其主要结构包含27,395细胞核心、13,582低复杂度细胞（中位481 counts/278基因，87.9%来自GSE242299）、866血管周混合细胞及703增殖细胞，不能作为4个巨噬生物亚型。该结果用于QC回溯并完整保留。

纯化后重新计算HVG/PCA/Harmony和邻接图。3000个study-aware seurat_v3高变基因，候选基因至少20细胞表达，聚类特征排除线粒体/核糖体/Hb；完整基因表达仍保留用于标记鉴定。使用全库计数CP10K、log1p、标准化PCA50；Harmony按study整合（theta=2、nclust=50、maxiter=20），使用30维、50邻居；Leiden比较0.8/0.6/0.4/0.2/0.15/0.1/0.05及种子262147、42、31415，另用80%按研究分层抽样重建图。

预设全局门槛：最低种子对ARI>=0.85、重抽样ARI>=0.80、小于100细胞簇占比<=1%。run_v3没有非平凡完整划分通过；脚本按设计报STABILITY_UNRESOLVED并保存结果。在新一轮marker查看前，按已有网格最高最低种子ARI选择0.8进行探索性鉴定；17簇最低种子ARI=0.80060316056，重抽样ARI=0.6328556116922294。

run_v4_exploratory评估各候选核心与其他种子/重抽样最佳匹配的Jaccard，阈值分别0.75和0.70。9/17核心通过：C2、C4、C5、C6、C10、C11、C13、C14、C15。核心算法稳定与细胞身份是不同证据，混合谱系簇也能算法稳定。C4两项为0.857388和0.758912。最终使用run_v3冻结connectivity进行聚类与UMAP；显示用重建PCA/Harmony存在0.0043447相对图差异，因此恢复冻结图而未悄然更换分析邻接关系。详见GRAPH_PROVENANCE.md。

## 第三阶段先鉴定全部候选，再在巨噬候选内复核SPP1

完整17候选注释及稳定性如下，未将DC/单核边界或混合谱系群删除出审计目录。

{annotation_table}

成熟巨噬候选比较排除C1、C2、C13、C14、C16，保留12组、22,444细胞、39患者、70有效文库。70文库来自已审计75文库的细胞子集，不是更换输入范围。完整17组结果在04_markers_annotation/run_v2；最终巨噬内统计以04_markers_annotation/run_v4_macrophage_only为准，不能混用两种对照范围的P/FDR。

C4巨噬骨架TYROBP检出100%、CD68约98.9%、C1QC约98.4%；GPNMB检出97.6%，SPP1检出47.8%。因此该候选不能写成“所有细胞SPP1阳性”或“C1QC阴性”。其最大研究贡献为GSE207493约51.0%，最大患者约10.2%，覆盖五研究。原始SPP1信号更高的C14为上皮/髓系混合（CD24/KRT18/CRYAB/PDZK1IP1等），已排除出巨噬比较，不能仅凭SPP1命名为巨噬亚型。

最终C4标记证据：

{evidence_table}

C1QC仍普遍表达，但相对其他巨噬候选的患者配对差中位数为-0.107770，20/26患者负向，P=0.0147604、FDR=0.0322570。这一结果不支持把C4与C1QC表达定义为互斥谱系。

按研究描述SPP1患者配对方向：GSE156632为4/4正向、中位差1.192548；GSE159115为3/3、2.229312；GSE178481为2/2、4.790183；GSE207493为10/12、1.097950；GSE242299为5/5、2.713788。五项研究中位方向均正。各研究样本量小，尤其2例研究只报告方向与效应；不宣称五项研究分别达统计显著或属于独立外部验证。完整表见05_patient_study_validation/run_v4_macrophage_only/C4_KEY_MARKER_STUDY_EFFECTS.tsv。

C5/C6等研究集中、C7/C11应激、C10身份未解且仅2名有效配对患者、C9研究偏倚以及C12边界/低复杂度信息均保留。完整候选身份、组成、全基因P/FDR和配对效应见关联目录；不能把这些次级候选统称为跨研究确认亚型。

## 统计方法以患者为单位，FDR在每个簇的同类基因家族内计算

按患者聚合原始基因counts，形成簇内与其他纳入候选的配对pseudobulk。双方均至少20细胞方可纳入该簇比较；效应为log2(CPM+1)配对差及其中位数。基因过滤要求至少max(3,ceil(0.2*n))名患者在配对任一聚合侧达到10 counts。对符合表达条件且配对差非全零的基因使用双侧Wilcoxon近似检验并作连续性校正，BH在每个cluster-versus-rest同类基因家族内计算，完整37,795基因表保留测试资格与缺失值。不使用细胞级P值。少于3名有效配对患者不作推断性检验，3–9名按探索性口径，达到10名报告效应/P/FDR。图中每点为一名患者、粗横线为中位数。

实现环境：anndata0.12.16、scanpy1.11.5、numpy2.4.6、scipy1.17.1、pandas2.3.3、harmonypy2.0、igraph1.0、leidenalg0.11、skmisc0.5.2。实现细节与参数以版本化scripts和logs为准。
官方方法接口：[HVG](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html)、[Harmony](https://scanpy.readthedocs.io/en/stable/generated/scanpy.external.pp.harmony_integrate.html)、[Leiden](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.leiden.html)。

## 讨论：已补齐SPP1候选发现链，整体亚型体系及下游机制仍有明确边界

当前已建立“来源/QC→不按SPP1筛选的无监督聚类→完整marker与谱系审查→巨噬候选内患者比较→研究方向与候选核心稳定性”的可复核发现链。SPP1是事后识别的多基因状态特征，不是定义初始SPP1-high细胞后再验证同一阈值。SPP1本身未被刻意从全转录组特征中屏蔽；无监督不等于基因盲法。

发现与marker对比使用同一数据，P/FDR属于给定数据驱动分群后的条件性发现证据，未消除聚类选择的影响。研究分层为内部一致性检查。算法核心重复性不等同于留研究重新聚类或独立队列复现。复杂度阈值和谱系纯化影响细胞构成，完整排除审计保留；全局不稳定与若干研究偏倚限制完整分类体系的外推，但不抹去C4的已测量富集方向与效应。未建立细胞转化方向、因果机制或临床可用性。

本次完成授权三步，不代表第二部分全部轨迹、regulon、通信分析已重建。未自动运行这些下游模型，未改写第三/四部分，未接管客户湿实验。C4细胞ID和状态定义已冻结供后续依赖评估使用。

## 复核、修正和交付

独立数值核验provenance/three_step_validation_v2为385/385通过；provenance/macrophage_only_validation_v1为278/278通过，共663项技术/数值检查，包含ID、原始计数再聚合、患者配对、BH家族和图一致性等，不代表663项独立生物学验证。

以下历史失败不覆盖：初始分辨率网格失败、纯化后完整划分稳定性失败；macrophage_only run_v3因沿用40患者断言在统计前中止，run_v4按实际39患者/70文库修正；首次验证的2项float32均值精度问题通过float64独立重算解决，未放宽容差；首次绘图读取object dtype基因索引失败，改用明文基因索引；run_v2图05漏画已计算的重抽样点，run_v3补齐0.632856且保留前版。

最终七图位于10_figures/run_v3，逐图可视核查通过；PNG/SVG共14文件与本地交付SHA256一致，见VISUAL_REVIEW.json和provenance/FIGURE_DELIVERY_RECEIPT.json。本地目录：C:/Users/31937/Downloads/26210/15_part2_rebuild_discovery_20260914_v1。完整审计及进度留在服务器。

主要入口：
- 03_macrophage_clustering/run_v4_exploratory/CLUSTERING_GRAPH.h5ad（冻结图与候选标签，表达取01_input_audit的原始对象）
- 04_markers_annotation/run_v2/CLUSTER_ANNOTATIONS.tsv（17候选完整身份）
- 04_markers_annotation/run_v4_macrophage_only/（12组巨噬内最终统计及全基因表）
- 05_patient_study_validation/run_v4_macrophage_only/（患者/研究证据）
- 09_downstream_crosswalk/C4_CANDIDATE_CELL_IDS.tsv.gz及C4_STATE_DEFINITION.json
- 11_reports/run_v1/ALL_CANDIDATES_AUDIT.tsv（完整注释与核心稳定性合表）
'''
evcols=['gene','eligible_patients','median_patient_log2CPM_difference','n_patients_positive','n_patients_negative','patient_wilcoxon_p','BH_FDR_cluster_gene_family']
report=report.replace('{annotation_table}',mdtable(summary[['cluster','cells','short_label','annotation_cn','minimum_seed_Jaccard','subsample_Jaccard','cluster_core_reporting_gate','interpretation']])).replace('{evidence_table}',mdtable(ev[evcols]))
(O/'THREE_STEP_REPORT_CN.md').write_text(report)
caption='''第二部分重构三步发现图件（最终远程10_figures/run_v3）
01：29,674候选细胞的冻结无监督UMAP、研究来源与谱系复核标记；17候选为探索性目录。
02：各候选细胞数及研究构成。细胞数不等于独立统计样本量。
03：全部17候选的完整谱系/状态标记展示；颜色为基因跨簇标准化平均log1p CP10K，点大小为检出比例。
04：左图为全部候选SPP1表达（包括上皮/髓系混合C14）；右图为仅12组巨噬候选之间的患者配对SPP1差，比较目录共39患者，各簇有效n见横轴。
05：全局最低种子ARI及选定0.8分辨率的80%重抽样ARI；全局门槛均未通过，绿色方点为0.632856。其他分辨率未执行/无有效重抽样值时不补造曲线。
06：各候选核心最佳匹配Jaccard；9/17通过0.75种子与0.70重抽样门槛。算法稳定性不等于巨噬身份或外部复现。
07：C4六基因在26患者内相对其他巨噬候选的配对log2(CPM+1)差。每点一患者、横线中位数，q为完整C4基因家族BH FDR，并非只对图示六基因校正。SPP1为24/26正向，效应中位1.448645，P=6.6778e-05，FDR=0.000544714。
共同口径：同一数据驱动分群后的条件性发现证据；五研究方向核查为内部一致性。C4包含2,674细胞、37名出现患者，26名满足双方各至少20细胞的配对条件。全候选为40患者/75文库；12组巨噬比较子集为22,444细胞/39患者/70文库。旧第二部分未覆盖。
'''
(F/'图件说明.txt').write_text(caption)
status={'status':'COMPLETE_THREE_STEPS_EXPLORATORY_DISCOVERY','generated_at':now,'completed_steps':['input_lineage_QC','unsupervised_fine_clustering','postcluster_marker_patient_study_review'],'source_cells':258120,'source_patients':40,'source_libraries':75,'purified_candidate_cells':29674,'candidate_clusters':17,'algorithmically_supported_cores':9,'global_partition_stability_pass':False,'macrophage_comparison_cells':22444,'macrophage_comparison_patients':39,'macrophage_comparison_libraries':70,'macrophage_candidate_groups':12,'C4':state,'numeric_checks_passed':663,'figures_reviewed':7,'figure_run':'10_figures/run_v3','report':'11_reports/run_v1/THREE_STEP_REPORT_CN.md','downstream_models_run':False,'old_part2_overwritten':False}
old=R/'STATUS.json';shutil.copy2(old,O/'BRANCH_STATUS_BEFORE_FINALIZATION.json');writej(old,status)
for name,value in {'CURRENT_CLUSTERING_RUN.txt':'03_macrophage_clustering/run_v4_exploratory','CURRENT_MARKER_RUN.txt':'04_markers_annotation/run_v2','CURRENT_VALIDATION_RUN.txt':'05_patient_study_validation/run_v2','CURRENT_PRIMARY_MACROPHAGE_MARKER_RUN.txt':'04_markers_annotation/run_v4_macrophage_only','CURRENT_PRIMARY_MACROPHAGE_VALIDATION_RUN.txt':'05_patient_study_validation/run_v4_macrophage_only','CURRENT_FIGURE_RUN.txt':'10_figures/run_v3','CURRENT_REPORT.txt':'11_reports/run_v1/THREE_STEP_REPORT_CN.md'}.items():(R/'provenance'/name).write_text(value+'\n')
ag=R/'AGENTS.md';t=ag.read_text();t=t.replace('- 当前仅完成目录建立，未运行聚类。后续先形成输入清单与可审阅的分析方案，再执行相应阶段；不为一般可逆分析额外引入审批门槛。','- 当前已完成用户授权的输入/QC、无监督精细聚类、分群后鉴定三步；最新结果入口见STATUS.json。全局分类仍探索性。一般可逆分析不额外引入审批门槛。');ag.write_text(t)
read=R/'README.md';t=read.read_text().replace('当前状态：DIRECTORY_READY，尚未启动分析。','当前状态：COMPLETE_THREE_STEPS_EXPLORATORY_DISCOVERY。').replace('本次仅建立目录、规则和状态记录。聚类参数、最终入选细胞和统计家族尚未冻结；空阶段目录不代表已完成分析。','输入/QC、无监督精细聚类及分群后鉴定三步已完成。C4为GPNMB/CCL18高、SPP1富集候选；17簇整体不稳定。完整报告见11_reports/run_v1/THREE_STEP_REPORT_CN.md；下游空目录不代表完成。').replace('实际巨噬/单核入选范围需在01_input_audit确认。','实际入选与纯化范围已记录在01_input_audit及02_myeloid_qc。');read.write_text(t)
entry='''## 2026-09-14 Phase147 三步完成：发现SPP1富集巨噬候选C4

独立part2_rebuild_v1完成输入/QC→无监督精细聚类→完整marker/患者研究核查。29,674候选、40患者、75文库；17候选全局ARI门槛未过，9核心算法支持，不称17稳定生物亚型。排除单核/DC及混合谱系后，最终巨噬内比较为12组、22,444细胞、39患者、70文库。C4为2,674细胞、37名出现患者；26名有效配对患者中SPP1有24名正向，配对log2(CPM+1)差中位1.448645，P=6.6778148e-05、BH FDR=0.0005447144，五研究中位方向同向。GPNMB/CCL18等多基因支持；C14上皮/髓系混合虽SPP1高，未当作巨噬亚群。C4核心Jaccard0.857/0.759，整体分类及同数据发现保持探索性边界。

最终巨噬内统计：04_markers_annotation/run_v4_macrophage_only及05_patient_study_validation/run_v4_macrophage_only；17候选完整注释仍在04_markers_annotation/run_v2。报告：part2_rebuild_v1/11_reports/run_v1/THREE_STEP_REPORT_CN.md。技术核验663/663通过，7图完成视觉检查，14个PNG/SVG本地交付哈希一致；最终图10_figures/run_v3。本地Downloads/26210/15_part2_rebuild_discovery_20260914_v1。C4细胞ID/定义已冻结。所有失败和排除审计保留，旧第二部分未覆盖，第三/四部分与湿实验未改动，未运行新轨迹/regulon/通信。

以下较早记录是执行历史，最新状态以上述完成记录及分支STATUS.json为准；早期UMI措辞实际指原始counts。

'''
for fn in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:
 p=P/fn;shutil.copy2(p,O/(fn+'.before_phase147_final'));p.write_text(entry+p.read_text())
artifacts=[O/'THREE_STEP_REPORT_CN.md',O/'ALL_CANDIDATES_AUDIT.tsv',D/'C4_CANDIDATE_CELL_IDS.tsv.gz',D/'C4_STATE_DEFINITION.json',M/'C4_KEY_MARKER_EVIDENCE.tsv',V/'C4_KEY_MARKER_STUDY_EFFECTS.tsv',R/'provenance/three_step_validation_v2/STATUS.json',R/'provenance/macrophage_only_validation_v1/STATUS.json',F/'VISUAL_REVIEW.json',R/'provenance/FIGURE_DELIVERY_RECEIPT.json',F/'图件说明.txt',R/'STATUS.json']
manifest=[]
for i,p in enumerate(artifacts):
 assert p.is_file();manifest.append({'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size})
writej(O/'FINAL_ARTIFACT_MANIFEST.json',{'generated_at':now,'artifacts':manifest})
with (P/'RESULTS_MANIFEST.tsv').open('a') as f:
 for i,p in enumerate(artifacts):f.write('\t'.join(['phase147_complete_'+str(i+1),'phase147','report_or_evidence',str(p),now,'COMPLETE_EXPLORATORY_DISCOVERY','663 technical checks; global partition gate false','New Part2 only; see final report'])+'\n')
print('FINALIZED',str(O/'THREE_STEP_REPORT_CN.md'),'C4',len(c4),'artifacts',len(manifest))
