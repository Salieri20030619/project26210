from pathlib import Path
import json,csv,hashlib,datetime,shutil,re,gzip
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=P/'phase145_macrophage_provenance_trace_remote_v1';Q=Path('/mnt/BioDisk/ludingyu/project26210');R=Path('/mnt/BioDisk/wangmingyuan/project26210_part4_sensitivity_pilot_v1')
# Preserve exact small source files and hash the large RDS in place, with before/after file metadata.
sources=[Q/'object/02_macrophage_final_subtypes.rds',Q/'run_macro_subtype_cellchat.R',Q/'p4_02_build_reference.py',Q/'plot_macro_subtype_figs.R',Q/'run_liana_tam.R',Q/'logs/macro_subtype_cellchat.log',Q/'docs/client_reports/part3_cell_communication/第三部分_细胞通讯分析报告_含巨噬细分 (2).docx',R/'build_patient_balanced_reference_v2.py',R/'reference_v2/run_metadata_v2.json',R/'reference_v2/reference_v2_cell_assignments.tsv.gz',R/'build_cell2location_reference_v3.py']
rows=[]
for i,f in enumerate(sources):
 if not f.exists():continue
 st=f.stat();h=hashlib.sha256()
 with f.open('rb') as z:
  for b in iter(lambda:z.read(8*1024*1024),b''):h.update(b)
 after=f.stat();assert st.st_size==after.st_size and st.st_mtime_ns==after.st_mtime_ns
 rows.append({'source':str(f),'bytes':st.st_size,'sha256':h.hexdigest(),'mtime_utc':datetime.datetime.fromtimestamp(st.st_mtime,datetime.timezone.utc).isoformat(),'unchanged_during_hash':True})
 if st.st_size<2000000:shutil.copy2(f,O/'evidence'/('source_'+str(i)+'_'+f.name))
with (O/'SOURCE_FILE_MANIFEST.tsv').open('w') as h:
 w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
report='''# 巨噬细胞精细分群及SPP1来源详细追溯

日期：2026-09-14。项目26210，Phase145。用户说明此前这部分不是其本人完成，因此本轮依据服务器对象、代码、日志和报告恢复事实，不把文件目录所属账号当作分析作者。

## 结论

历史上确实做过巨噬细胞无监督聚类，并产生七亚型；不能再概括为“以前没有精细分群”。但旧分群高度集中于不同研究，其完整生成脚本、原始marker筛选表和删去4,204细胞的依据未能恢复。它目前不足以直接作为跨队列稳定SPP1亚群的发现证据。

当前主线使用另一套定义：2026年8月5日为第四部分空间参考敏感性分析，在新版QC对象上以患者内SPP1/C1QC程序划定状态。旧七亚型与当前状态不是同一标签体系，不存在已验证的直接继承关系。已有轨迹和通信主要描述后一套预定义状态，不能补上前一套分群的独立发现与验证环节。

因此，用户希望的“巨噬细胞精细分群→由完整marker识别SPP1亚群→功能/轨迹/通信”证据链，有历史分群对象基础，但仍未在当前重建数据中完整建立。

## 追溯范围与证据层级

本轮只读检查旧巨噬Seurat对象、旧注释H5AD元数据、当前重建H5AD元数据、reference-v2细胞标签及当前QC表；扫描旧项目1,984份脚本/文本文档（命中7个文件、25条引用），提取10份已有Word文档文本；检索17个以旧project26210为工作目录的历史会话（18个匹配事件），以及8月5日3个相关会话（124个匹配事件）。本地Downloads及D盘项目进行了相关文件名补查，发现旧第三部分报告副本，未发现独立的巨噬生成脚本。

Seurat对象证据与逐细胞匹配优先于报告叙述；报告仅作为历史声称，不能作为独立验证。保留全文命令记录、元数据、事件索引、原文件SHA256和分析脚本。本轮未重新聚类、未做差异基因推断检验、未改写任何原始分析对象或既有图件。

## 一 旧对象确实来自聚类后命名

对象：`/mnt/BioDisk/ludingyu/project26210/object/02_macrophage_final_subtypes.rds`。文件1,445,650,334字节，Seurat v5.4.0，37,795基因×41,105细胞；五研究、75文库。旧sample字段有44个水平，不能当作44名患者；映射当前curated_patient_id后共有40名患者。

对象中的直接记录：

- 2026-07-18 10:47：FindNeighbors(macro, dims=dims_use)，RNA pca的前30维，k=20，SNN图。
- 10:48：FindClusters，分辨率0.05、0.10、0.15、0.20；algorithm=1、random.seed=0。不能因全局元数据存在leiden字段就声称此次宏观亚群采用Leiden。
- 所有41,105个细胞的macro_cluster_final与RNA_snn_res.0.1完全一致；七个非空簇与七个亚型一一对应。
- 10:56–10:58：对macro_final再次LogNormalize(10000)、vst 3000特征、回归percent.mt、计算pca_macro_clean（50 PCs）、umap_macro_clean（前30 PCs，n_neighbors=30，min.dist=0.3，seed=123）。

后半段是清理后重新降维的记录。现存邻居图/聚类命令仍指向此前的pca，不能把较晚生成的“clean UMAP”当作重新聚类或已消除研究效应的证据。原聚类PCA的完整前处理命令已无法从这组覆盖后的commands中完全还原。

| 簇 | 旧标签 | 细胞数 | 最大来源研究 | 最大研究占比 |
|---|---|---:|---|---:|
|0|SPP1+Macro|11,290|GSE242299|93.92%|
|1|CX3CR1+Macro|10,823|GSE207493|95.76%|
|2|APOE+C1QC+Macro|6,030|GSE242299|93.33%|
|3|IL1B_IL8+Macro|5,911|GSE178481|99.71%|
|4|HLA-II+Macro|3,427|GSE156632|99.91%|
|5|CCL3+Macro|3,284|GSE159115|99.94%|
|8|Cycling+Macro|340|GSE207493|67.94%|

这种高度研究集中是直接观察结果，足以要求进一步验证；本轮没有据此将全部差异断言为纯技术效应，也不能称这些簇已经获得跨研究稳健验证。

## 二 输入同源及细胞清理链

旧`01_annotation_adata_annotated.h5ad`共有275,150细胞，标注Macrophage的细胞45,309个。最终七亚型对象的41,105个barcode全部存在于旧注释对象，并且全部在旧对象中标为Macrophage。因此“旧亚型对象与旧注释同源”已获细胞级确认。

旧注释巨噬细胞中有4,204个没有进入最终七亚型对象。最终对象的res.0.1字段保留12个factor水平，其中6、7、9、10、11为空；这与中间清理阶段相容，但不足以确定每个被删细胞原来属于哪一簇，更不能凭空补写删簇理由。未定位原始完整marker表、明确排除清单及生成saveRDS命令。

## 三 SPP1命名依据的证据尚不足

SPP1旧标签与cluster 0一一对应，说明它是一个已有聚类的命名，而非当前患者内上四分位阈值生成的标签。但是，执行“cluster 0→SPP1+Macro”的原始命名代码和全部marker筛选表没有找到。

本轮仅从保存的RNA data层复算描述性标记表达：旧SPP1群的SPP1平均表达约1.043、阳性比例30.59%；旧HLA-II群分别约1.526、66.33%。这说明对原七亚群而言，命名为SPP1的群并非SPP1表达最高的唯一群。该复核是原对象的描述性检查，不是新的亚群发现或独立生物学验证。

旧报告中的三节点DotPlot把其余五亚群合并为Other，比较对象已经改变。三节点图可以展示SPP1群相对合并Other的表达，但不能代替原七簇完整marker谱、患者一致性和研究效应核查。旧报告称SPP1群“跨全部五套数据集、最为稳健”，前半句仅表示存在：五研究数量依次28、47、44、567、10,604；后半句不能由这种存在性直接推出。

## 四 旧标签如何进入第三和第四部分

`run_macro_subtype_cellchat.R`读取p3_01_cellchat_input.rds和旧七亚型对象，按barcode回填宏细胞标签，构建七群版与三节点版。三节点保留SPP1、APOE/C1QC，其余亚型与未匹配巨噬细胞并入Other。代码在合并群体后运行CellChat，并不是按患者分别运行的新版Phase35/109分析。

旧`p4_02_build_reference.py`从旧注释H5AD的counts层出发，经临时label_cellchat_group.csv与label_macro_subtype.csv合并标签，以同样三节点规则建立空间参考；RegressionModel以study为batch。该study协变量存在于后续参考训练，不能倒推旧巨噬聚类已完成批次校正，也不会自动验证亚群命名。

这些代码和输出证明旧标签被复用到通信/空间参考，不能补回标签最初形成时的marker证据。

## 五 8月5日发生了明确的定义切换

历史会话与文件元数据相互吻合：

- 2026-08-05 18:45（+08:00）提出“第一部分重建结果对第四部分旧结论的参考体系敏感性分析”，明确不继承旧SPP1+ TAM标签；髓系改用连续SPP1与APOE/C1QC程序。
- 18:53:59创建build_patient_balanced_reference_v2.py；18:55:24完成reference-v2。源记录保存在reference_transition_trace/event_083.json、event_085.json、event_092.json。
- 输入是新版258,120细胞、75文库重建对象；生成229,032参考细胞。源码预设SPP1/GPNMB/LGALS3/FABP5/CD9/LPL与C1QA/C1QB/C1QC/APOE，用每患者程序分数百分位划定互斥上四分位状态。它的初始用途明确为参考敏感性试跑。
- 后续reference-v3、Phase8及Phase21等沿用这套状态定义。Phase21直接读reference_v2_cell_assignments.tsv.gz，没有先进行巨噬细胞无监督再分群。

这次切换有记录、有明确敏感性分析目的，并非旧七亚型被原样迁移。现有“SPP1-state”与旧“SPP1+Macro”同名相近但定义不同；应分别引用。

## 六 与当前重建主线逐细胞核对

旧41,105细胞中，38,692个在当前258,120细胞对象中精确匹配；其余2,413个均可在当前原始QC表定位且analysis_keep=False。其中2,367个有low_quality标记，55个有预测doublet标记，二者重叠9个；不能将两数直接相加。

旧SPP1群11,290细胞中：

- 10,276个保留于当前对象，1,014个已被当前QC排除。
- 保留细胞中9,232个当前仍标为Macrophage，1,044个标为其他类别，其中Epithelial cell 506个、T cell 459个。这是新旧注释不一致，不把当前标签自动当成实验真值。
- 9,232个当前巨噬细胞进一步按现行规则划为SPP1-state 1,889、C1QC-state 1,629、Other 5,714；即只有20.46%对应当前SPP1-state。

这个20.46%是两种定义的重叠率，不是以现行规则为真值计算出的旧标签准确率。它证明旧簇与当前预定义状态不能直接等同，也不能把旧分群图与当前轨迹/通信统计拼接为同一批同一定义细胞的连续证据。

## 七 能够恢复与仍缺失的内容

| 环节 | 当前判定 | 证据 |
|---|---|---|
|旧细胞来源同源性|确认|41,105/41,105匹配旧注释巨噬细胞|
|旧聚类确实运行|确认|Seurat commands及四分辨率字段|
|最终分辨率与簇标签对应|确认|res.0.1与最终簇41,105/41,105一致|
|七亚型完整marker命名链|未恢复|未见完整FindAllMarkers输出、命名代码|
|4,204细胞排除依据|未恢复|有数量差额及空factor水平，无完整排除记录|
|跨研究稳健性|不能确认|多数亚型被单一研究占据93%–99%以上|
|旧标签参与旧通信/空间|确认|下游脚本和日志直接引用|
|当前状态的定义与切换时间|确认|8月5日会话、源码和run_metadata一致|
|旧标签直接等同当前状态|否|细胞集合、注释及SPP1定义明显不同|
|原始操作者/完整执行过程|无法确认|目录账号、文件mtime和导入对象不能证明作者|

Git没有跟踪旧巨噬对象或相关macro文件；检索到的旧project会话主要是后续引用，并无对象生成过程。这里的“未恢复”限于本次列明的可访问资料，不等于证明其他电脑、其他工具或未保存文件中从未存在源代码。

## 八 对后续工作的具体影响

目前可以保留旧对象作为历史候选和标签敏感性对照；保留当前规则状态下已经完成的轨迹、regulon和通信结果，并明确它们的条件定义。旧七群结果尚不能直接升级为当前数据的正式发现链，当前程序状态也不能改名冒充无监督亚群。

若要完成用户要求的叙事，下一步应在当前确认的巨噬细胞宇宙中建立可追溯的精细分群：先检查巨噬/单核边界、QC和患者组成，进行不以SPP1阈值直接分组的聚类；比较完整marker谱、患者/研究分布及稳定性，再决定是否存在可命名的SPP1富集群；以患者级、尽量不重用定义基因的证据评价功能，最后核对并衔接下游分析。现阶段不需要推翻全部已有结果，但不能以补一张SPP1表达图替代这个环节。

本轮只完成记录与详细追溯；没有启动上述新分析，也未向原分析人员发送消息。

## 复核附件

- CURRENT_STATE_BEFORE_TRACE_CN.md：按用户要求先保存的原始断点。
- legacy_object/：对象结构、commands、软件环境、metadata和描述性marker表。
- identity/：逐细胞对应、研究组成、QC去向、旧/新标签列联表。
- SOURCE_SEARCH_FILES.json、SOURCE_SEARCH_HITS.json：源码/日志搜索边界与结果。
- project_session_trace/、reference_transition_trace/：相关历史事件与来源行号。
- SOURCE_FILE_MANIFEST.tsv：关键源文件SHA256，旧RDS原位校验。
- scripts/、logs/：本轮只读追溯代码与执行日志。
'''
(O/'MACROPHAGE_SPP1_PROVENANCE_REPORT_CN.md').write_text(report)
# Independently cross-check exported tables using standard-library CSV counts.
with gzip.open(O/'legacy_object/metadata.tsv.gz','rt') as h:a=list(csv.DictReader(h,delimiter='\t'))
with gzip.open(O/'identity/cell_label_crosswalk.tsv.gz','rt') as h:j=list(csv.DictReader(h,delimiter='\t'))
checks=[]
def ck(k,v,e):checks.append({'check':k,'actual':v,'expected':e,'pass':v==e})
ck('legacy_cell_count',len(a),41105);ck('unique_legacy_ids',len({r['cell_id'] for r in a}),41105);ck('old_libraries',len({r['sample_id'] for r in a}),75);ck('old_studies',len({r['study'] for r in a}),5)
ck('resolution01_matches',sum(r['macro_cluster_final']==r['RNA_snn_res.0.1'] for r in a),41105)
ck('seven_nonempty_subtypes',len({r['macro_subtype'] for r in a}),7)
ck('current_join_count',len(j),38692);ck('current_join_unique',len({r['cell_id'] for r in j}),38692)
spp=[r for r in a if r['macro_subtype']=='SPP1+Macro'];ck('old_spp1',len(spp),11290);ck('old_spp1_GSE242299',sum(r['study']=='GSE242299' for r in spp),10604)
x=[r for r in j if r['macro_subtype']=='SPP1+Macro'];ck('retained_old_spp1',len(x),10276);ck('retained_old_spp1_current_macro',sum(r['cell_type_broad_v1']=='Macrophage' for r in x),9232)
for state,n in [('SPP1-state Macrophage',1889),('C1QC-state Macrophage',1629),('Other Macrophage',5714)]:ck('old_spp1_to_'+state,sum(r['state_v2']==state for r in x),n)
with gzip.open(O/'identity/OLD_CELLS_ABSENT_CURRENT_QC.tsv.gz','rt') as h:m=list(csv.DictReader(h,delimiter='\t'))
ck('missing_in_QC',len(m),2413);ck('all_excluded',sum(r['analysis_keep']=='False' for r in m),2413);ck('low_quality',sum(r['low_quality']=='True' for r in m),2367);ck('doublets',sum(r['scrublet_predicted_doublet']=='True' for r in m),55);ck('overlap_low_quality_doublet',sum(r['low_quality']=='True' and r['scrublet_predicted_doublet']=='True' for r in m),9)
ck('parent_overlap',json.loads((O/'identity/LEGACY_PARENT_STATUS.json').read_text())['exact_parent_overlap'],41105)
ck('old_project_sessions',len(json.loads((O/'project_session_trace/SCANNED_PROJECT_SESSIONS.json').read_text())),17)
ck('reference_transition_sessions',len(json.loads((O/'reference_transition_trace/SCANNED_PROJECT_SESSIONS.json').read_text())),3)
assert all(c['pass'] for c in checks)
validation={'status':'PASS_TRACE_WITH_EXPLICIT_PROVENANCE_GAPS','checks':len(checks),'passed':sum(c['pass'] for c in checks),'new_clustering_run':False,'raw_results_modified':False,'remaining_gaps':['legacy full generator and marker table','legacy 4204-cell exclusion rationale','original author identity','current-data unsupervised discovery chain'],'details':checks}
(O/'TRACE_VALIDATION.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2))
print(json.dumps({k:v for k,v in validation.items() if k!='details'},ensure_ascii=False))
