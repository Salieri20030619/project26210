from pathlib import Path
import json,csv,hashlib,datetime,shutil
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=P/'phase144_parts2_3_4_audit_remote_v1';S=P/'local_records_migration_20260908_150455/project';N=O/'numerical_v2'
extra={
'S17_status':S/'phase17_part4_tymp_reset_local_v1/final_audit_sync_v1/PHASE17_MAPPING_QUEUE_STATUS_v1.json',
'S18_status':S/'phase18_part4_composition_niche_local_v1/server_results_sync_v1/PHASE18_COMPOSITION_NICHE_STATUS_v1.json',
'S19':S/'phase19_tymp_spp1_spatial_mainline_local_v1/server_results_v2_sync_v1/tables/primary_across_section_inference_v1.tsv',
'S97_plot_source':S/'phase97_spp1_vascular_spatial_bridge_local_v1/scripts_v1/make_spp1_vascular_bridge_v1.py',
'S97_validation':S/'phase97_spp1_vascular_spatial_bridge_local_v1/qc_v1/SPP1_VASCULAR_BRIDGE_VALIDATION_v1.json',
'Y03_analysis_source':S/'phase34_part2_y03_myeloid_pyscenic_local_v1/05_analyze_y03_candidate_regulons_v1.py'}
rows=list(csv.DictReader((N/'SOURCE_SHA256.tsv').open(),delimiter='\t'))
for k,f in extra.items():
 raw=f.read_bytes();(N/'evidence'/f'{k}{f.suffix}').write_bytes(raw);rows.append(dict(key=k,path=str(f),bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
with (O/'SOURCE_SHA256_ALL.tsv').open('w') as h:
 w=csv.DictWriter(h,fieldnames=['key','path','bytes','sha256'],delimiter='\t');w.writeheader();w.writerows(rows)
coverage=[
('Y01','第二部分','髓系四状态细分','已完成，可用','Phase21/33','50,761髓系细胞；40患者；五研究分开；不与Phase8的71,045巨噬细胞池混为同一分母'),
('Y02','第二部分','DPT与Monocle3轨迹','已完成，可用','Phase21/33','支持横断面状态排序；不是谱系转换；多数单研究n<10'),
('Y03','第二部分','髓系JUN/JUNB/CEBPD regulon','分析完成，可用；图件另审','Phase34','9/12端点FDR显著；CEBPD对照依赖、JUNB跨研究存在反向'),
('C01','第三部分','TYMP高低恶性细胞通信','已完成，可用','Phase35','14患者；3/4总通信方向FDR显著；不是SPP1专属通信'),
('C02','第三部分','受体映射与可用性','已完成，可用','Phase35','ITGAL/ITGB2为锁定主要映射；ITGAM/ITGB2仅敏感性'),
('C03','第三部分','NicheNet','已完成，可用','Phase35','完整靶基因响应正向；HIF1A端点未获FDR支持'),
('C04','第三部分','受体与TYMP/F1关联','已完成，以阴性/探索性呈现','Phase35','无FDR支持的主要受体关联；n=9敏感性P=0.01171875,q=0.09375'),
('C05','第三部分','TCGA联合生存','已完成，独立临床小节','Phase36','OS 508/168事件、PFS 506/156事件；分组预后不等于协同或治疗预测'),
('C109','第三部分补充','SPP1细分通信','已完成，可用','Phase109 analysis_v3','SPP1与内皮双向特异性支持；TYMP高低恶性端特异性不支持'),
('S01','第四部分','Cell2location映射','已完成，可用','Phase17','14切片，13状态；本轮核对完成记录，未重新训练/全量H5AD校验'),
('S02','第四部分','四通路空间代谢评分','已完成，方向性结果','Phase37','8项0项FDR显著；嘧啶10/12正向；不是实测通量或MTA浓度'),
('S03','第四部分','生态位与共定位','已完成，可用','Phase18-20','41,688 spots，42切片内生态位；空间主链异质'),
('S04','第四部分','空间受体互作','计算完成，部分不可估计','Phase37 v3 gate/v2 block means','6行完整保留；无支持；MST1R/IGDCC4缺失，PRTG仅2切片'),
('S05','第四部分','肿瘤-NAT病理边界','草稿敏感性完成；正式未完成','Phase38','已定位多边形未获病理审核；两切片方向相反'),
('S06','第四部分','未治疗/治疗后机制对照','仅描述性完成','Phase19/20','两切片不能提供可靠治疗效应推断，不作为配对疗效验证'),
('S70','第四部分补充','GSE175540空间关联','已完成，分制备类型可用','Phase70/97','SPP1状态程序与血管重塑FFPE/frozen均支持；不合并为24独立患者')]
with (O/'MODULE_COVERAGE.tsv').open('w') as h:
 w=csv.writer(h,delimiter='\t');w.writerow(['module','part','requirement','audit_state','formal_version','interpretation']);w.writerows(coverage)
report='''# 26210 第二、第三、第四部分已有结果审计

日期：2026-09-14。权威根目录：`/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1`。本审计为Phase144。

## 审计结论

第二部分已经形成“髓系状态细分—横断面轨迹排序—调控程序”的可用结果。第三部分已经形成患者级通信关联、SPP1与内皮的细分通信证据及独立的TCGA预后关联。第四部分完成主要空间计算，其中SPP1状态程序与血管重塑的关联具有分制备类型的一致性。

现有结果可以组织成这三部分，但当前数据不支持将其串成“肿瘤TYMP升高→SPP1极化→免疫排斥/治疗耐受”的普遍机制闭环。正式病理边界验证仍未完成。报告需要更新旧版完成度、明确空间端点名称，并保留关键阴性、反向及不可估计结果。

## 方法与覆盖范围

按历史清单Y01–Y03、C01–C05、S01–S06对应第二部分髓系、第三部分通信/临床、第四部分空间，加入改变结论的Phase109细分通信及Phase70/97空间补充。不是把第一部分报告中的2.2/2.3/2.4小节误作三大部分。

盘点旧part3/part4目录、现有重建阶段及服务器迁移快照，共2,167个文件。冻结46份来源，保留路径、文件大小和SHA256。使用患者级效应表复算主要轨迹/regulon/通信的样本数、方向、中位效应、双侧Wilcoxon和分家族BH；复核Cox表的HR、标准误、z、P及置信区间一致性；核对TCGA事件数及主要log-rank统计、空间汇总BH。修正审计读取精度并显式重现历史检验口径后，1,387/1,387项数值检查通过。检查数是字段级核验数，不是独立生物学假设数。

独立统计单位按患者或已登记的切片身份；细胞、spot、邻域、技术重复不替代样本量。五研究、GLMF、TJRCC、GSE175540不同制备类型及TCGA分别记录。历史75文库结论不自动外推为所有参考分支的共同分母。

## 第二部分：状态排序和调控差异可用，SPP1替代C1QC的叙事不成立

- Y01/Y02覆盖50,761髓系细胞、40患者、五个独立建图研究。Monocle3患者级SPP1-minus-Monocyte拟时序效应中位数在五研究均为正。GSE207493为12/12患者正向，中位差0.5347，P=0.0004883，FDR=0.0006510。
- 其余研究有效配对n=5–7，属于探索性结果：GSE156632仅3/5方向正、P=0.3125；GSE178481虽5/5方向正，P=0.0625，不能写显著。GSE159115为6/6正向；GSE242299为7/7正向。五研究方向一致不等于全部患者一致或五研究均显著。
- 排除SPP1本身的状态程序随Monocle3排序仅4/5研究方向正；SPP1/C1QC相对程序在五研究、两种轨迹方法均下降。保留“Monocyte到SPP1状态的横断面排序”；删除“SPP1持续替代C1QC”的单一路径表述。
- Y03三个候选TF×四端点中9/12通过对应端点三TF BH。CEBPD随Monocle3排序正向（患者内rho中位数0.2721，29/40正向，P=1.0405e-5，q=3.1215e-5）；SPP1相对Monocyte的AUCell中位差+0.003551，29/35正向，P=6.5051e-5，q=0.0001952；相对C1QC却40/40降低，中位差−0.005027，P=1.8190e-12，q=5.4570e-12。不能称为SPP1特异激活。
- JUNB在四个总体端点均负向并通过FDR，但GSE178481的轨迹及SPP1-Monocyte比较均反向为正（该研究n=5）；正文保留总体下降，补充材料展示研究异质性。JUN只有SPP1相对C1QC降低获支持，其他端点当前数据不支持。

## 第三部分：通信结果已补齐，受体机制与临床协同需分开判定

**C01：TYMP高低恶性细胞通信。** 14患者中，内皮→恶性细胞为13/14正向、中位高低差0.1559、P=0.0006104、q=0.002441；巨噬→恶性12/14正向、差0.07994、P=0.001709、q=0.003418；恶性→巨噬11/14正向、差0.2713、P=0.008545、q=0.01139。恶性→内皮仅8/14正向、P=q=0.07849。上述为总通信兼容性，不等价于某一配体的实验作用。

**Phase109补充已经解决旧记录中的细分通信缺口。** 正式结果为analysis_v3，前两版失败保留。SPP1状态→内皮相对C1QC/Other的特异性中位差0.21065，14/14正向，P=0.0009815，q=0.0030466；内皮→SPP1中位差0.14642，13/14正向，P=0.0015233，q=0.0030466，两个方向的留一患者中位数均保持正向。恶性→SPP1和SPP1→恶性的TYMP高低特异性均未通过（q=0.2642、0.3967）。

**肿瘤轴与SPP1丰度。** 20名符合门槛患者的四个预设端点均未获支持：严格TYMP轴与SPP1比例/相对C1QC比值估计分别−0.1911、−0.1776；TYMP RNA对应+0.2281、+0.0404；q分别0.6775、0.6775、0.6775、0.8822。保留方向，不能将内皮通信阳性用于替代这条上游阴性结果。

**C02–C04：受体与下游。** ICAM5–ITGAL/ITGB2是项目资源锁定的主要映射；ITGAM/ITGB2单列敏感性。NicheNet完整靶基因响应与ICAM5/MST1预测权重正相关（中位r=0.06761/0.04025；14/14、13/14正向；P=0.001097/0.003173，q=0.002194/0.003173）；HIF1A加权响应多数患者负向且均未获FDR支持。C04主要受体关系没有FDR支持；ITGAM/ITGB2敏感性与TYMP的n=9结果为8/9正向，rho中位数0.02385，P=0.01171875，q=0.09375，只能称名义显著。F1–ITGB2 n=17，P=0.04477、q=0.35815，同样不能称校正后显著。

**C05：TCGA预后。** 510名唯一患者；OS 508人/168事件，PFS 506人/156事件。预设六受体轴×两终点的分组log-rank有11/12通过各终点内六轴BH，MST1R的PFS不支持。该结果说明分组预后差异，并不证明受体与TYMP正向协同。

年龄/性别/分期调整后的Both-high对Both-low：ICAM5–ITGAL/ITGB2的OS HR=1.5096，P=0.03999，但原Cox表没有该组比较的多重校正。本审计另外按“终点×模型×对比、六轴”计算审计性BH，q=0.05999，因此该特定HR只能标名义显著。CD209的OS、PRTG的OS/PFS、IGDCC4的OS/PFS在这一审计口径下仍通过BH（完整表附后）。这些新q是事后审计标注，不冒充原预设主要检验。

连续交互大多不支持；IGDCC4的OS交互HR=0.7679，P=0.001843，审计q=0.01106，为负向；PFS交互HR=0.8519，P=0.04694、审计q=0.2816。不得据分组生存图声称正协同、已验证治疗靶点或ICI疗效预测。

## 第四部分：基础空间结果完成，可靠阳性位于SPP1状态与血管重塑

**S01/S03基础。** Phase17既有完成记录14/14映射成功；Phase18为14切片、41,688 spots、13参考状态、42个切片内生态位。TYMP空间表达来自混合spot，只能称恶性区域TYMP暴露，不能独占归因于恶性细胞。F1为辅助连续变量。

**TYMP–SPP1主链不获支持。** TJRCC12中邻域TYMP与SPP1/C1QC表达极化仅3/12正向，rho中位数−0.03382，P=0.17627、q=0.26440；与SPP1组成极化7/12正、rho=0.08971、P=q=0.79102；与排除SPP1的状态程序4/12正、rho=−0.27206、P=0.09570、q=0.26440。不得恢复旧参考版本的普遍极化主张。

**血管及免疫层。** TJRCC12的SPP1极化–血管组分仅6/12方向正、q=0.7910；生态位–EMT为8/12正，rho=0.2147、P=0.020996、q=0.10498，保留名义显著；细胞毒性关系仅3/12符合预期负方向，当前数据不支持统一免疫抑制解释。

GSE175540分制备类型的组成调整分析支持SPP1状态程序（排除SPP1基因）–血管重塑：FFPE中位rho=0.21206，frozen=0.20835，各12切片、12/12正向，P均0.00048828125，q分别0.0021973、0.0043945。两个制备类型分开报告，不能当作24个新增独立患者。TYMP–SPP1状态程序在同一调整模型中均不支持：FFPE rho=0.00995、P=q=0.063965；frozen rho=−0.00305、P=0.85010、q=0.96973。SPP1–T细胞功能障碍仅FFPE支持，不能泛化为跨制备免疫排斥。

**S02空间代谢。** 自定义四轴与KEGG四轴分为两个家族，8项均未获FDR支持。自定义嘧啶轴10/12正向、中位rho=0.43235、P=0.092285、q=0.36914，是有方向一致性的候选线索。评分由表达与邻域平滑构建，不是实测代谢通量，也不是MTA浓度。

**S04空间LR。** 正式采用v2空间块均值及v3稀疏性完成门禁，不能引用最初median-block结果。六条预设记录完整保留，没有受支持的空间LR：ICAM5–ITGAL/ITGB2有效n=10、rho=0.07846、P=0.49219；CD209 n=7、rho=0.10728、P=0.765625；PRTG仅n=2、方向均负；MST1R及IGDCC4缺失，n=0。不可估计行的P=1或q=1是校正/占位值，不能解释为完成了真实阴性检验。

**S05/S06。** 已定位的边界仍为not_pathologist_reviewed草稿。GLMF1/2距离边界1000 μm内各922/1010 spots，嘧啶空间评分朝肿瘤方向在GLMF1上升（rho=0.577），GLMF2下降（rho=−0.501）。正式病理验证未完成。两切片的处理状态对照只能描述，不能估计可靠治疗效应。

## 必须落实的审计处理

1. **更新完成度。** Phase24旧表中的Y02、Y03、C01–C05、S02/S04未完成记录已被后续阶段替代；Phase105所说SPP1细分通信缺口已由Phase109补齐。新报告不得继续复制这些旧待办。
2. **区分统计支持标签。** 小样本表的statistical_support=True可以仅指P<0.05。引用时必须同时查看n、reporting_tier、P和FDR。保留原表，审计提供解释，不改写源结果。
3. **Cox校正口径。** 主log-rank显著性不能转移给某个调整后HR或交互项；本审计附单独的事后BH敏感性。
4. **空间图端点要写清。** Phase97空间底图使用q05 SPP1巨噬细胞丰度；右侧血管统计来自SPP1_state_no_SPP1表达程序的组成调整模型。两者是不同展示/检验端点，图注不能声称气泡rho就是左侧丰度共定位的相关系数。Phase97的TYMP邻域统计与Phase70组成调整统计也不是同一端点，数值不同不构成矛盾；不得混引。
5. **软件复现细节。** 初轮2项P差异已解释：CEBPD来自审计CSV读取精度引入的并列，round_trip读取解决；JUN的SPP1-C1QC比较原P=0.0005080，SciPy1.17.1自动方法给P=0.0007784，源于并列秩时算法选择。显式历史exact可重现原值，现代auto敏感性仍通过同家族BH；原记录与敏感性均保留。
6. **完成计算不等于证实假设。** S04稀疏端点及S05正式病理应保留未估计/未完成状态；不以技术PASS升级科学结论。

## 讨论与审计边界

本轮为已有结果、原脚本与统计表审计，没有重跑CellChat、pySCENIC、Monocle3、Cell2location或Cox模型；Cox检验复算属于参数表内部一致性及队列核对，不是从原始表达重新拟合。对大型H5AD未重新做全量hash和重新注释，对旧Word/PDF未做逐页视觉审计。Phase17/18的大对象正确性引用既有验证并核对当前记录；不能把本轮字段检查说成全流程复现。

轨迹不证明谱系转换，AUCell不测量TF结合，CellChat/NicheNet不测量配体结合，空间共定位不证明作用方向。临床队列为观察性bulk，不能将基因表达指定给某类细胞或直接推断治疗响应。以上边界不否定已建立的状态差异及患者级方向关联。

客户湿实验由客户负责。本轮不重分析、不审计客户原始数据；按最新说明保留shTYMP使MTA等下降、Thymidine上升的背景，不沿用历史MTA升高或44倍说法。本次并非外部免疫治疗、公共干预、细胞系代谢全部分支的重新审计；这些独立数据宇宙不混入本三部分的主要统计家族。全表包括关键阴性与反向记录。

## 附件与复核入口

- `MODULE_COVERAGE.tsv`：14项原清单及2项关键补充的逐项完成状态。
- `SOURCE_SHA256_ALL.tsv`：46份冻结来源与服务器路径。
- `numerical_v2/FILE_INVENTORY.tsv`：2,167个盘点文件元数据。
- `numerical_v2/NUMERICAL_CHECKS.tsv`：1,387项字段级数值核验。
- `numerical_v2/Y03_SCIPY_AUTO_SENSITIVITY.tsv`：历史exact与现代auto差异。
- `numerical_v2/C05_cox_AUDIT_BH.tsv`及`C05_interaction_AUDIT_BH.tsv`：审计性多重校正。
- `numerical_v2/evidence/`：完整源统计表（含P/FDR和患者级内部审计数据，仅服务器保存）。
- `01_audit.py`、`02_audit_precision_and_extended_checks.py`：审计代码；初轮结果和日志保留。
'''
(O/'PARTS2_3_4_AUDIT_CN.md').write_text(report,encoding='utf-8')
# Render simple Markdown tables of every primary summary without requiring document software.
keys=['Y02','Y03','C01','C03','C04','C109','A109','C05_logrank','S19','S20','S02','S04','S70_meta']
text=['# 完整主要统计表\n\n以下保留各来源的原字段。不可估计占位值及小样本标志按审计正文解释。\n']
for k in keys:
 f=N/'evidence'/f'{k}.tsv';rr=list(csv.reader(f.open(),delimiter='\t'))
 text += ['\n## '+k+'\n', '| '+' | '.join(rr[0])+' |','| '+' | '.join(['---']*len(rr[0]))+' |']
 text += ['| '+' | '.join(v.replace('|','/') for v in r)+' |' for r in rr[1:]]
(O/'ALL_PRIMARY_STATISTICS.md').write_text('\n'.join(text),encoding='utf-8')
status=json.loads((N/'AUDIT_COMPUTATION_STATUS.json').read_text());assert status['checks']==1387 and status['passed']==1387
assert len(coverage)==16 and len(rows)==46
for r in rows:
 assert hashlib.sha256(Path(r['path']).read_bytes()).hexdigest()==r['sha256']
validation={'status':'PASS_RESULT_TABLE_AUDIT_WITH_DOCUMENTED_INTERPRETATION_FINDINGS','checks':1387,'passed':1387,'source_hashes_rechecked':46,'scope_rows':16,'model_rerun':False,'source_results_modified':False,'local_records_written':False}
(O/'REPORT_VALIDATION.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2))
print(json.dumps(validation,ensure_ascii=False))
