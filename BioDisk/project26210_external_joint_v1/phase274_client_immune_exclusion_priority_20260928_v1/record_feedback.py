from pathlib import Path
import pandas as pd,json,hashlib,datetime,shutil,csv
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');p=r/'phase274_client_immune_exclusion_priority_20260928_v1'
for s in ['delivery','qc','source_snapshots']:(p/s).mkdir(parents=True,exist_ok=True)
a=r/'phase246_spatial_marker_structure_20260924_v1/tables'
src=a/'MARKER_STRUCTURE_COHORT_SUMMARY_WITH_EFFECTIVE_N.tsv';x=pd.read_csv(src,sep='\t')
x[x.endpoint.str.contains('CD8|Fibro',case=False,na=False)].to_csv(p/'delivery/EXISTING_CD8_STROMA_EVIDENCE.tsv',sep='\t',index=False)
src2=a/'COMPOSITE_RECHECK_COHORT_SUMMARY_WITH_EFFECTIVE_N.tsv';y=pd.read_csv(src2,sep='\t')
y[y.endpoint.str.contains('Tcell|T_cell|T.cell|T_balance|balance',case=False,na=False)].to_csv(p/'delivery/EXISTING_TCELL_BALANCE_EVIDENCE.tsv',sep='\t',index=False)
src3=r/'phase234_part4_replan_execution_20260921_v1/tables/TECHNICAL_SCOPE_AND_UNESTIMABLE.tsv'
shutil.copy2(src3,p/'delivery/EXISTING_IDENTITY_BOUNDARY_GAPS.tsv')
src4=r/'phase234_part4_replan_execution_20260921_v1/tables/EXPRESSION_BOUNDARY_BIN_PROFILES.tsv';z=pd.read_csv(src4,sep='\t')
boundary={'frozen_partition':'malignant_fraction >=0.20; expression-defined, not audited histology','original_scales_um':[200,400,800],'profile_bin_um':200,'existing_features':sorted(z.feature.unique()),'missing_from_prior_boundary_profiles':['T-cell abundance','CD8 lineage signal'],'cohort_patient_coverage':z.groupby('dataset').patient.nunique().to_dict()}
(p/'delivery/EXISTING_BOUNDARY_SCOPE.json').write_text(json.dumps(boundary,ensure_ascii=False,indent=2))
report="""# 客户反馈：免疫空间排斥的证据梳理与实施顺序
日期：2026-09-28。反馈登记及资产核查；本轮未新增统计、未改稿或绘图，正式稿仍Phase273。

客户评价原文：
> 这个要跟免疫扯上关系的
> 屏障什么的应该可以看吧
> 比如把T细胞挡在外面什么的

## 对主线的影响
优先检验SPP1相关巨噬富集区域是否伴随T细胞在肿瘤外侧积聚及向内部浸润减少。临床关联继续保留为转化背景；免疫空间关系承担更直接的组织层联系。免疫“排斥样分布”和由巨噬“造成物理屏障”的因果机制分层解释。

## 已有证据及明确缺口
1. Phase246的GSE175540表达区域距离结果：有效16患者，12正/3负/1零，SPP1高比例区域到CD8表达富集区的距离差中位数+300 μm，P=0.03516，FDR=0.05273，属于名义显著。它支持进一步研究区域分离，尚未区分T细胞位于肿瘤外还是肿瘤内其他位置。
2. 同来源固定400 μm邻域中，SPP1巨噬内部比例与CD8谱系程序20患者16负，rho中位-0.2814，P=0.01182、FDR=0.03545；完整组成/深度等调整后rho=-0.09162，P=0.2632、FDR=0.3948。调整后减弱必须保留。
3. TJ-RCC直接T细胞绝对丰度未显示一致降低；T细胞相对平衡降低不能等同数量减少。其CD8距离端点仅5例可评、方向未重复。成纤维/ECM邻近没有一致支持，不能预设SPP1-CAF双层屏障。
4. 已有Phase234冻结的肿瘤表达分区（恶性构成≥0.20）及有符号组织图距离，可复用；原界面曲线只包含SPP1四种读出、总巨噬、POSTNCAF及ECM，没有T-cell或CD8。这是本轮找到的具体补充位置，不需要从头重建全部空间流程。
5. 当前可用H&E图像并不等于已有审核后的病理边界；原资格表仍将病理边界和独立CD8细胞定位列为未获资格。总T细胞有既有映射，CD8为谱系RNA区域读出。旧空间SPP1-state不能直接改名为新C4。

## 下一阶段的实施顺序
第一步：复核表达分区、组织图连通性和H&E对应。保留原恶性构成阈值，分区不能用SPP1或T细胞信号来圈定；记录边界碎裂、切片截断、内外侧覆盖不足。若缺乏可信病理标注，明确采用“肿瘤表达界面”，仍可开展表达区域层面的探索，不冒称病理肿瘤巢边界。

第二步：在既有界面框架中补入总T细胞丰度和CD8谱系信号。主宽度沿用400 μm，200/800 μm为固定敏感性；以200 μm有符号距离分箱呈现从外侧到内侧的连续变化，不跨组织空洞。各侧需达到既定区域覆盖，不按结果修改分区、宽度或纳入名单。总T细胞丰度内外差和CD8表达梯度分别解释，不能把表达分数换称细胞密度。

第三步：联合检验界面SPP1富集与T细胞进入程度。比较边界内外的T细胞分布，并评价SPP1界面富集与内部浸润的患者层联系；同时保留总巨噬丰度、组织组成和深度等控制。CAF/胶原只作为具备资格时的辅助结构证据，不作为必须“做出阳性”的部分。正式检验前按实际可用性锁定主端点、协变量和同类FDR家族。

第四步：先做GLMF1/2结构核对及技术展示，再按合格患者扩展GSE175540与TJ-RCC；两个GLMF病例只作描述。各来源独立统计，同患者多切片不增加样本量，n<10不作推断检验。案例按边界和覆盖资格选择，不按T细胞是否在外侧选择；展示总体结论时同时保留反向案例。

建议图件由三层组成：H&E/表达界面与SPP1及T细胞叠图；界面外到内的患者等权梯度；患者层SPP1富集与T细胞浸润效应图。现有三元距离图不替代这组内外侧证据。

## 结果措辞
若空间数据支持，可写“SPP1相关巨噬在肿瘤表达界面富集，并伴随T细胞外侧积聚和内部浸润减少，提示免疫排斥样空间组织”。只有空间分离或相对构成变化时，按对应层级陈述；不能仅凭spot图称单细胞接触、物理阻挡或已证实免疫屏障机制。

本轮交付包括既有CD8/基质完整相关端点、T细胞相对平衡端点及边界资格记录，不覆盖原统计。上述实施顺序是后续工作优先级，尚未作为已完成的新结果写入正文。
"""
(p/'delivery/客户反馈_免疫排斥证据与实施顺序.md').write_text(report,encoding='utf-8')
sources=[src,src2,src3,src4,r/'phase246_spatial_marker_structure_20260924_v1/REPORT_CN.md',r/'phase234_part4_replan_execution_20260921_v1/scripts/freeze.py',r/'phase234_part4_replan_execution_20260921_v1/scripts/narrow.py']
sourceinfo=[dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in sources]
(p/'qc/SOURCE_MANIFEST.json').write_text(json.dumps(sourceinfo,indent=2))
now=datetime.datetime.now().astimezone().isoformat()
val={'status':'FEEDBACK_AND_ASSET_REVIEW_COMPLETE','time':now,'new_tests':0,'manuscript_unchanged':'Phase273','core_gap':'T/CD8 not included in existing expression-boundary profiles','new_claims_of_barrier':False}
(p/'qc/FINAL_VALIDATION.json').write_text(json.dumps(val,indent=2));shutil.copy2(p/'qc/FINAL_VALIDATION.json',p/'delivery/FINAL_VALIDATION.json')
manifest=[]
for f in sorted((p/'delivery').iterdir()):
 if f.is_file():manifest.append(dict(file=f.name,bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
(p/'delivery/DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
record=f"""## Phase274（2026-09-28）：客户反馈转向免疫空间排斥，资产核查完成
客户新增评价要求联系免疫，考虑屏障及“T细胞挡在外面”。已登记为下一优先问题：SPP1相关巨噬界面富集与T细胞外侧积聚/内部浸润。临床图保留转化背景，不将临床关联替代屏障空间证据。
已复核Phase234/246/250/260：GSE175540 CD8距离差+300um，n16，12正3负1零，P0.03516/q0.05273；400um邻域CD8负关联调整后减弱，TJ-RCC未形成重复，直接T绝对丰度未一致下降，CAF/ECM未支持普遍邻近。关键具体缺口：Phase234已冻结malignant_fraction>=.20表达分区及200um距离分箱，但边界读出仅SPP1/总巨噬/POSTNCAF/ECM，没有T或CD8。可复用其骨架补界面内外免疫分布，避免重建既有流程。
当前无审核病理边界/独立CD8定位；总T映射可用，CD8属于谱系RNA区域，新C4与旧SPP1-state仍分开。后续先核对界面资格/H&E，再补400um主宽度+200/800um敏感性，患者单位，控制总巨噬/组成/深度，保留反向，不按结果选例。可评价免疫排斥样结构，不能预先宣称物理阻挡因果。
目录phase274_client_immune_exclusion_priority_20260928_v1，delivery含反馈原文、证据梳理、实施顺序和既有端点表。本轮只登记反馈/盘点，不新增检验、图件或改稿，主稿仍Phase273。无后台任务。
"""
bk=p/'qc/root_record_backup';bk.mkdir(exist_ok=True)
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 f=r/name;shutil.copy2(f,bk/name);f.write_text(record+'\n'+f.read_text(),encoding='utf-8')
f=r/'RESULTS_MANIFEST.tsv';shutil.copy2(f,bk/f.name)
with f.open('a') as out:
 wr=csv.writer(out,delimiter='\t',lineterminator='\n')
 for it in manifest:wr.writerow(['phase274_'+it['file'],'phase274','client_feedback_and_existing_evidence_audit',str(p/'delivery'/it['file']),now,'PASS',str(p/'qc/FINAL_VALIDATION.json'),'No new test; next priority immune exclusion at tumor expression interface'])
print('REGISTERED',len(manifest),'files; CD8/stroma rows',len(x[x.endpoint.str.contains('CD8|Fibro',case=False,na=False)]))
