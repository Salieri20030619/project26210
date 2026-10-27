from pathlib import Path
import json
import pandas as pd,numpy as np
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');P=R.parent;O=R/'09_downstream_crosswalk/remaining_analysis_v1';O.mkdir(exist_ok=False)
assert (R/'08_regulons/statistics_v1/STATUS.json').exists()
new=pd.read_csv(R/'07_trajectory/input_v1/CELL_METADATA.tsv.gz',sep='\t',dtype={'candidate':str},index_col=0,low_memory=False);old=pd.read_csv(P/'phase34_part2_y03_myeloid_pyscenic_v1/input_v1/cell_metadata_v1.tsv.gz',sep='\t',index_col=0,low_memory=False);assert new.index.isin(old.index).all();m=new[['candidate','study','curated_patient_id','sample_id','source_tissue']].copy();m['legacy_state']=old.loc[m.index,'state_v2'].astype(str);m.to_csv(O/'NEW_TO_LEGACY_CELL_LABELS.tsv.gz',sep='\t',index_label='cell_id')
pd.crosstab(m.candidate,m.legacy_state).to_csv(O/'NEW_CANDIDATE_BY_LEGACY_STATE.tsv',sep='\t');pd.crosstab(m.loc[m.candidate=='4','study'],m.loc[m.candidate=='4','legacy_state']).to_csv(O/'C4_LEGACY_STATE_BY_STUDY.tsv',sep='\t')
rows=[]
for state in sorted(set(old.state_v2.astype(str))):
 A=set(m.index[m.candidate=='4']);B=set(old.index[old.state_v2.astype(str)==state]);hit=A&B;rows.append({'legacy_state':state,'new_C4_cells':len(A),'legacy_cells_all50761':len(B),'intersection':len(hit),'C4_fraction_in_legacy_state':len(hit)/len(A),'legacy_state_fraction_recovered_as_C4':len(hit)/len(B),'Jaccard':len(hit)/len(A|B)})
pd.DataFrame(rows).to_csv(O/'C4_LEGACY_OVERLAP_METRICS.tsv',sep='\t',index=False)
deps=[['Part2_old_trajectory','Old graph and root results remain historical; new cells require new graph','completed_new_partitioned_model','07_trajectory/monocle3_v2'],['Part2_old_regulon','Old rule-stratified discovery remains sensitivity; new network learned from new catalogue','completed_new_network','08_regulons/network_v1'],['Part3_broad_malignant_macrophage_endothelial','Coarse cell communication remains valid for its original broad-cell endpoint','retain_original_no_automatic_rerun','phase35'],['Part3_SPP1_fine_communication','New C4 identity differs from legacy SPP1 score state; need C4 patient-specific interaction re-estimation before using old effect as C4','future_dependency_only','phase109'],['Part3_TCGA_receptor_TYMP','Bulk clinical endpoint independent of new macrophage cluster labels','retain_original_no_automatic_rerun','phase36'],['Part4_legacy_cell2location','Old SPP1-state posterior cannot be renamed C4; new C4 reference/coverage audit and deconvolution required if exact new abundance requested','future_dependency_only','phase17_18_78'],['Part4_expression_programs','Old published/fixed programs remain their original endpoint; new C4 program transfer must be frozen and tested separately','future_dependency_only','phase70_97'],['Client_wetlab','Client-owned validation and direction, not reanalysed','outside_analysis_scope','user materials']]
pd.DataFrame(deps,columns=['module','decision','status','source']).to_csv(O/'DOWNSTREAM_DEPENDENCY_MATRIX.tsv',sep='\t',index=False)
(O/'INTERPRETATION_CN.md').write_text('''# 新旧标签与第三/四部分衔接
新C4是冻结无监督候选，旧SPP1-state是患者内程序上四分位规则；两者不能改名等同。逐细胞重叠、研究组成、Jaccard及分母均在本目录。
原有第三部分广义细胞通信和独立bulk临床分析可保留其原端点。要把“旧SPP1状态—内皮”结论换成“新C4—内皮”，需按新C4标签重新估计患者级细分通信；当前仅登记依赖。
空间旧SPP1后验丰度不能直接改为C4。若使用新C4精确丰度，需要新reference覆盖审查和映射；若保留旧表达程序，可保留原名称和原统计，不能声称验证新C4细胞丰度。
第二部分轨迹为横断面排序，regulon为推断调控程序，均不证明肿瘤TYMP诱导C4。新状态与原主线的衔接通过下一阶段明确端点检验，不使用下游阳性替代上游关键阴性。此次没有运行第三/四部分新模型，没有绘图。
''')
(O/'STATUS.json').write_text(json.dumps({'status':'COMPLETE_CROSSWALK_AND_DEPENDENCY_REVIEW','new_cells':len(m),'new_C4_cells':int((m.candidate=='4').sum()),'legacy_cells':len(old),'downstream_models_rerun':False,'figures_created':False},indent=2));print('CROSSWALK_DONE',flush=True)
