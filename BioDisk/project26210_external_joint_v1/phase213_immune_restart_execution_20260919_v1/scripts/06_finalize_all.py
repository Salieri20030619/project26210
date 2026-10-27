from pathlib import Path
import json,datetime,hashlib,csv,shutil
import numpy as np,pandas as pd
H=Path.home();R=H/'project26210_external_joint_v1';O=R/'phase213_immune_restart_execution_20260919_v1';F=R/'phase211_immune_receiver_refreeze_20260919_v1/freeze_v2';stages=['I1_patient_inputs','I2_main_statistics','I3_robustness','I4_directed_communication','I5_external_description']
statuses={s:json.loads((O/s/'STATUS.json').read_text()) for s in stages};assert all(v['status']=='COMPLETE_VALIDATED' for v in statuses.values());seal=json.loads((F/'FREEZE_SHA256.json').read_text());assert all(hashlib.sha256((F/k).read_bytes()).hexdigest()==v for k,v in seal.items())
main=pd.read_csv(O/'I2_main_statistics/PRIMARY_12_HYPOTHESES.tsv',sep='\t');lr=pd.read_csv(O/'I4_directed_communication/ALL_DIRECTED_HYPOTHESES.tsv',sep='\t');ext=pd.read_csv(O/'I5_external_description/EXTERNAL_STRATIFIED_DESCRIPTIVE_SUMMARY.tsv',sep='\t');infl=pd.read_csv(O/'I3_robustness/INFLUENCE_SUMMARY.tsv',sep='\t');sens=pd.read_csv(O/'I3_robustness/PRESET_SENSITIVITY_EFFECTS.tsv',sep='\t');sources=pd.read_csv(O/'I4_directed_communication/C4_VS_OTHER_SOURCE_PATIENT_DIFFERENCES.tsv',sep='\t');pat=pd.read_csv(O/'I4_directed_communication/DIRECTED_PATIENT_METRICS.tsv',sep='\t')
# Add sample-information accounting without altering the original effect/P/FDR tables.
info=[]
for row in lr.itertuples():
 g=pat[pat.hypothesis_id==row.hypothesis_id];valid=[]
 for study,a in g.groupby('study'):
  if len(a)>1 and a.x.nunique()>1 and a.y.nunique()>1:valid.append((study,len(a)))
 info.append(dict(hypothesis_id=row.hypothesis_id,n_patients=row.n_patients,n_studies=g.study.nunique(),informative_studies=len(valid),informative_patients=sum(n for _,n in valid),study_counts=json.dumps(g.study.value_counts().to_dict())))
pd.DataFrame(info).to_csv(O/'I4_directed_communication/DIRECTED_STUDY_INFORMATION.tsv',sep='\t',index=False)
source_summary=[]
for li,g in sources.groupby('lineage'):source_summary.append(dict(receiver_scope=li,n_patients=len(g),median_C4_minus_other_SPP1_log2CPM=float(g.C4_minus_other.median()),positive_patients=int((g.C4_minus_other>0).sum()),negative_patients=int((g.C4_minus_other<0).sum()),interpretation='existing source-state difference, not evidence of receiver suppression'))
source_summary=pd.DataFrame(source_summary);source_summary.to_csv(O/'I4_directed_communication/C4_SOURCE_DESCRIPTIVE_SUMMARY.tsv',sep='\t',index=False)
pred_cn={'mature_SPP1_pseudobulk':'整体成熟巨噬SPP1','C4_fraction':'C4比例','C4_SPP1_pseudobulk':'C4内SPP1'};ep_cn={'CD8_cytotoxicity':'CD8细胞毒','CD8_dysfunction':'CD8功能障碍相关程序','Treg_CD4':'Treg/CD4','NK_cytotoxicity':'NK细胞毒'}
def num(v):return 'NA' if pd.isna(v) else f'{v:.5g}'
def table(headers,rows):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rows])
main_table=table(['预测量','免疫端点','独立患者n','研究调整rho','95% bootstrap区间','双侧置换P','FDR'],[[pred_cn[r.predictor],ep_cn[r.endpoint],r.n_patients,num(r.rho),f'[{num(r.bootstrap_CI_low)}, {num(r.bootstrap_CI_high)}]',num(r.p_value),num(r.FDR) if r.adjustment=='BH' else '不适用：预设单终点'] for r in main.itertuples()])
source_table=table(['接收端配对范围','患者n','C4−其他巨噬SPP1中位差','正向患者数'],[[r.receiver_scope,r.n_patients,num(r.median_C4_minus_other_SPP1_log2CPM),r.positive_patients] for r in source_summary.itertuples()])
comm_table=table(['接收端','受体/先验靶表达','来源预测量','n','调整rho','P','FDR'],[[r.lineage,r.measure,r.predictor,r.n_patients,num(r.rho),num(r.p_value),num(r.FDR)] for r in lr[lr.n_patients>0].itertuples()])
ext_table=table(['治疗/取材层','端点','患者n','描述性rho','最高−最低SPP1患者的端点差'],[[r.treatment,ep_cn[r.endpoint],r.n_patients,num(r.rho_descriptive_only),num(r.outcome_difference_highest_vs_lowest_SPP1)] for r in ext.itertuples()])
report='''# 26210 免疫分析重启：I1–I5完整结果

阶段：Phase213。分析日期：2026-09-19。执行依据：Phase212预设计划及用户“按顺序逐个展开，出去期间不用汇报”的明确授权。I1患者指标、I2主统计、I3有限检查、I4定向通信、I5外部描述性对应已依次完成并验收。

## 摘要：当前数据不支持统一的SPP1相关免疫抑制关联

在冻结接收端与同肿瘤文库配对范围内，整体成熟巨噬SPP1与CD8细胞毒、CD8功能障碍相关程序、Treg/CD4及NK细胞毒均未获得研究调整后的统计支持。12项预设主要/次要关联的调整后P均>0.05。

C4比例与CD8功能障碍存在未经研究调整的名义正相关，并在留一患者/研究时保留正方向；主模型调整研究来源后该关联减弱且不显著。此项作为方向性候选线索保留，不提升为已证实免疫抑制。C4的SPP1富集来源特征仍在，不因下游免疫端点未支持而删除。

受体表达支持有限的SPP1接收可能性，但没有形成获得FDR支持的受体协同或先验靶表达联系。唯一名义受体关联为负向，完整保留。外部各治疗层1–4患者只作描述，没有正式重复验证。

## Methods

### 数据与统计单位

采用Phase211/freeze_v2逐细胞身份、患者和文库名单；成熟巨噬使用Phase154继承的净化定义，排除candidate1及Monocyte。原始75文库发现对象为258120细胞；本轮不改变源对象、前两部分或C4的2674个冻结成员。CD8两端点纳入16患者/19文库/4研究，Treg/CD4和NK纳入13患者/16文库/2研究。C4群内表达分别使用11和10名患者的独立文库集合。

SPP1按成熟巨噬或C4原始计数伪汇总计算log2(CPM+1)。C4比例为C4/成熟巨噬。功能端点在固定谱系内计算细胞等权、基因等权log1p(CP10K)平均值；CD8和NK分别计算。Treg比例分母为CD4_Tcon+CD4_Treg，不作每细胞FOXP3阳性筛选。预测量与结局使用完全相同的端点合格文库，患者内文库等权汇总；每患者为一个独立统计单位。

细胞毒集合NKG7、GNLY、GZMB、PRF1、CTSW；功能障碍集合PDCD1、CTLA4、LAG3、HAVCR2、TIGIT、TOX、CXCL13、ENTPD1。活化FOS、JUN、JUNB、NR4A1、CD69、IL2RA单独汇总。

### 关联模型与多重检验

对X/Y患者值分别取平均秩，回归研究分类变量与截距，计算秩残差相关。主要双侧P采用研究内置换，固定99999次、种子262212、加一修正；本次各假设的独特允许置换数均超过枚举门槛。5000次研究内患者分层bootstrap提供效应区间；所有主要端点5000次均有效。未经调整的Spearman及其渐近P另表保存。

12项假设分为预设家族2/1/1/4/2/2：整体SPP1的CD8两功能端点BH；单个Treg组成及次要NK端点不额外校正；两类C4预测量分别在CD8四项、Treg两项和NK两项家族内BH。留一患者、留一研究、单个复杂度协变量、患者集中/高置信派生集合及活化调整仅作预设效应检查，没有增加机会性显著性检验。

### 定向通信与外部对应

使用锁定的人CellChatDB资源中全部9个SPP1受体映射；完整复合体须在同一接收细胞中联合检出，达到至少5细胞及10%比例。受体读出为最低亚基平均log1pCP10K；与SPP1分开计算，不将包含SPP1的乘积分数再与SPP1相关。9受体×3谱系为27项独立登记的受体家族。

NicheNet采用既有human ligand_target_matrix_nsga2r_final资源，按SPP1正权重选择前50靶基因，预先排除全部身份/功能面板、TYMP/SPP1及受体基因；同一固定权重计算先验靶表达。至少40基因在特征索引内、至少10基因在5个以上接收细胞中检出，并有合格完整受体，才进入比较。整体SPP1及C4内SPP1×3谱系共6项，单独BH。小于10患者的通信条目只给描述效应，不计算P；不可估计项P/FDR显示NA，校正内部按原家族大小保留占位。

外部GSE294109按取材/治疗层分开，逐样本使用同样原始计数定义，再在同患者、同时间点和治疗内对合格样本等权汇总。1/2/4/2患者层不计算推断性P/FDR；仅4例层给描述性rho，2例层给原始患者指标与最高/最低SPP1对应的端点差。

## Results

### 1. 整体成熟巨噬SPP1与四个免疫端点均未获支持

研究调整后，SPP1与CD8细胞毒rho=0.01710，P=0.95527，FDR=0.95527；与CD8功能障碍rho=−0.31025，P=0.25593，FDR=0.51186。SPP1与Treg/CD4 rho=−0.11389，P=0.73739；与NK细胞毒rho=0.02193，P=0.95206。后两项为预设单端点，不把未作FDR误标为校正显著。

当前数据不支持把整体成熟巨噬SPP1升高解释为统一的CD8/NK杀伤降低、CD8功能障碍增强或Treg比例增加。相关效应的区间宽且跨零，完整数值如下。

'''+main_table+'''

### 2. C4丰度保留方向性线索，研究调整后的证据尚未成立

C4比例与CD8功能障碍的未经调整rho=0.59706、P=0.01461，为名义相关；研究调整后rho=0.32526、P=0.35629、FDR=0.84644。16次留一患者及4次留一研究均保持正方向，效应范围分别0.0133–0.5751和0.1402–0.3071。该方向性事实保留，但不替代预设主检验。

C4比例与Treg/CD4 rho=0.36645、P=FDR=0.26053；13次留一患者和2次留一研究均正向，同样只作为未获统计支持的候选线索。C4比例与NK细胞毒为正向rho=0.21802，与预期抑制方向相反，P=0.54024、FDR=0.78011，不能隐藏。

C4内SPP1的四个关联均未获支持。其中Treg/CD4为正向rho=0.44136、P=0.23784、FDR=0.26053。该组研究调整信息主要来自GSE207493的9例；总n=10–11不能误写成多研究充分重复。

### 3. 技术与样本影响检查没有改变主结论

SPP1与CD8细胞毒的留一患者效应可跨零；SPP1与NK细胞毒也跨零。SPP1与CD8功能障碍在留一患者时保持负向，但加入接收端检出基因数后变为正向rho=0.20318，显示该方向对技术因素敏感，不能据未经支持的负相关推断保护作用。活化调整后rho=−0.33585，未改变主要结论。

C4比例与CD8功能障碍在去患者集中群、高置信集合、复杂度调整和活化调整后方向仍正，作为一致的方向性线索记录；这些检查没有新的P值，不能据较大的某个调整后rho升级显著性。全部158次留一患者、34次留一研究及36项预设检查已保存，未选择性删除。

### 4. C4的SPP1来源富集保留，但受体/靶表达未形成FDR支持的下游链

在同时具备C4和其他成熟巨噬的固定配对范围，C4−其他巨噬的SPP1伪汇总差值如下。该项是既有来源定位的描述性复核，不是接收端功能验证。

'''+source_table+'''

9个资源受体中，CD44与ITGA4/ITGB1具有达到联合检出门槛的接收端范围；其余整合素复合体未达到本轮门槛，属于不可估计，不称有效阴性或“生物学不存在”。27项登记受体问题中只有4项达到本轮推断资格，另2项为9患者描述，其余未获可估计样本。

CD8的SPP1—ITGA4/ITGB1受体表达容量关联为负向rho=−0.60132，n=12，P=0.02413、FDR=0.65151，只能称名义反向联系。CD4和NK相应9例的描述性rho分别−0.76667和−0.83333，未计算P。CD44对应亦未获支持。所有先验靶表达比较均未通过FDR。

'''+comm_table+'''

### 5. 外部对应仍属小样本描述，没有形成正式重复验证

GSE294109合格范围为9名唯一患者、15个样本、27个患者×端点值；各治疗层始终分开。Ipilimumab/Nivolumab的4例中，SPP1与功能障碍及Treg/CD4描述性rho均为−0.6，未支持预期的正向方向；全部端点和治疗层如下。两例的高低差是具体患者对照，不是估计稳定总体效应。

'''+ext_table+'''

GSE264586未加入正式验证，原因仍是冻结时的身份及肿瘤取材性质缺口；GSE207493不作为独立外部队列重复计数。外部没有用旧Macrophage_SPP1代替冻结C4，没有用旧NK混合群充当常规NK。

## Discussion

本轮完成了预设免疫分支，而不是完成了免疫抑制假设的确认。整体巨噬SPP1相关的主要免疫端点均不支持统一抑制叙事；C4丰度与CD8相关程序/Treg比例的方向线索可以保留，但必须同时保留研究调整后P/FDR及未支持结论。受体可测性和负向名义联系不能填补功能证据缺口。

发现队列16或13例来自不均衡研究，CD8只有3个研究提供研究内变化信息，单例研究不构成独立重复。C4内表达虽然总计10–11例，调整后的变化信息集中在9名同研究患者；外部每层1–4例仅描述。冻结标签仍含明确未定集合，受体复合体联合检出门槛会使稀疏端点不可估计。程序RNA不等于实测杀伤、配体结合或因果作用；NicheNet先验靶表达不是直接SPP1干预响应或签名激活。

前两部分已由数据建立的TYMP相关状态、代谢对应及C4表型不因此被反向否定。历史Phase72采用不同细胞和程序定义，原正负结果完整保留；本次不能表述为独立推翻旧结果，也不能借旧空间阳性改写本轮免疫阴性。当前不再通过换分群、换端点、改阈值或扩展评分追求阳性。

## 验收与交付索引

- I1_patient_inputs：文库指标、158条主要患者假设输入、完整派生集合、原始巨噬SPP1计数、基因覆盖、数据字典和冻结名单匹配。
- I2_main_statistics：12项全部效应、区间、P/FDR、研究内描述、置换/自举分布及独立QR核验。
- I3_robustness：患者/研究影响、36项有限检查、活化对应、三预测量共同C4文库对照。
- I4_directed_communication：资源SHA256、9受体及50靶基因预冻结、每文库联合检出、全部33项登记问题、源区室比较及验证。
- I5_external_description：45个样本×端点指标、27个患者×端点值、12个治疗层端点摘要，无推断性检验。
- EVIDENCE_LEDGER.tsv、FINAL_VALIDATION.json、OUTPUT_SHA256.json：证据账本、整体验收及文件身份。

所有源数据和Phase211冻结文件SHA256未改变；前两部分及C4未覆盖。无新增图、全量CellChat、空间或临床模型。用户离开期间的授权工作I1–I5全部完成，无需以某项显著作为完成条件。
'''
(O/'IMMUNE_RESULTS_CN.md').write_text(report)
# Preserve every registered biological hypothesis in a cross-stage ledger without pooling FDR families.
ledger=[]
for r in main.itertuples():ledger.append(dict(stage='I2',evidence_id=r.hypothesis_id,predictor=r.predictor,endpoint=r.endpoint,n=r.n_patients,effect=r.rho,P=r.p_value,FDR=r.FDR,family=r.family,classification=r.evidence_grade,source=str(O/'I2_main_statistics/PRIMARY_12_HYPOTHESES.tsv')))
for r in lr.itertuples():ledger.append(dict(stage='I4',evidence_id=r.hypothesis_id,predictor=r.predictor,endpoint=r.lineage+':'+r.measure,n=r.n_patients,effect=r.rho,P=r.p_value,FDR=r.FDR,family=r.family,classification='not_estimable' if r.n_patients==0 else ('descriptive_only' if pd.isna(r.p_value) else ('nominal_only' if r.p_value<.05 and r.FDR>=.05 else ('FDR_supported' if r.FDR<.05 else 'not_supported'))),source=str(O/'I4_directed_communication/ALL_DIRECTED_HYPOTHESES.tsv')))
for j,r in enumerate(ext.itertuples()):ledger.append(dict(stage='I5',evidence_id=f'EXT{j+1:02d}',predictor='mature_SPP1',endpoint=r.treatment+':'+r.endpoint,n=r.n_patients,effect=r.rho_descriptive_only,P=np.nan,FDR=np.nan,family='external_description_no_tests',classification='not_assessable_n1' if r.n_patients==1 else 'descriptive_only',source=str(O/'I5_external_description/EXTERNAL_STRATIFIED_DESCRIPTIVE_SUMMARY.tsv')))
pd.DataFrame(ledger).to_csv(O/'EVIDENCE_LEDGER.tsv',sep='\t',index=False)
# Final limited acceptance: required stage completion, complete registered families and immutable freeze.
checks={'all5_stages_complete_validated':True,'frozen_files_sha256_unchanged':True,'all12_main_hypotheses_retained':len(main)==12,'all27_receptor_and6_target_questions_retained':len(lr)==33,'all12_external_strata_endpoints_retained':len(ext)==12,'main_pooled_family_not_used':set(main.groupby('family').size())=={1,2,4},'no_external_inferential_P':ext.inferential_P_not_computed.all(),'no_new_figures':not any(p.suffix.lower() in ['.pdf','.png','.svg','.jpg'] for p in O.rglob('*')),'report_and_complete_ledger_present':(O/'IMMUNE_RESULTS_CN.md').stat().st_size>0 and len(ledger)==57};checks={k:bool(v) for k,v in checks.items()};assert all(checks.values());(O/'FINAL_VALIDATION.json').write_text(json.dumps(checks,indent=2))
now=datetime.datetime.now().astimezone().isoformat();(O/'EXECUTION_STATE.json').write_text(json.dumps({'status':'COMPLETE_VALIDATED_I1_TO_I5','completed_at':now,'stages':{k:v['status'] for k,v in statuses.items()},'report':str(O/'IMMUNE_RESULTS_CN.md'),'primary_conclusion':'Current data do not support a uniform SPP1-associated immunosuppression relationship; C4 directional clues and nominal inverse receptor association retained','no_jobs_remaining':True,'user_notified_during_work':False},indent=2))
backup=O/'final_root_record_backups';backup.mkdir(exist_ok=False)
for n in ['PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','CODEX_HANDOFF.md','AGENTS.md']:shutil.copy2(R/n,backup/n)
note='''\n## 2026-09-19 Phase213：免疫I1-I5顺序执行全部完成\n用户已授权逐步执行、离开期间无需汇报。权威phase213_immune_restart_execution_20260919_v1/IMMUNE_RESULTS_CN.md；I1-I5均COMPLETE_VALIDATED。I1冻结名单与等权文库患者聚合一致；I2完成全部12项、99999次研究内置换及5000次自举；I3保留158个留一患者、34个留一研究和36项有限检查；I4完成9个锁定SPP1受体、50靶基因先验及全部33项问题；I5外部9例/15样本分治疗描述，无P。\n主要结果：成熟巨噬SPP1—CD8细胞毒rho0.0171/P0.95527/q0.95527；—CD8功能障碍rho−0.31025/P0.25593/q0.51186；—Treg/CD4 rho−0.11389/P0.73739；—NK细胞毒rho0.02193/P0.95206。12项调整后P均>0.05。C4比例—CD8功能障碍未经调整rho0.59706/P0.01461名义，调整后rho0.32526/P0.35629/q0.84644；留一患者16/16及研究4/4方向正，保留方向线索，不升级为已证实抑制。C4内表达的研究内变化信息主要来自9例同研究患者，按探索性解释。\nI4唯一名义受体联系为CD8-SPP1/ITGA4_ITGB1负向rho−0.60132 n12/P0.02413/q0.65151；CD4/NK对应n9仅描述；其他稀疏复合体未达资格不是有效阴性。6个先验靶表达均未获FDR支持。外部Ipi/Nivo4例功能障碍/Treg描述rho均−0.6，其他1/2例层不作正式验证。完整正负/NA保留57条证据账本。\n结论：当前数据不支持统一SPP1相关免疫抑制主张；C4表型与方向线索保留，不因下游阴性否定前两部分。不换端点/分群追阳性。全部授权分析完成，无运行任务；未新增图、空间或临床模型，Phase211冻结及C4 SHA256不变。\n'''
for n in ['PROJECT_STATUS.md','DECISIONS.md']:
 with (R/n).open('a') as f:f.write(note)
p=R/'CODEX_HANDOFF.md';p.write_text(note+'\n'+p.read_text())
with (R/'AGENTS.md').open('a') as f:f.write('\n## 2026-09-19 免疫分析I1-I5执行完成\n用户授权按Phase212顺序自主执行，已完成Phase213。主要免疫端点未获研究调整后支持，C4方向线索及负向名义受体结果必须保留。后续不得以旧空间/临床或挑选新阈值替代本轮结果，不覆盖Phase211冻结；新增图或新分析按后续任务推进。当前无运行任务。\n')
# Build artifact identities, excluding mutable session logs and historical root backups.
artifacts=[]
for p in sorted(O.rglob('*')):
 if p.is_file() and 'final_root_record_backups' not in p.parts and p.name not in ['OUTPUT_SHA256.json','finalize.log','handoff_before_execution.md'] and '__pycache__' not in p.parts:
  artifacts.append(dict(path=str(p.relative_to(O)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
(O/'OUTPUT_SHA256.json').write_text(json.dumps(artifacts,indent=2))
with (R/'RESULTS_MANIFEST.tsv').open('a') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n')
 for p in [O/'IMMUNE_RESULTS_CN.md',O/'EVIDENCE_LEDGER.tsv',O/'FINAL_VALIDATION.json',O/'OUTPUT_SHA256.json',O/'EXECUTION_STATE.json']+[O/s/'STATUS.json' for s in stages]+[O/'I2_main_statistics/PRIMARY_12_HYPOTHESES.tsv',O/'I4_directed_communication/ALL_DIRECTED_HYPOTHESES.tsv',O/'I5_external_description/EXTERNAL_STRATIFIED_DESCRIPTIVE_SUMMARY.tsv']:
  w.writerow(['phase213_'+p.parent.name+'_'+p.stem,'phase213','immune_analysis_result',str(p),now,'COMPLETE_VALIDATED',str(O/'FINAL_VALIDATION.json'),'All authorized I1-I5 complete; negative and opposing evidence retained'])
(O/'COMPLETE_I1_I5.txt').write_text('COMPLETE_VALIDATED\n'+now+'\n')
print('ALL_I1_I5_COMPLETE',len(artifacts),'artifacts',len(ledger),'evidence rows');print(source_summary.to_string(index=False));print(json.dumps(checks))
