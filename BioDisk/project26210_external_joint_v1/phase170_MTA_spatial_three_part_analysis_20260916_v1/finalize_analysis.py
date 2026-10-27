from pathlib import Path
import json,hashlib,datetime,csv,shutil
import numpy as np,pandas as pd
O=Path(__file__).resolve().parent;R=O.parent
s=pd.read_csv(O/'PATIENT_SUMMARY.tsv',sep='\t');d=pd.read_csv(O/'SECTION_EFFECTS.tsv',sep='\t');v=json.loads((O/'VALIDATION_STATUS.json').read_text())
assert v['checks']==v['passed']
main=s[(s.model=='extended_primary')&(s.stratum=='all')&(s.policy=='equal_fisher_mean')].copy()
def fmt(x):return 'NA' if pd.isna(x) else f'{x:.6g}'
labels={'source':'源区TYMP—','receiver_link':'源区候选—邻域SPP1：'}
lines=['| 关系 | 独立患者n | 正/负 | 患者rho中位数 | 中位数95%区间 | P | FDR |','|---|---:|---:|---:|---|---:|---:|']
for x in main.sort_values(['role','candidate']).itertuples(index=False):
 lines.append(f'| {labels[x.role]}{x.candidate} | {x.n_patients} | {x.positive}/{x.negative} | {fmt(x.median_rho)} | [{fmt(x.median_ci_low)}, {fmt(x.median_ci_high)}] | {fmt(x.P)} | {fmt(x.FDR)} |')
table='\n'.join(lines)
report=f"""# 三部分空间分析完成：多胺桥接仅有弱方向线索，SLC16A3上游联系较一致

## 总结
当前24张GSE175540切片的患者层分析显示：TYMP—SRM/SMS上游与SRM/SMS—邻域SPP1均呈弱正向总体趋势，但尚未获得跨患者一致的统计支持。接收端信号在FFPE/BIONIKK较强，冷冻或去除BIONIKK后不能稳定复现，并在坐标调整后明显减弱。本轮没有为TYMP—多胺/MTA—邻域SPP1完整空间链提供稳健支持。
TYMP—SLC16A3上游呈小效应、较一致正向关联，20/23患者为正，患者相关中位数0.03890，P=0.000488、FDR=0.001465；其接收端未得到支持。MTAP RNA背景调整对上述效应影响很小。
以上是本轮空间数据的结论；此前第一部分SRM程序、公共SMS—实测MTA及SMS—组织SPP1的分段事实继续保留，不被本轮覆盖或删除。

## Results

### 1. TYMP与SRM/SMS仅呈弱正向空间联系
SRM和SMS均为15/23患者正向；中位partial rho分别0.008741和0.019980，P均0.21004、FDR均0.21004。两者患者中位数bootstrap区间均跨0。留一患者后总体中位方向仍正，但方向稳定与达到统计支持分开表述。
SMS在FFPE为9/12正、中位0.02121，在frozen为7/12正、中位0.004261；SRM分别7/12、8/12正。主线定位为弱方向线索，不称名义显著或校正后显著。

### 2. SRM/SMS与邻域SPP1总体偏正，跨制备及空间调整后支持不足
SRM接收端11/19正，中位0.02726，P=0.64761、FDR=0.64761；SMS为13/19正，中位0.02495，P=0.16707、FDR=0.50121。主检验未达到名义显著。
探索性FFPE子集（9患者）：SRM7/9正、中位0.06621；SMS8/9正、中位0.03990。该小样本子集只报告方向和效应，不计算推断性P/FDR。frozen子集（11患者）：SRM4/11正、中位-0.003720；SMS6/11正、中位0.002664。不能写成跨制备一致的巨噬接收支持。
FFPE全部来自BIONIKK，制备与队列来源部分重合，不能将差异单独归因制备技术。去BIONIKK后仅剩8患者：SRM3/8正、中位-0.009337；SMS4/8正、中位-0.002072，作为主要相反证据保留。
坐标扩展模型因GSM5924041条件数46.9–47.4超过30而不拟合该片接收端。相同18患者配对比较：SRM主模型中位0.02119，坐标调整后0.002573；SMS由0.01994降至0.003502。调整后均10/18正，P/FDR=0.81453。该衰减在相同患者集合可见，并非单纯样本数变化。坐标调整结果提示现有弱邻域关联对较大尺度空间结构敏感，不能据此宣称空间特异的跨细胞作用。

### 3. MTAP调整没有明显改变效应；SLC16A3上游较一致但接收端不支持
MTAP背景端点：TYMP—MTAP为14/23正，中位0.008034，P=0.40487、FDR=0.64761；源区MTAP—邻域SPP1为11/19正，中位0.007208，P/FDR=0.64761。当前空间RNA数据不支持一致的MTAP关联。
加MTAP调整后：上游SRM/SMS中位0.008299/0.01966，接收端0.02639/0.02446。配对患者rho变化中位数绝对值均<0.001。保留MTAP作为清除背景，但本轮未显示它解释了主要空间关系。
SLC16A3上游20/23正，中位0.03890，95%患者bootstrap区间[0.03079,0.05792]，P=0.000488、FDR=0.001465。基础组成、MTAP、坐标与最大连通区敏感性均保留20/23正；坐标调整中位0.03177。c_2仅FFPE、仅frozen、完全排除的FDR分别0.001465、0.007798、0.002566。
FFPE及frozen上游各10/12正，分别中位0.03775/0.04184，均P=0.03857、FDR=0.11572，仅名义显著；两组共享c_2，不能当两个完全独立验证队列。
SLC16A3接收端仅8/19正，中位-0.01398，P/FDR=0.64761；留一患者后中位数始终为负（-0.01497至-0.008227）。所以可保留上游RNA状态联系，当前数据不支持其作为已串联至邻域SPP1的空间桥梁，也不用于替换MTA假说。

### 主模型完整统计
{table}

rho为先逐片计算、再按患者合并后的partial Spearman效应；表中中位数不是把所有spot池化得到的相关。P为患者方向一致性的双侧精确符号检验；FDR在预设上游3候选、接收端3候选、MTAP背景2端点三个家族内分别计算。区间为独立患者bootstrap的中位数区间。

## Methods
输入使用Phase165全部24片原始counts/坐标、Phase78冻结q05丰度、Phase167的源区/接收端成员及图结构，患者映射为Phase168。每片表达按log1p(counts/总counts×10000)计算，18目标/背景基因保留，float64计算。
源区为恶性归一化q05≥20%且占优；接收区合并成熟巨噬q05≥10%且不低于片内中位数。邻接最多6个非自身第一圈邻居，距离≤最近邻中位距离×1.3；每个源spot至少3个合格接收邻居。技术门槛为50源spot、满列秩、标准化设计条件数<30、有限且非零残差方差。
主模型先对每个端点区域的变量分别取平均秩，再控制深度、恶性与合并巨噬组成、基质与缺氧评分；接收端同时控制相应邻域变量及源区TYMP。基质和缺氧评分沿用原冻结基因，不含TYMP/SPP1/SRM/SMS/SLC16A3。分别残差化预测变量与结局后计算相关。
每张切片贡献一个效应；c_2两种制备用Fisher z等权平均后反变换，只计一位患者。其余患者各一片，所有患者等权；不按spot量加权。上游24片/23患者，接收端20片/19患者。患者中位数、四分位数、正负方向、10000次患者bootstrap及双侧精确符号检验使用预冻结方案与随机种子20260916。n<10仅描述；队列分层与留一队列均按预计划只描述。完整P/FDR及不可估计端点保留。
敏感性包括仅组成、MTAP背景、坐标x/y/x²/xy/y²、最大连通区、GSM5924048两个合格分区、分制备、分队列、留一队列、留一患者以及三种c_2处理替代方案；各敏感性P/FDR单独保留，不替换主模型统计。不拟合任何spot层P值或中介模型。

## Discussion与证据边界
本轮支持的是条件化空间RNA关联。MTA没有被直接测量；SRM/SMS表达不能换算MTA浓度或净生成，MTAP RNA不能代表酶活性或清除通量。这里也未测试生物学因果中介。
源区和接收区为组成富集区域，spot并非纯肿瘤细胞或纯巨噬细胞；邻域SPP1不能完全归属巨噬来源。旧13状态参考包含目标基因，目标基因剔除后重新解卷积未在本次范围内实施，丰度调整及区域选择仍条件于该参考。旧SPP1-state未改名新C4；整体巨噬SPP1联系与C4单细胞定位分工保留。
患者层检验避免将spot量当生物学n，但不能排除共享组织梯度制造的区域共表达；坐标敏感性用于检查该问题。GSM5924048的原分离区域未连接，分区结果及相反方向完整保留。全部样本经原研究TLS选择，制备与队列部分混杂，逐例治疗状态未确认；不据此作治疗反应或ICI耐药推断。
因此，MTA多胺轴继续保留已有分段证据及客户实验背景；这轮空间数据只补充有限方向线索，不能宣称完整空间证据链得到确认。较一致的SLC16A3上游联系是并行发现，下游缺口仍明确保留。

## 补充材料与可复核性
- SECTION_EFFECTS.tsv：780个切片/端点/模型组合，713可估计，所有失败原因保留。
- PATIENT_EFFECTS.tsv / PATIENT_SUMMARY.tsv：患者效应和132行主/敏感性总结。
- LEAVE_ONE_PATIENT.tsv：逐患者删除稳定性。
- MTAP_ADJUSTMENT_PAIRED.tsv / COORDINATE_ADJUSTMENT_PAIRED.tsv：相同患者的调整前后配对效应。
- design_cache/：每个模型实际使用的秩设计、结局和spot索引，便于独立复核。
- INPUT_MANIFEST.tsv：101个输入/冻结文件的执行前哈希；执行后全部一致。
- FROZEN_ANALYSIS_PLAN.md / PLAN_FREEZE.json：查看关联结果前冻结的方案和时间。
- VALIDATION.tsv / QR_RECOMPUTED_EFFECTS.tsv：{v['checks']}/{v['passed']}检查通过；全部713可估计相关经独立QR投影复算，最大误差{v['max_qr_error']:.3g}；独立组合公式复核患者精确符号P值，SciPy复核BH。
- 软件：Python3.11.15、numpy1.26.4、scipy1.17.1、pandas2.3.3、anndata0.12.16、h5py3.16.0。执行OPENBLAS_NUM_THREADS=1、OMP_NUM_THREADS=1。
新分析图0；客户原始实验数据入模0。既有阶段均只读；所有进度记录在服务器。
"""
(O/'THREE_PART_SPATIAL_RESULTS_CN.md').write_text(report,encoding='utf-8')
main.to_csv(O/'PRIMARY_RESULTS.tsv',sep='\t',index=False)
now=datetime.datetime.now().astimezone().isoformat()
status=json.loads((O/'COMPUTE_STATUS.json').read_text())
status.update(status='PASS_THREE_PART_ANALYSIS_COMPLETE',validation=v,independent_source_patients=23,independent_receiver_patients=19,primary_findings=dict(MTA_polyamine_complete_spatial_chain='NOT_SUPPORTED_ROBUSTLY',SRM_SMS='WEAK_POSITIVE_DIRECTION_HETEROGENEOUS',MTAP_adjustment='MINIMAL_RNA_EFFECT_CHANGE',SLC16A3_source='SMALL_EFFECT_DIRECTIONALLY_CONSISTENT_FDR_SIGNIFICANT',SLC16A3_receiver='NOT_SUPPORTED'),completed_at=now)
(O/'STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
(O/'EXECUTION.json').write_text(json.dumps(dict(command='OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /mnt/BioDisk/wangmingyuan/project26210_envs/part4_advanced_py_v1/bin/python -u '+str(O/'run_analysis.py'),validation='QR projection + independent binomial formula + scipy BH + file hashes',association_plan='FROZEN_ANALYSIS_PLAN.md',failures='Intermittent SSH connection closures retried; analysis and validation completed without computational errors'),indent=2))
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
pd.DataFrame([dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(O.rglob('*')) if p.is_file() and 'authority_backups' not in p.parts and p.name not in ['ARTIFACT_HASHES.tsv','FINAL_STATE_VALIDATION.json']]).to_csv(O/'ARTIFACT_HASHES.tsv',sep='\t',index=False)
B=O/'authority_backups';B.mkdir(exist_ok=False)
statefiles=['PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','CODEX_HANDOFF.md']
for n in statefiles:shutil.copy2(R/n,B/n)
notice=f"""
## Phase170 MTA空间三部分分析完成（{now}）
用户授权三部分：TYMP—SRM/SMS上游、邻域SPP1接收端、MTAP背景/SLC16A3并行。方案先冻结，c_2两制备Fisher等权合并；主模型扩展partial Spearman，患者符号检验和三家族BH，无spot推断/图。
主结果：SRM/SMS上游均15/23正，rho中位0.008741/0.01998，P/FDR0.21004；接收端SRM11/19正rho0.02726/P0.64761/FDR0.64761，SMS13/19正rho0.02495/P0.16707/FDR0.50121。FFPE/BIONIKK较强、frozen和去BIONIKK不稳，坐标调整同18患者效应趋近0。当前不支持稳健的完整MTA空间链，既有分段事实保留。
MTAP加入前后配对rho改变中位绝对值<0.001；MTAP背景两端均不支持一致关联。SLC16A3上游20/23正，rho0.03890/P0.000488/FDR0.001465，敏感性保留；接收端8/19正rho-0.01398/P=FDR0.64761，不能作已完成的替代桥。
报告：{O/'THREE_PART_SPATIAL_RESULTS_CN.md'}
713效应QR复算通过，6626检查全PASS，101输入哈希不变。患者n23/19，治疗未知保留；无MTA实测、无中介、无新C4投射、无目标基因剔除解卷积。三部分已完成，待用户审核结果，不自动追加模型或画图。
"""
for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']:
 p=R/n;t=p.with_name(n+'.phase170.tmp');t.write_text(notice+p.read_text(),encoding='utf-8');t.replace(p)
p=R/'RESULTS_MANIFEST.tsv'
with p.open(newline='') as f:rows=list(csv.reader(f,delimiter='\t'))
assert len(rows[0])==8 and not any(r and r[0].startswith('phase170_') for r in rows[1:])
reg=['THREE_PART_SPATIAL_RESULTS_CN.md','PRIMARY_RESULTS.tsv','SECTION_EFFECTS.tsv','PATIENT_EFFECTS.tsv','PATIENT_SUMMARY.tsv','LEAVE_ONE_PATIENT.tsv','MTAP_ADJUSTMENT_PAIRED.tsv','COORDINATE_ADJUSTMENT_PAIRED.tsv','FROZEN_ANALYSIS_PLAN.md','PLAN_FREEZE.json','INPUT_MANIFEST.tsv','INPUT_VALIDATION.tsv','QR_RECOMPUTED_EFFECTS.tsv','VALIDATION.tsv','VALIDATION_STATUS.json','STATUS.json','EXECUTION.json','run_analysis.py','validate_analysis.py','summarize_sensitivity.py','finalize_analysis.py','ARTIFACT_HASHES.tsv']
for n in reg:
 path=O/n;assert path.exists()
 rows.append(['phase170_'+n,'Phase170','MTA_spatial_three_part',str(path),now,'PASS','23/19 independent patients; 6626 validation checks pass','sha256='+sha(path)])
t=p.with_name(p.name+'.phase170.tmp')
with t.open('w',newline='') as f:csv.writer(f,delimiter='\t',lineterminator='\n').writerows(rows)
t.replace(p)
gate=dict(stage='Phase170',four_backups=all((B/n).exists() for n in statefiles),registered_artifacts=len(reg),registered_paths_exist=all((O/n).exists() for n in reg),latest_state_updated=all((R/n).read_text().startswith(notice) for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']),authority_sha256={n:sha(R/n) for n in statefiles})
assert gate['four_backups'] and gate['registered_paths_exist'] and gate['latest_state_updated']
(O/'FINAL_STATE_VALIDATION.json').write_text(json.dumps(gate,indent=2))
print(json.dumps(status,ensure_ascii=False,indent=2));print(json.dumps(gate,indent=2))
