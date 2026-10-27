from pathlib import Path
import json,datetime,shutil,hashlib
import pandas as pd,numpy as np
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');R=P/'part2_rebuild_v1';O=P/'phase157_ordered_bridge_reassessment_20260915_v1';D=O/'06_report';D.mkdir(exist_ok=False)
v=json.loads((O/'05_validation/STATUS.json').read_text());assert v['status']=='PASS'
def table(df):
 def fmt(x):
  if isinstance(x,float):return 'NA' if np.isnan(x) else f'{x:.6g}'
  return str(x).replace('|','/')
 return '| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(fmt(x) for x in row)+' |' for row in df.itertuples(index=False,name=None))
co=pd.read_csv(O/'01_public_TYMP_MTA/cohort_associations_v2.tsv',sep='\t')
up=co[co.endpoint=='A_B_TYMP_MTA_MARGINAL'][['dataset','n','rho','p','q_within_family']]
mta=pd.read_csv(O/'03_MTA_to_new_C4/MTA_PROXY_C4_BRIDGE_RESULTS.tsv',sep='\t')
state=pd.read_csv(O/'03_MTA_to_new_C4/C4_MTA_RESPONSE_SUMMARY.tsv',sep='\t')
rank=pd.read_csv(O/'04_other_X_screen/CANDIDATE_EVIDENCE_RANKING.tsv',sep='\t')
selected=['program:Cholesterol Homeostasis','program:Glycolysis','ligand:IGFBP3','ligand:CXCL16','ligand:GRN']
focus=rank.set_index('candidate').loc[selected].reset_index()
focus.to_csv(D/'PRIORITY_CANDIDATE_EVIDENCE.tsv',sep='\t',index=False)
brief=focus[['candidate','a_rho','a_P','XY_rho','XY_P','b_partial_rho','b_P']]
interval=focus[['candidate','indirect_ab','ab_low','ab_high','study_ab','study_ab_low','study_ab_high','LOPO_ab_positive','LOSO_ab_positive']]
receiver=focus[focus.family=='secreted_ligand'][['candidate','receptor_label','receptor_genes','n_C4_patients','patients_all_subunits_detected_10pct','median_min_subunit_detection','median_limiting_subunit_C4_minus_other']]
report="""# TYMP—X—新C4按顺序重新考察：实测MTA、MTA相关指标、其他X

日期：2026-09-15。分析根目录：phase157_ordered_bridge_reassessment_20260915_v1。全程以公共资料和冻结项目数据为依据，不借用客户湿实验方向、效应或代谢物变化作为公共结果。

## 三步均完成，已有可继续考察的X候选

1. 公共组织RNA—实测MTA：四队列255份肿瘤样本的TYMP—MTA合并关联接近零，但CPTAC单队列的名义正关联保留。
2. MTA相关指标—新C4：已按新C4重新计算。当前直接MTA响应转录代理在C4中相对其他巨噬降低，肿瘤侧响应代理没有形成通向C4增加的正向桥接；MTAP与C4有原始名义关联，但研究调整后减弱。MTA仍保留为候选，不能用这一转录代理的结果排除真实MTA作用。
3. 其他X：完成116个具备表达覆盖的候选筛查。程序层候选以胆固醇稳态、糖酵解值得继续关注；具体分泌候选以IGFBP3值得进一步核查。CXCL16与C4关联较强，但当前C4受体覆盖缺口明显；GRN有受体相容性，条件桥接关联较弱。

以上是候选发现，不是已确认唯一中介。FDR完整保存但不作结果去留硬门槛；保留方向、效应、患者/研究敏感性和关键反向结果。116个候选中未出现原始模型a、条件b两段同时P<0.05的候选，这不阻止讨论方向一致或某一段已得到支持的线索。

## 第一步：公共TYMP—实测MTA仍有队列差异

{up}

四队列合并rho=0.039467、P=0.538188；MTAP调整后rho=0.044810、P=0.484562。RC18/RC20方向负，CPTAC/CPTAC_val方向正。CPTAC n=50，rho=0.298055、P=0.035526、原四队列家族FDR=0.142106，为名义正向信号，不应隐藏。

同一公共资料中，实测MTA与bulk SPP1的合并rho=0.138937、P=0.029264，原meta多端点FDR=0.097548，四队列方向正。这个结果是组织层SPP1表达，不能改名为C4关联。旧SPP1程序结果也不等同新C4。

本轮复核了原队列表、样本数、既有验证状态，并用Fisher-z重新计算关键合并效应与P值，结果一致。没有从客户实验补入数值，没有把255份bulk组织样本与16例单细胞配对患者合并。组织层TYMP也不自动等同恶性细胞特异TYMP。

## 第二步：MTA相关指标已对新C4重算

现有资料没有同一批样本的实测MTA及新C4标签，因此“实测MTA—C4”的直接关联仍不可估计。本轮完成的是明确标注来源的转录层替代指标，而不是把RNA分数当成MTA浓度。

### C4相对其他巨噬的公共MTA响应相似性较低

复用公共人类成纤维细胞MTA干预的89基因签名，当前输入覆盖84个（68上调、16下调），排除TYMP/SPP1。按研究及细胞区室内逐基因标准化，用既有有符号权重评分，再进行患者内C4—其他成熟巨噬比较：

{state}

肿瘤来源24名患者中19名负向，差中位−0.014108、P=0.000375；全组织26名患者中20名负向。结果说明C4在该评分下不呈更高的公共MTA响应相似性。该签名来自成纤维细胞，不能据此断言MTA不能影响巨噬或新C4；原Phase77跨情境签名资格和迁移边界保持。

### 肿瘤侧MTA相关指标尚未形成正向中间环节

同16名同肿瘤文库患者，以恶性细胞TYMP为A、新C4比例为Y，分别检验四类M：

{mta}

其中a表示TYMP—M关联，XY表示M—C4边际关联，b表示控制TYMP后的M—C4关联。研究调整另列，完整P/FDR和探索性间接项均保留。

肿瘤TYMP与公共MTA响应转录代理的a路径为正（rho=0.617647、P=0.010788），但代理与C4边际关联rho=0.055882，控制TYMP后的方向为负。多胺合成核心ODC1/SRM/SMS及加入AMD1/MAT2A的扩展程序，也未形成同向通向C4增加的桥接。

MTAP表达与C4原始边际rho=0.600、P=0.014007，控制TYMP后rho=0.540809、P=0.037376；研究调整后该关系接近零。它是清除相关基因的转录关联，不能写成“MTAP缺失引起MTA积累并促进C4”。

MTA响应代理原始模型的间接项为负，bootstrap诊断区间为−0.560992至−0.031113；研究调整后区间跨零。即使保留这条负向信号，也不能将其改写为解释TYMP—C4正向关系的正中介。

## 第三步：其他X的候选发现

公共候选资源为冻结CellChat表的3,233个交互记录。主要分泌轨道限定为资源中Secreted Signaling类、可直接映射到单基因的配体；跨复合物或非蛋白标签不按单基因配体强行解释。候选基因需在至少12/16个恶性细胞文库中有至少10个伪汇总counts；程序至少保留10个合格基因，并排除TYMP/SPP1。

最终纳入56个该类配体基因、50个Hallmark程序、10个代谢/转运基因，共116个X。代谢列表中的SLC29A2、CHKA、PCYT2未达覆盖条件，记录为当前不可充分评估，不作为阴性检验。完整候选、排除原因和全部结果保存，未仅保存排名靠前者。

### 优先关注的候选及两端证据

{brief}

所有行n=16。a为TYMP—X；XY为X—C4边际关系；b为控制TYMP后的X—C4条件关系。上表P均为原始P；相应家族FDR保存在完整候选表。不能把XY的显著性直接转给b或间接效应。

**胆固醇稳态程序**：TYMP—程序和程序—C4两个边际端点均有名义正关联。控制TYMP后b仍正但P=0.106742。原始与研究调整的间接项均正，16次留一患者乘积均正；去掉GSE207493后剩6名患者，乘积接近零且略负。因此可作为程序层候选X，保留研究依赖。该指标来自恶性细胞Hallmark Cholesterol Homeostasis，不与此前巨噬C4中的KEGG Cholesterol metabolism混称同一结果，也不是胆固醇实测浓度。

**糖酵解程序**：与C4边际rho=0.564706、P=0.022663；TYMP—程序及条件b均为正向趋势。16次留一患者乘积为正。虽然4次留一研究乘积均为正，但剔除GSE207493后a和b均变为负、乘积接近零；不能仅凭乘积符号称正向路径跨研究稳健。

**IGFBP3**：与C4的边际rho=0.714706、P=0.001863；控制TYMP后rho=0.657786、P=0.007690，配体家族b路径FDR=0.107655，属于名义关联。TYMP—IGFBP3为正向趋势（rho=0.376471、P=0.150648）。原始模型16次留一患者和4次留一研究的间接项均为正，剔除任一研究的两段方向也均正；研究调整后条件关系减弱。它适合作为下一轮具体分泌分子候选，但尚未形成两端均获充分支持的中介。

**CXCL16**：与C4边际rho=0.747059、P=0.000882；控制TYMP后rho=0.724335、P=0.002258，配体家族FDR=0.063231；研究调整后的条件关系也保留（rho=0.722959、P=0.007892）。不过TYMP—CXCL16一端较弱，且当前C4的CXCR6表达覆盖不足。它首先是C4相关的候选生态指标，不能直接写成CXCL16对C4的已确认作用。

**GRN**：有公开受体配对和C4接收端表达相容性；TYMP一端方向为正，研究调整后名义支持，但原始条件X—C4关系较弱。作为辅助候选保留，不优先宣称中介。

### 接收端相容性帮助区分候选角色

{receiver}

上述接收端检查来自24名具备足量C4及其他巨噬的肿瘤患者，不是额外独立验证队列。门槛是每个受体复合物全部所需亚基在C4中均达到10%细胞检出；它用于覆盖判断，不等同生物学作用阈值。

IGFBP3的公共配对受体TMEM219在24/24患者达到该覆盖标准，C4内检出比例中位约84.17%，但相对其他巨噬并不富集。GRN的SORT1在22/24患者达到覆盖标准，C4相对其他巨噬的限制亚基表达差中位约+0.483。CXCL16对应CXCR6为0/24达到门槛，检出中位数为0；这是本轮直接接收端证据不足，不是证明任何细胞或条件下都无该作用。

### 间接项及其不确定性

{interval}

a×b为标准化秩线性模型的路径乘积；区间为固定秩资料、研究内患者重抽样2,000次的探索性bootstrap区间。它不是多重校正后的中介显著性，也不能识别因果。原始模型胆固醇稳态和糖酵解的区间位于正侧，研究调整区间均跨零；其余候选区间亦需保留。不计算或宣传“解释了多少百分比的中介效应”。

## 当前优先级

若下一步需要具体分子X，优先核查IGFBP3，GRN作辅助；CXCL16保留为关联较强但接收端待解释的候选。
若允许X是肿瘤状态程序，优先保留胆固醇稳态与糖酵解，进一步区分研究背景与可复现的肿瘤—C4联系。
MTA仍留在候选图谱中，但当前这些公共转录指标未形成解释TYMP—C4正向趋势的证据链。后续若有同样本实测MTA与新C4信息，属于有新增判别力的检验；不能把现有转录代理当成已完成这项实测检验。

## 统计方法与执行记录

所有跨细胞区室候选关联使用冻结16患者/16肿瘤文库：GSE156632 4、GSE159115 1、GSE178481 1、GSE207493 10。癌细胞资格、C4成员和成熟巨噬分母均未改变。共提取34,394个恶性/成熟巨噬细胞，形成154个文库×群体伪汇总；原始文库counts逐块校验。细胞数、签名基因数和检查数均不作为统计样本量。

a、XY和条件b均基于秩转换；每候选单独建模，不将116个X同时拟合。原始及研究调整分开；MTA四指标与其他三类候选分别在各家族/检验边/调整方式内BH。FDR不作为继续分析或正文去留硬门槛。研究内两个单患者队列不能提供各自研究内变异；去主要研究后的小样本仅解释方向。

独立复核{checks}项全部通过：关键MTA meta重算、冻结C4身份、表达指标与Phase154一致性、MTA签名覆盖与方向、患者配对差、所有候选的独立秩回归/偏相关/P值/线性路径分解和BH。验证不是外部生物学复现。

数据来源：SOURCE_MANIFEST.json及01_public_TYMP_MTA、02_inputs；MTA新C4结果：03_MTA_to_new_C4；全部其他候选：04_other_X_screen；核查：05_validation；本报告和优先表：06_report。全部脚本位于当前分析根目录。

本轮未调用客户湿实验、未重训SCENIC/CellChat/Cell2location、未生成图件。当前完成的是按顺序重新考察和候选发现；尚未完成独立队列/干预验证。
"""
report=report.replace('{up}',table(up)).replace('{state}',table(state)).replace('{mta}',table(mta[['metric','adjustment','a_rho','a_P','XY_rho','XY_P','b_partial_rho','b_P']])).replace('{brief}',table(brief)).replace('{receiver}',table(receiver)).replace('{interval}',table(interval)).replace('{checks}',str(v['checks']))
(D/'ORDERED_TYMP_X_C4_REASSESSMENT_CN.md').write_text(report,encoding='utf-8')
summary="""# 按顺序重新考察的简明结果

第一步：四公共配对队列255份肿瘤样本，TYMP—实测MTA合并rho0.039、P0.538；CPTAC单队列rho0.298、P0.0355为保留的名义正向信号。该组织层数据不能直接代表新C4。

第二步：已对新C4重算MTA相关指标。公共MTA干预响应转录代理在C4相对其他巨噬降低（肿瘤24患者19负，P0.000375）。恶性TYMP与该响应代理正相关，但代理未形成通向C4增加的正向桥。MTAP与C4原始有名义关联，研究调整后减弱。没有同样本实测MTA+C4，所以不能排除真实MTA作用，也不能把转录分数当MTA浓度。

第三步：完成116个候选X筛查（56公共配体、50程序、10代谢/转运基因）。程序层重点是胆固醇稳态、糖酵解；具体分子优先候选是IGFBP3，GRN作为辅助。CXCL16与C4关联强，但C4的CXCR6接收端覆盖不足。

IGFBP3—C4 rho0.715、P0.00186，控制TYMP后rho0.658、P0.00769；TYMP—IGFBP3仅正向趋势rho0.376、P0.151。公共配对受体TMEM219在C4中有表达，尚不能称已验证中介。
胆固醇稳态程序两端边际关联均为正且名义显著：TYMP—程序rho0.591、P0.0159；程序—C4 rho0.585、P0.0172；条件b和研究敏感性需保留。

所有候选效应、P/FDR、相反方向、留一结果和接收端依据完整归档。FDR不作为去留硬门槛；本轮得到候选线索，没有确定唯一中介。独立复核1246项通过，无新图，无客户湿实验借用。
"""
(D/'SUMMARY_CN.md').write_text(summary,encoding='utf-8')
B=O/'state_before_update';B.mkdir(exist_ok=False)
for p in [P/'CODEX_HANDOFF.md',P/'PROJECT_STATUS.md',P/'DECISIONS.md',P/'RESULTS_MANIFEST.tsv',R/'STATUS.json']:shutil.copy2(p,B/p.name)
note="""## 2026-09-15 Phase157 按顺序完成TYMP—MTA、MTA—新C4、其他X重新考察

用户指定顺序已执行。权威报告phase157_ordered_bridge_reassessment_20260915_v1/06_report/ORDERED_TYMP_X_C4_REASSESSMENT_CN.md；独立复核1246/1246 PASS，无新图、无客户湿实验使用。
1 公共组织RNA—实测MTA：4队列255肿瘤样本合并rho0.03947/P0.5382，保留CPTAC n50名义正向rho0.2981/P0.0355。实测样本没有新C4，不代替下游。
2 新C4：公共成纤维MTA响应签名84/89覆盖，C4相对其他巨噬肿瘤n24中19负，P0.000375。肿瘤TYMP与响应RNA代理正相关，但代理未形成正向C4桥；MTAP原始C4 rho0.6/P0.014，研究调整后弱。保留MTA候选；无同样本实测MTA+C4，不能排除真实MTA机制。
3 其他X：筛查116项（56公共单基因配体、50Hallmark程序、10代谢/转运）。程序层候选胆固醇稳态/糖酵解；具体分子优先IGFBP3（X-C4 rho0.7147/P0.00186，条件b rho0.6578/P0.00769；A-X弱趋势rho0.3765/P0.1506；C4 TMEM219覆盖24/24）。GRN接收端SORT1相容但b弱。CXCL16-C4关联强，C4 CXCR6未达覆盖标准，不写直接作用。
116项无a与条件b同时名义显著者；FDR不是行动门槛，候选按方向/稳定性/作用依据分级。胆固醇和糖酵解受GSE207493影响，原始bootstrap正区间不能升级因果中介。全部负向、FDR、留一及表达排除均保留。
下一阶段可围绕IGFBP3及肿瘤胆固醇/糖酵解候选做独立资料和机制对应核查；不自动宣称找到了唯一X。长期SSH已建立用于本轮，操作记录和脚本全部服务器归档。
"""
for fn in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 p=P/fn;p.write_text(note+'\n'+p.read_text(),encoding='utf-8')
p=R/'STATUS.json';d=json.loads(p.read_text());d['remaining_analysis'].update(current_step='ordered_bridge_reassessment_complete',latest_report=str(D/'ORDERED_TYMP_X_C4_REASSESSMENT_CN.md'),MTA_retained_as_candidate=True,confirmed_mediator=False);p.write_text(json.dumps(d,ensure_ascii=False,indent=2))
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
with (P/'RESULTS_MANIFEST.tsv').open('a') as f:
 for path in [D/'ORDERED_TYMP_X_C4_REASSESSMENT_CN.md',D/'PRIORITY_CANDIDATE_EVIDENCE.tsv',O/'03_MTA_to_new_C4/MTA_PROXY_C4_BRIDGE_RESULTS.tsv',O/'04_other_X_screen/ALL_CANDIDATE_BRIDGE_RESULTS.tsv',O/'05_validation/STATUS.json']:
  f.write('\t'.join(['phase157_'+path.stem,'phase157','ordered_reassessment',str(path),now,'COMPLETE_REVIEWED','Exploratory evidence; no FDR action gate','No figures; no client data'])+'\n')
status={'status':'COMPLETE_ORDERED_PUBLIC_BRIDGE_REASSESSMENT','steps_complete':3,'validation_checks':v['checks'],'MTA_retained':True,'confirmed_mediator':False,'candidate_count':116,'figures_created':False,'client_data_used':False,'time_utc':now}
(O/'STATUS.json').write_text(json.dumps(status,indent=2))
receipt=[]
for p in [D/'ORDERED_TYMP_X_C4_REASSESSMENT_CN.md',D/'SUMMARY_CN.md']:
 receipt.append(dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
(D/'DELIVERY_MANIFEST.json').write_text(json.dumps(receipt,indent=2))
print(json.dumps(status));print(json.dumps(receipt))

