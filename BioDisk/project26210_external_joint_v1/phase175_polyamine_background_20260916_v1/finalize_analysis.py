from pathlib import Path
import json,csv,hashlib,datetime,shutil
import pandas as pd
O=Path(__file__).resolve().parent;R=O.parent
S=pd.read_csv(O/'MALIGNANT_PROGRAM_SUMMARY.tsv',sep='\t');E=pd.read_csv(O/'PUBLIC_PROGRAM_EDGES_META.tsv',sep='\t');M=pd.read_csv(O/'PUBLIC_PROGRAM_SHORT_META.tsv',sep='\t');V=json.loads((O/'VALIDATION_STATUS.json').read_text());assert V['checks']==V['passed']
def f(v):return '-' if pd.isna(v) else f'{v:.6g}'
tables=[]
for title,data,cols in [
('恶性细胞全部结果',S,['module','score','model','n','positive','median_rho','P','FDR2','low','high']),
('公共条件边全部结果',E,['module','adjustment','edge','n','positive_cohorts','rho','P','FDR2','historical_FDR2','historical_FDR3','I2']),
('公共短路径全部结果',M,['module','adjustment','n','positive_cohorts','estimate','low','high','P','FDR2','bootstrap_low','bootstrap_high'])]:
 lines=['## '+title,'','| '+' | '.join(cols)+' |','|'+'---|'*len(cols)]
 for _,r in data.iterrows():lines.append('| '+' | '.join(f(r[c]) if isinstance(r[c],(int,float)) else str(r[c]) for c in cols)+' |')
 tables.append('\n'.join(lines))
(O/'ALL_STATISTICS_CN.md').write_text('\n\n'.join(tables),encoding='utf-8')
report="""# 多胺代谢背景补充：广义代谢状态获得分段支持，MTA短路径仍未确认

## 总结
采用早已定义的两个模块，未搜索新的基因组合：合成核心ODC1/SRM/SMS；SAM相关供给与脱羧程序AMD1/MAT2A（源列名SAM_supply）。MTAP保持清除背景，不加到评分；TYMP/SPP1不参与评分。
本轮支持将“TYMP伴随蛋氨酸/SAM及部分多胺相关状态”作为正文代谢背景。较清楚的新支持位于恶性细胞SAM相关程序，以及公共组织TYMP—合成核心的关联。恶性细胞合成核心调整后不稳定，核心—实测MTA主要模型较弱，两个模块的合并短路径均未确认。
因此，该结果可以补充MTA候选选择的代谢理由，不能升级为已证明“TYMP通过多胺程序使MTA积累”。分段事实保留；本轮不重新组合模块寻找阳性。

## 结果一：恶性细胞的SAM相关程序保留小效应正向联系
16位患者中，主要基因z均值评分经深度、增殖和缺氧调整后12/16正，中位rho=0.01810，P=0.02496、两模块FDR=0.04993，患者bootstrap中位数区间[0.006961,0.03448]。该FDR靠近0.05，应与效应大小和敏感性一起呈现，不称强关联。
加入检出基因数后12/16正，中位rho=0.01551，P=0.03864、FDR=0.07727，只保留名义显著。预设平均秩评分敏感性在主调整为14/16正、rho=0.02190、P=0.01099、FDR=0.02197；进一步加检出基因数为13/16正、rho=0.01745、P=0.02139、FDR=0.04279。两种评分分别保留，不将敏感性替换主结果。
各研究效应有差异：GSE1566324/4正，GSE2074937/10正；另两个研究各1例，一正一负。去GSE207493后6例5正、中位0.01567，属于小样本描述。
该程序代表MAT2A/AMD1的RNA状态，并不等同实测SAM供给量、dcSAM浓度或酶活性。

## 结果二：合成核心原始关联较清楚，但条件关联减弱
14位合格恶性细胞患者中，ODC1/SRM/SMS基因z均值模块未调整为12/14正、中位rho=0.07083，P=0.00305、FDR=0.00305。仅深度调整为9/14正、中位0.01147，P/FDR=0.46313；主调整为8/14正、中位0.01244，P/FDR=0.62573。平均秩评分主调整同为8/14正、中位0.01339，P/FDR=0.35754。
同14患者直接对照，主模型SMS中位rho=0.001869、合成核心0.01244，患者配对差的中位数0.01670；这是描述性对照，不能声称两者差异显著或合成核心已获稳定支持。
逐基因删除仅作构成诊断：主调整下去SMS中位0.03542、11/14正；去ODC1中位0.005598、9/14正；去SRM中位-0.01067、7/14正。不能据此把去SMS的两基因组合改成新主模块，亦不新增该组合P值。
研究差异需保留：GSE156632三例均负，GSE207493十例七正；去GSE207493只剩四例，一正三负，中位-0.01555。原始关联与调整后关联回答不同问题，本轮不将后者弱化解释成原始状态事实不存在。

## 结果三：公共组织TYMP—合成核心得到支持，至MTA仍较弱
四个配对RNA—实测代谢组队列共255肿瘤患者，复用既有模块列和样本。
主要MTAP+髓系调整下：
- TYMP—合成核心：3/4队列正，随机效应rho=0.22792，P=0.003844、本轮两模块FDR=0.007689；历史三候选FDR=0.005767并列保留，I²=31.69%。
- 合成核心—MTA（再控制TYMP）：3/4队列正，rho=0.13889，P=0.08646、原/本轮两模块FDR=0.12729，未达名义显著。
加缺氧后，上游减弱至rho=0.15158、P=0.12060、FDR=0.24120；下游4/4同向，rho=0.16245、P=0.04679、FDR=0.09358，为名义正向。不能将主模型较好的上游与另一敏感性模型较好的下游拼成一条同时成立的链。
既有髓系调整、未控制TYMP/MTAP的核心—MTA边际关联rho=0.19821、P=0.01589、历史两模块FDR=0.03178保留为背景。它与主要条件边是不同问题，不因条件边较弱删除，也不替代该条件边。

## 结果四：两个模块均未形成稳定的MTA短路径
主要模型中，合成核心的合并a×b=0.01471，随机效应95%区间[-0.02964,0.05906]，P=0.51557、两模块FDR=0.85427；固定汇总权重的患者bootstrap区间[-0.03110,0.05446]。
SAM相关程序的合并a×b=-0.006069，区间[-0.07083,0.05869]，P/FDR=0.85427；bootstrap区间[-0.06296,0.04474]。其公共TYMP上游及MTA下游主要模型也均未达名义显著。
全部8个模块×调整合并模型的正态近似和bootstrap区间均跨零；合成核心主要模型留一队列亦均跨零。当前不支持将任何一个模块写成已经确认的TYMP—MTA中介。

## 对统一文稿的建议
可以写：
“TYMP相关恶性上皮状态伴随蛋氨酸及部分多胺相关程序；补充分析显示，SAM供给与脱羧相关RNA程序保留小效应正向联系，公共组织数据也支持TYMP与预定义多胺合成核心的关联。这些分段结果为将MTA作为相关候选代谢物提供了代谢背景。”

同时保留关键边界：
“多胺核心在恶性细胞中的调整后关联、其与实测MTA的条件关联及短路径分解尚不稳定，因此不能据此确认TYMP通过该程序引起MTA积累。”
原第一部分已经建立的通路结果及客户实验继续保留，本轮不将转录评分当成MTA测量，不新增胞外释放、巨噬响应或空间屏障结论。下一步若推进主线，重点应转向MTA与髓系/SPP1状态的具体依据，避免继续更换上游评分追求闭环；本轮未执行该下游部分。

## Methods
细胞输入来自冻结11950恶性细胞目标表达元数据；合格患者需至少50细胞且TYMP与模块每一基因各至少10个检出细胞。核心14患者、SAM16患者；不以效应选人。评分基因为既有production_core和SAM_supply定义，主要评分为患者内log1pCP10K逐基因z标准化后的等权均值，敏感性为逐基因标准化秩均值。协变量复用Phase174并逐细胞ID核对，包含深度、检出基因数、排除原八目标的G2M和Hypoxia。
在患者内进行偏Spearman，患者之间双侧Wilcoxon；每评分/调整两模块BH；5000次患者bootstrap给中位数区间。研究及逐基因删除只描述，不新增检验或替换评分。
公共模块为原Phase162队列内基因log2TPM z均值，已逐队列核对无变更。短路径拟合X~TYMP+C、MTA~TYMP+X+C及MTA~TYMP+C；各变量在队列内标准化秩。32个队列/模块/调整模型各5000次患者bootstrap，共160000次，抽样内重新取秩，模块定义与原评分常数冻结。各队列先算a×b，再以bootstrap方差随机效应合并。条件边及合并乘积按每调整两模块BH，历史FDR另列。
保留raw、MTAP、MTAP+髓系、MTAP+髓系+缺氧四种调整，给随机效应正态近似区间及固定汇总权重的bootstrap区间，后者不重估每次抽样的异质性。未计算中介比例，不要求总效应先显著；SPP1和客户数据不进入模型。

## Discussion与解释范围
细胞内程序共变与公共组织跨患者关系处在不同尺度，不能视作同一患者同一途径的直接复现。RNA评分不代表实际多胺产量或SAM浓度；组织MTA也不代表肿瘤细胞的胞外暴露。背景关联、局部条件关联与中介分解属于不同证据层级。
这是已有数据上的定向追加，不是新独立验证。程序成员预先固定，逐基因删除与两种评分的所有结果保留；不能把最有利子集重新命名为主结果。协变量调整可能刻画共同状态或部分作用过程，因此原始及调整结果并列解释。
SAM主要FDR靠近0.05，更多技术调整后仅名义显著；需要保持小效应、模型敏感的定位。合成核心研究间相反方向及短路径不支持均属于主要结果，不隐藏。

## 验证
1344/1344检查通过。患者相关以独立QR复核，Wilcoxon由独立精确符号秩分布复算；32模型各20次bootstrap独立重算，共640次，最大误差1.67e-15。公共短路径恒等式、条件边精度矩阵、合并统计、FDR、来源哈希和冻结方案均通过。
无新图、无细胞/spot层P值、无客户数据入模。原始结果只读，完整表与失败/可估计性、构成诊断、同患者SMS对照和留一队列结果均保存。全文继续按素材模块拆分整合，免疫功能与空间结构仍为两条独立后续终点。
"""
(O/'POLYAMINE_BACKGROUND_RESULTS_CN.md').write_text(report,encoding='utf-8')
module="""# 可并入统一文稿：TYMP相关SAM及多胺代谢背景

在第一部分已观察到的蛋氨酸/多胺代谢特征基础上，我们以预先定义的ODC1/SRM/SMS合成核心及AMD1/MAT2A供给与脱羧程序作定向补充。16位恶性细胞患者中，SAM相关程序在深度、增殖和缺氧调整后12/16正向，中位rho=0.0181，P=0.0250、FDR=0.0499；进一步调整检出基因数后保留名义正向（rho=0.0155，P=0.0386、FDR=0.0773）。合成核心的未调整关联为12/14正向、rho=0.0708，P/FDR=0.00305，但主调整后降至rho=0.0124、P/FDR=0.626。

公共四队列255患者中，TYMP与合成核心在MTAP及髓系背景调整后保留正向关联（rho=0.2279，P=0.00384、本轮FDR=0.00769）；核心与实测MTA的条件关联较弱（rho=0.1389，P=0.0865、FDR=0.1273），对应短路径合并区间跨零。

这些结果支持保留TYMP相关的广义蛋氨酸/SAM及部分多胺代谢背景，可作为MTA候选的引出依据；当前尚不能据此确认该程序介导TYMP对MTA积累的作用。所有效应按其细胞或组织尺度解释，未将转录评分当作代谢物浓度。
"""
(O/'POLYAMINE_MANUSCRIPT_MODULE_CN.md').write_text(module,encoding='utf-8')
now=datetime.datetime.now().astimezone().isoformat();status=json.loads((O/'COMPUTE_STATUS.json').read_text())
status.update(status='PASS_POLYAMINE_BACKGROUND_ANALYSIS_COMPLETE',validation=V,conclusion='segmented_metabolic_background_supported; core_conditional_cell_association_and_MTA_short_path_not_robust',completed_at=now)
(O/'STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
pd.DataFrame([dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(O.iterdir()) if p.is_file() and p.name not in ['ARTIFACT_HASHES.tsv','FINAL_STATE_VALIDATION.json']]).to_csv(O/'ARTIFACT_HASHES.tsv',sep='\t',index=False)
B=O/'authority_backups';B.mkdir(exist_ok=False);names=['PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','CODEX_HANDOFF.md']
for n in names:shutil.copy2(R/n,B/n)
note=f"""
## Phase175 多胺背景补充完成（{now}）
沿既有ODC1/SRM/SMS合成核心、AMD1/MAT2A SAM供给/脱羧程序，MTAP清除单列。恶性core14、SAM16患者，原评分与秩评分敏感性；无新基因搜索。
SAM主调整rho0.01810、12/16正、P0.02496/FDR0.04993，加检出基因数P0.03864/FDR0.07727名义；core原始rho0.07083/P=FDR0.00305，主调整rho0.01244/P=FDR0.62573，研究差异保留。公共255患者TYMP—core主rho0.22792/P0.003844/FDR0.007689（历史三候选0.005767另存）；core—MTA主rho0.13889/P0.08646/FDR0.12729。
合成核心主短路径a*b0.01471，区间[-0.02964,0.05906]，P0.51557/FDR0.85427；SAM短路径亦不支持。可保留广义代谢背景，不写已确认TYMP经该程序驱动MTA积累。不将各模型较好的一段拼成同一模型整链。
报告：{O/'POLYAMINE_BACKGROUND_RESULTS_CN.md'}
正文模块：{O/'POLYAMINE_MANUSCRIPT_MODULE_CN.md'}
1344检查全PASS，640bootstrap独立复算，无图。背景补充已完成；若继续转向MTA—髓系/SPP1具体证据，不继续更换上游模块追阳性。免疫与空间仍分开，后续未执行。
"""
for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']:
 p=R/n;t=p.with_name(n+'.phase175.tmp');t.write_text(note+p.read_text(),encoding='utf-8');t.replace(p)
p=R/'RESULTS_MANIFEST.tsv'
with p.open(newline='') as f:rows=list(csv.reader(f,delimiter='\t'))
assert len(rows[0])==8 and not any(r and r[0].startswith('phase175_') for r in rows[1:])
reg=['POLYAMINE_BACKGROUND_RESULTS_CN.md','POLYAMINE_MANUSCRIPT_MODULE_CN.md','ALL_STATISTICS_CN.md','PLAN_BEFORE_ANALYSIS_CN.md','MALIGNANT_PROGRAM_SUMMARY.tsv','MALIGNANT_PROGRAM_EFFECTS.tsv','MALIGNANT_ELIGIBILITY.tsv','GENE_LEAVE_ONE_DIAGNOSTIC.tsv','MATCHED_CORE_VS_SMS.tsv','MALIGNANT_STUDY_SENSITIVITY.tsv','PUBLIC_PROGRAM_SHORT_COHORT.tsv','PUBLIC_PROGRAM_SHORT_META.tsv','PUBLIC_PROGRAM_EDGES_META.tsv','LEGACY_PROGRAM_ASSOCIATIONS.tsv','PUBLIC_LEAVE_ONE_COHORT.tsv','PUBLIC_BOOTSTRAPS.npz','INPUT_MANIFEST.tsv','VALIDATION.tsv','VALIDATION_STATUS.json','STATUS.json','run_analysis.py','stats_helpers.py','validate_analysis.py','finalize_analysis.py','ARTIFACT_HASHES.tsv']
for n in reg:rows.append(['phase175_'+n,'Phase175','polyamine_background',str(O/n),now,'PASS','1344 checks; two frozen modules; no new plots','sha256='+sha(O/n)])
t=p.with_name(p.name+'.phase175.tmp')
with t.open('w',newline='') as f:csv.writer(f,delimiter='\t',lineterminator='\n').writerows(rows)
t.replace(p)
gate=dict(stage='Phase175',registered=len(reg),paths_exist=all((O/n).exists() for n in reg),state_updated=all((R/n).read_text().startswith(note) for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']))
assert gate['paths_exist'] and gate['state_updated'];(O/'FINAL_STATE_VALIDATION.json').write_text(json.dumps(gate,indent=2))
D=O/'delivery';D.mkdir()
for n in ['POLYAMINE_BACKGROUND_RESULTS_CN.md','POLYAMINE_MANUSCRIPT_MODULE_CN.md','ALL_STATISTICS_CN.md','MALIGNANT_PROGRAM_SUMMARY.tsv','PUBLIC_PROGRAM_EDGES_META.tsv','PUBLIC_PROGRAM_SHORT_META.tsv']:shutil.copy2(O/n,D/n)
pd.DataFrame([dict(name=p.name,sha256=sha(p)) for p in sorted(D.iterdir())]).to_csv(O/'DELIVERY_HASHES.tsv',sep='\t',index=False)
print(json.dumps(status,ensure_ascii=False));print(json.dumps(gate))
