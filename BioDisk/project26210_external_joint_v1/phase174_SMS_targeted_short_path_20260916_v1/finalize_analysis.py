from pathlib import Path
import json,hashlib,csv,shutil,datetime
import pandas as pd
O=Path(__file__).resolve().parent;R=O.parent
V=json.loads((O/'VALIDATION_STATUS.json').read_text());assert V['checks']==V['passed']
S=pd.read_csv(O/'MALIGNANT_SMS_SUMMARY.tsv',sep='\t');F=pd.read_csv(O/'PUBLIC_SHORT_PATH_COHORT.tsv',sep='\t');M=pd.read_csv(O/'PUBLIC_SHORT_PATH_META.tsv',sep='\t');E=pd.read_csv(O/'PUBLIC_EDGE_EVIDENCE_HISTORY.tsv',sep='\t')
def f(x):return '-' if pd.isna(x) else f'{x:.6g}'
st=['| 模型 | 患者n | 正/负 | 中位rho | 中位数95%区间 | P | 原六基因FDR | 新四模型FDR |','|---|---:|---:|---:|---|---:|---:|---:|']
for x in S.itertuples(index=False):st.append(f'| {x.model} | {x.n_patients} | {x.positive}/{x.negative} | {f(x.median_rho)} | [{f(x.median_low)}, {f(x.median_high)}] | {f(x.P)} | {f(x.historical_FDR_axis6)} | {f(x.supplement_FDR_4)} |')
ct=['| 队列 | n | TYMP—SMS偏rho / P | SMS—MTA偏rho / P | a×b | 患者bootstrap 95%区间 |','|---|---:|---|---|---:|---|']
for x in F[F.adjustment=='MTAP_myeloid'].itertuples(index=False):ct.append(f'| {x.cohort} | {x.n} | {f(x.a_rho)} / {f(x.a_P)} | {f(x.b_rho)} / {f(x.b_P)} | {f(x.ab)} | [{f(x.ab_low)}, {f(x.ab_high)}] |')
mt=['| 调整 | 合并a×b | 随机效应95%区间 | P | 三敏感性FDR | 固定权重患者bootstrap区间 |','|---|---:|---|---:|---:|---|']
for x in M.itertuples(index=False):mt.append(f'| {x.adjustment} | {f(x.estimate)} | [{f(x.low)}, {f(x.high)}] | {f(x.P)} | {f(x.sensitivity_FDR3)} | [{f(x.bootstrap_low)}, {f(x.bootstrap_high)}] |')
report="""# SMS专门补充：SMS—MTA联系保留，TYMP—SMS及整体短路径尚不稳定

## 总体判断
本轮把TYMP—SMS—MTA从先前包含SPP1的长串联中单独拆出，并直接补查第一部分恶性细胞内的TYMP—SMS关系。结果显示，上游在技术与状态调整后接近零；四公共队列中SMS—MTA仍有支持，但TYMP—SMS存在队列异质性，短路径合并区间跨零。
因此，SMS目前可保留为“有局部公共队列支持的候选合成环节”，不宜写为已确定的TYMP至MTA必经中介。这一判断不单纯依据FDR门槛，而同时考虑调整后效应接近零、方向一致性和跨队列差异。不会因本轮结果删除SMS—MTA事实、客户实验或整个MTA候选。
本轮仅SMS，无SRM/多胺模块新检验，无SPP1进入短路径，无新图。

## 结果一：恶性细胞中TYMP—SMS的直接关系调整后接近零
使用同一冻结16名合格患者、每患者至少50恶性细胞且TYMP/SMS各至少10个检出细胞，不只选择双阳性细胞。原48个患者×模型效应（raw/depth/depth_cycle）得到精确复现。
未调整13/16患者正向，中位rho=0.03113，P=0.09344、原六基因FDR=0.09344，保留弱正向线索。仅深度调整后8/16正，中位rho=-0.002037，P=0.63217、原FDR=0.75861。主要追加模型（深度+G2M+缺氧）9/16正，中位rho=0.001869，P=0.74356、四追加模型FDR=0.74356；进一步加检出基因数后8/16正，中位rho=-0.000955。
四个新增模型的效应均接近零。既有TYMP阳性/阴性SMS表达比较及SMS—F1、轨迹结果均不被本轮替换；SMS—F1不能当作TYMP—SMS直接证据。

## 恶性细胞完整模型表
""" + '\n'.join(st)+"""

主要模型的分研究描述为：GSE156632 2/4正、中位0.004243；GSE159115 1例正0.10057；GSE178481 1例负-0.18546；GSE207493 6/10正、中位0.001869。去GSE207493剩6例，3/6正、中位0.004243，只作探索性描述。
原始模型的bootstrap中位数区间为正，而原Wilcoxon P未过0.05，两者保留；前者针对中位数抽样分布，后者利用带符号秩，不能将该原始区间替代调整后主模型或改换检验寻求显著。

## 结果二：公共组织SMS—MTA保留，上游TYMP—SMS存在明显异质性
四个配对肿瘤队列共255患者；主要调整为MTAP和髓系RNA代理。合并TYMP—SMS偏rho=0.21730，P=0.06146、历史三候选FDR=0.06146，3/4队列正向，I²=68.08%。未达到名义显著，保留探索性方向。
控制TYMP及上述背景后，SMS—实测MTA四队列均正，rho=0.24993、P=0.000578。原Phase162六基因家族FDR=0.002888，Phase163三候选FDR=0.001733；本轮追加二边FDR=0.001155。三个口径完整保存，正文优先沿用原六基因FDR，不将更小校正池描述成新独立验证。

## 结果三：短路径在CPTAC有探索性支持，合并尚不稳健
主要模型先在各队列拟合SMS~TYMP+C和MTA~TYMP+SMS+C，得到各队列标准化秩系数a、b与乘积a×b；未将两个合并相关系数相乘。c=cprime+a×b恒等式核验通过。
""" + '\n'.join(ct)+"""

CPTAC的正向乘积及区间不跨零是本轮需要保留的单队列线索，不等于四队列共同确认中介。其他三个队列的乘积bootstrap区间均跨零，RC18点估计为负。逐患者删除时CPTAC乘积范围0.16569—0.21921仍为正，RC18范围-0.03632至-0.002115仍为负，说明该差异不能简单归因于单个患者。

四队列随机效应合并a×b=0.04484，95%区间[-0.02692,0.11659]，P=0.22069，I²=39.02%。固定点估计汇总权重的患者bootstrap区间为[-0.01453,0.09395]，同样跨零。该bootstrap正向比例0.9268只是抽样方向比例，不是P值，也不是“92.68%概率存在因果中介”。
留一队列后全部正态近似区间仍跨零；去CPTAC后合并乘积降至0.01551，P=0.56976。其他调整亦未建立稳定的合并短路径。

## 短路径合并与敏感性
""" + '\n'.join(mt)+"""

## 对正文线路的处理
1. 保留第一部分的核苷/蛋氨酸/多胺代谢结果。
2. 保留SMS与实测MTA的分段联系，以及CPTAC短路径探索性线索。
3. 不写“TYMP通过上调SMS导致MTA积累”，也不将CPTAC单队列提升为普遍机制。
4. 可写：“SMS与实测MTA存在稳定分段联系，但其与TYMP的上游关联在恶性细胞及跨队列调整中尚不稳定，暂不足以指定SMS为TYMP—MTA之间的必经环节。”
5. 按用户顺序，本轮SMS定向核查结束；如继续，下一步讨论更广的多胺代谢背景，而非继续为SMS追加模型。本轮未执行多胺整体补充，也未改变第一部分或客户实验结论。
6. 全文仍按统一文稿拆分整合，免疫功能和空间结构两个终点保持独立。

## Methods
恶性细胞沿用Phase162的11950细胞和SMS合格16患者；读取原始count矩阵，核对TYMP/SMS counts和总counts与冻结元数据一致。G2M和Hypoxia基因集来自已用Hallmark2020资源，排除原八个目标基因；患者内对有变异基因的log1p CP10K标准化后求均值。偏Spearman以患者内平均秩、协变量残差相关计算；患者间采用双侧Wilcoxon符号秩。四新增模型共同BH，旧模型保留六基因家族。患者中位数区间为5000次患者bootstrap。
公共数据使用原四队列实测MTA配对值，不含预测MTA或客户数据，未改变患者清单。A=TYMP、X=SMS、Y=MTA，各队列内标准化秩；拟合X~A+C、Y~A+X+C、Y~A+C，保存a、b、a×b、direct、total。每个队列每种调整5000次患者重抽样，重新取秩和标准化，共16模型80000抽样。raw、MTAP、MTAP+髓系、MTAP+髓系+缺氧全保留。
各条件边用Fisher-z随机效应合并；乘积在队列内先算，再以bootstrap方差按DL随机效应汇总，给正态近似区间/P。另用该模型固定汇总权重组合独立队列bootstrap抽样分布；后者不重新估计每次抽样的异质性。主要合并乘积为单个预定追加端点，未做多候选筛选；三个敏感性乘积P单独BH。各队列/模型及合并层的两条边作为追加家族，同时保留历史三/六候选FDR。
主要模型进行逐患者删除，各调整进行留一队列；小样本研究子集仅描述。未计算中介百分比，不要求总效应先显著。

## Discussion与解释范围
本轮为已有数据上的SMS定向追加，既有部分边已被查看，不是新独立队列或独立预注册验证。细胞内共表达、患者间组织相关和短路径系数回答不同层面的问题，不能相互替代。
横断面方向由拟定模型指定，a×b只表示关联分解，不证明时间先后、因果中介或MTA净生成。仅四队列，DL方差和正态近似区间的稳定性有限；保留患者bootstrap及留一队列结果。细胞RNA不能代表酶活，组织MTA不能替代胞外浓度。协变量可能刻画共同状态或部分作用路径，因此调整与未调整并列，不将调整后近零等同因果排除。
细胞资格使正式直接检验为16患者，不能把总26患者作为该端点n。SMS—F1并非TYMP—SMS，正向表达关联也不直接证明酶活上升。
第一部分已建立的代谢状态、既有SMS—MTA分段事实、客户shTYMP代谢变化均继续保留。当前不支持的是把SMS指定为已稳定贯通的必经中介。

## 复核及档案
786/786检查通过；恶性细胞偏相关经独立QR复算，患者Wilcoxon用独立精确符号秩分布复核。16个公共模型各独立复算25个bootstrap样本，共400次，最大误差1.89e-15；条件边用精度矩阵复核，短路径恒等式、合并乘积、区间、FDR及输入哈希均通过。
主要输出：MALIGNANT_SMS_SUMMARY.tsv、PUBLIC_SHORT_PATH_COHORT.tsv、PUBLIC_SHORT_PATH_META.tsv、PUBLIC_EDGE_EVIDENCE_HISTORY.tsv、两类留一结果及完整bootstrap/输入缓存。具体方案、样本、资源基因集、脚本、环境与来源哈希均在本目录。没有新图或SPP1模型，没有启动其他多胺候选。
"""
(O/'SMS_TARGETED_RESULTS_CN.md').write_text(report,encoding='utf-8')
short="""# 可并入全文的SMS补充结果

在第一部分恶性上皮细胞的16名合格患者中，TYMP与SMS的未调整相关呈弱正向趋势（13/16正，中位rho=0.0311，P=0.0934、原六基因FDR=0.0934），但深度、增殖及缺氧调整后中位rho接近零（0.00187，9/16正，P=0.7436、追加四模型FDR=0.7436）。当前数据不支持稳定的恶性细胞TYMP—SMS直接正关联。

公共四队列255患者中，控制TYMP、MTAP及髓系背景后，SMS与实测MTA仍四队列同向（合并rho=0.2499，P=0.000578、原六基因FDR=0.00289）。TYMP—SMS上游则为3/4队列正向，rho=0.2173，P=0.0615、历史三候选FDR=0.0615，存在队列异质性。独立拆出的TYMP—SMS—MTA短路径在CPTAC获得单队列正向探索性支持，但四队列合并乘积0.0448的95%区间为[-0.0269,0.1166]，P=0.2207；患者bootstrap敏感性区间亦跨零。

因此，SMS可以保留为与MTA合成背景相关的候选环节，但尚不足以指定为TYMP至MTA的稳定必经中介。该结果与第一部分已建立的更广代谢状态并列解释，不用单队列线索覆盖上游弱关联，也不否定客户已有TYMP干预后MTA变化。
"""
(O/'SMS_MANUSCRIPT_MODULE_CN.md').write_text(short,encoding='utf-8')
now=datetime.datetime.now().astimezone().isoformat()
status=json.loads((O/'COMPUTE_STATUS.json').read_text());status.update(status='PASS_SMS_TARGETED_ANALYSIS_COMPLETE',conclusion='SMS_MTA_supported_but_TYMP_SMS_and_pooled_short_path_not_robust',single_cohort_signal='CPTAC positive bootstrap interval; exploratory',polyamine_expansion_started=False,validation=V,completed_at=now)
(O/'STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
pd.DataFrame([dict(path=str(p),size=p.stat().st_size,sha256=sha(p)) for p in sorted(O.iterdir()) if p.is_file() and p.name not in ['ARTIFACT_HASHES.tsv','FINAL_STATE_VALIDATION.json']]).to_csv(O/'ARTIFACT_HASHES.tsv',sep='\t',index=False)
B=O/'authority_backups';B.mkdir(exist_ok=False);names=['PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','CODEX_HANDOFF.md']
for n in names:shutil.copy2(R/n,B/n)
notice=f"""
## Phase174 SMS定向补充完成（{now}）
用户要求先SMS、再考虑多胺。冻结同16恶性患者，复现raw/depth/depth_cycle并追加四模型；深度+周期+缺氧rho中位0.001869、9/16正，P/FDR0.74356。未经调整13/16正rho0.03113/P=原FDR0.09344保留。SMS—F1不代替TYMP—SMS。
短路径仅TYMP—SMS—实测MTA，不含SPP1。四队列255患者、16模型各5000bootstrap。主要MTAP+髓系：a偏rho0.2173/P0.06146、b偏rho0.24993/P0.000578（原六基因FDR0.002888）。CPTAC a*b0.1946、bootstrap区间[0.0171,0.3717]单队列探索性支持；RC18反向。合并a*b0.04484，区间[-0.02692,0.11659]，P0.22069；另一bootstrap区间亦跨零。当前不指定SMS为稳定必经中介，已有SMS—MTA事实保留。
报告：{O/'SMS_TARGETED_RESULTS_CN.md'}
正文模块：{O/'SMS_MANUSCRIPT_MODULE_CN.md'}
786检查全PASS，400bootstrap独立复算；无新图、无其他多胺候选分析。SMS定向检验结束，若继续可讨论多胺整体路线，不继续为SMS追加模型追求阳性。
"""
for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']:
 p=R/n;t=p.with_name(n+'.phase174.tmp');t.write_text(notice+p.read_text(),encoding='utf-8');t.replace(p)
p=R/'RESULTS_MANIFEST.tsv'
with p.open(newline='') as f:rr=list(csv.reader(f,delimiter='\t'))
assert len(rr[0])==8 and not any(x and x[0].startswith('phase174_') for x in rr[1:])
reg=['SMS_TARGETED_RESULTS_CN.md','SMS_MANUSCRIPT_MODULE_CN.md','PLAN_BEFORE_ANALYSIS_CN.md','MALIGNANT_SMS_SUMMARY.tsv','MALIGNANT_PATIENT_SMS_EFFECTS.tsv','MALIGNANT_STUDY_SENSITIVITY.tsv','MALIGNANT_LEAVE_ONE_PATIENT.tsv','LEGACY_SMS_REPRODUCTION.tsv','PUBLIC_SHORT_PATH_COHORT.tsv','PUBLIC_SHORT_PATH_META.tsv','PUBLIC_EDGE_EVIDENCE_HISTORY.tsv','PUBLIC_LEAVE_ONE_PATIENT.tsv','PUBLIC_LEAVE_ONE_COHORT.tsv','PUBLIC_AB_BOOTSTRAPS.npz','PUBLIC_INPUT_CACHE.npz','COVARIATE_GENESETS.json','INPUT_MANIFEST.tsv','VALIDATION.tsv','VALIDATION_STATUS.json','STATUS.json','run_analysis.py','validate_analysis.py','finalize_analysis.py','ARTIFACT_HASHES.tsv']
for n in reg:rr.append(['phase174_'+n,'Phase174','SMS_targeted_short_path',str(O/n),now,'PASS','786 checks; SMS targeted; no other candidate expansion','sha256='+sha(O/n)])
t=p.with_name(p.name+'.phase174.tmp')
with t.open('w',newline='') as f:csv.writer(f,delimiter='\t',lineterminator='\n').writerows(rr)
t.replace(p)
gate=dict(stage='Phase174',backups=4,registered=len(reg),paths_exist=all((O/n).exists() for n in reg),state_updated=all((R/n).read_text().startswith(notice) for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']))
assert gate['paths_exist'] and gate['state_updated']
(O/'FINAL_STATE_VALIDATION.json').write_text(json.dumps(gate,indent=2))
D=O/'delivery';D.mkdir()
for n in ['SMS_TARGETED_RESULTS_CN.md','SMS_MANUSCRIPT_MODULE_CN.md','MALIGNANT_SMS_SUMMARY.tsv','PUBLIC_SHORT_PATH_COHORT.tsv','PUBLIC_SHORT_PATH_META.tsv','PUBLIC_EDGE_EVIDENCE_HISTORY.tsv']:shutil.copy2(O/n,D/n)
pd.DataFrame([dict(name=p.name,sha256=sha(p)) for p in sorted(D.iterdir())]).to_csv(O/'DELIVERY_HASHES.tsv',sep='\t',index=False)
print(json.dumps(status,ensure_ascii=False));print(json.dumps(gate))
