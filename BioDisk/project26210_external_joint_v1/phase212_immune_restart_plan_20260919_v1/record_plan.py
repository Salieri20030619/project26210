from pathlib import Path
import json,csv,hashlib,datetime,shutil
H=Path.home();R=H/'project26210_external_joint_v1';O=R/'phase212_immune_restart_plan_20260919_v1';F=R/'phase211_immune_receiver_refreeze_20260919_v1/freeze_v2'
assert json.loads((F/'STATUS.json').read_text())['status']=='FROZEN_VALIDATED'
plan='''# 26210 免疫分析重启计划

Phase212；2026-09-19；状态：PLAN_ONLY_NOT_STARTED。用户要求“现在重启免疫分析，先列计划”。本轮只提交计划和预设假设，不启动评分、关联、通信或新图。第一、二部分及冻结C4不变。分析入口唯一为Phase211/freeze_v2。

## 1. 主问题与叙事

承接已完成的TYMP相关代谢发现、MTA候选及SPP1富集巨噬定位，回答：SPP1相关巨噬状态是否对应CD8功能变化、Treg组成变化及NK功能变化？

沿既定双层分工：整体成熟巨噬SPP1是主要预测量；C4比例和C4内SPP1用于说明该联系是否与已定位的C4状态相容。不能以SPP1或C4与某端点相关就确认免疫抑制机制，也不以空间结果代替本分支。

不将MTA转录代理重新作为免疫端点的必经预测量，不重做前两部分，不为形成完整串联中介追加模型。主要反向或未支持结果完整保留。

## 2. 预设端点、样本和分析层级

| 层级 | 预测量 | 结局 | 冻结患者/文库 | 预期方向 |
|---|---|---|---|---|
| 主要功能A1 | 整体成熟巨噬SPP1 | CD8细胞毒 | 16 / 19 | 负 |
| 主要功能A2 | 整体成熟巨噬SPP1 | CD8功能障碍相关程序 | 16 / 19 | 正 |
| 主要组成A3 | 整体成熟巨噬SPP1 | Treg / 已分类CD4 | 13 / 16 | 正 |
| 次要功能A4 | 整体成熟巨噬SPP1 | NK细胞毒 | 13 / 16 | 负 |
| 次要状态定位B | 冻结C4比例 | 同上四端点 | 各16或13 | 同上 |
| 次要群内表达C | C4内SPP1 | 同上四端点 | CD8各11；NK/Treg各10 | 同上 |

预期方向用于解释，不用于单侧检验、筛样或隐藏反向结果。C4比例、C4内表达与整体巨噬SPP1不是互相可替换的指标；不因其中一个不支持就更换主要预测量。C4内表达的患者名单和文库集合单独继承冻结表。

CD8研究构成为GSE159115 2例、GSE178481 4例、GSE207493 9例、GSE242299 1例；NK/Treg为GSE178481 4例、GSE207493 9例。不是4个均有充分样本的独立验证研究。

## 3. 第一阶段：生成可核对的患者层输入

1. 核对Phase211冻结清单SHA256及必需对象可读性；仅做必要完整性检查，不重跑注释、聚类、全项目审计。
2. 在固定合格肿瘤文库内分别汇总：成熟巨噬SPP1原始counts伪汇总log2(CPM+1)、C4/成熟巨噬比例、C4内SPP1伪汇总；不因SPP1为零删除样本。
3. 在固定接收谱系内按基因集等权计算每细胞log1p(CP10K)，再细胞等权得到文库程序均值。CD8与NK分开；不混合T+NK评分。
4. Treg分母固定为CD4_Tcon+CD4_Treg；保留符合冻结资格的零Treg。不得按每细胞FOXP3重新筛选。
5. 对每个端点，以该端点同一合格文库集合等权汇总预测量与结局，得到每患者一个值。C4内表达用其独立合格文库集合，不能拼接其他端点文库。比例同样先文库比例后患者等权，不让文库细胞数变相决定患者权重。
6. 不做全队列重新标准化或按期望方向反向筛基因；记录完整基因覆盖、接收端细胞数、计数深度、检出基因数及患者/文库名单。

功能集合继承freeze_v2：细胞毒NKG7/GNLY/GZMB/PRF1/CTSW；CD8功能障碍PDCD1/CTLA4/LAG3/HAVCR2/TIGIT/TOX/CXCL13/ENTPD1。另计算Phase211既定活化面板FOS/JUN/JUNB/NR4A1/CD69/IL2RA用于解释，不将其与功能障碍合成一个免疫抑制分数，不据它筛细胞。

本阶段验收：16/13及C4群内11/10的名单与冻结表一致；患者不重复计数；源counts抽取与聚合可复核；SPP1零值、Treg分母、固定评分成员正确。若关键输入确实缺失或常量，不补值、不换患者；登记为不可估计，保留原冻结身份。

交付：文库指标表、患者指标表、评分覆盖、数据字典和技术资格核对表。此时不画图。

## 4. 第二阶段：主分析与C4对应分析

### 4.1 主要统计

主要效应量为研究来源调整的偏Spearman相关：对患者层X/Y分别作平均秩，分别回归研究分类变量及截距，计算两组秩残差相关。双侧P来自研究内置换Y，保留各研究患者数及标签；随机种子262212。若允许的独特置换数不超过99999，枚举；否则固定99999次Monte Carlo置换，采用加一修正。对不可估计的残差变异输出NA及原因，不输出伪P=1。

报告独立患者总n及提供研究内变化信息的研究/患者数。CD8中仅1例的研究保留在固定样本账本，但不提供研究内变化信息；不能把它写作第四个内部重复验证。C4群内子集亦逐项登记这一问题。

效应区间采用患者层、研究内分层bootstrap 5000次，仅作小样本不确定性描述，不用于覆盖置换检验；有效重复不足80%则不报告数值区间，说明不可稳定估计。未经研究调整的Spearman rho及双侧P作为对照，不在两套结果中挑更显著者充当主结果。

分研究展示方向与散布；各研究n<10，作为探索性/描述性结果，不要求它们各自显著，也不对n=1计算相关。合并方向不等于每研究/每患者一致，不能把细胞数当统计n。

### 4.2 多重检验

完整12项预设关联登记于HYPOTHESIS_REGISTER.tsv，按同类问题分家族：
- 主要SPP1→CD8两个功能端点：BH，2项。
- 主要SPP1→Treg/CD4：单个预设组成端点，不另作FDR校正。
- 次要SPP1→NK：单个预设次要功能端点，不另作FDR校正，不升级为主要终点。
- 两种C4预测量→CD8：4项次要功能家族，共同BH。
- 两种C4预测量→Treg/CD4：2项次要组成家族，BH。
- 两种C4预测量→NK：2项次要功能家族，BH。

原假设家族保留不可估计条目身份；正式表P/FDR显示NA。计算家族BH时不可估计项作为p=1占据原家族大小，只参与校正，不将其解释为有效阴性检验。未经调整对照与各技术检查不混入主要家族，不能用它们替代主要P/FDR或缩小原家族。历史Phase72的P/FDR保持原定义。

### 4.3 结果解释

依次报告效应方向/大小、研究来源调整结果、患者层不确定性及P/FDR。有FDR的家族以FDR判定校正后支持；P<0.05而FDR未过只写名义显著。若正向CD8细胞毒同时伴随功能障碍升高，报告共存状态，不能概括成杀伤降低。组成端点只能支持Treg比例联系，不能直接证明抑制功能。

C4与整体巨噬SPP1不同样本量的相关强弱不得直接排大小；需要并排对照时，另列共同患者集合的描述性效应，保留原主要分析，不选最有利集合。

## 5. 第三阶段：有限稳健性及功能解释

所有预设检查不以主分析显著为启动门槛，也不因阴性无限增加模型。

1. 留一患者和留一研究：主要查看效应方向、范围、是否被单例/单研究支配，不为每次删除生成一套“显著性结论”。删除GSE207493后样本很小，明确描述性，不能以其不显著否定原效应。
2. 单个技术协变量：在研究调整基础上，分别加入接收端患者平均检出基因数，检查功能分数与测序复杂度关系；每次只加一个预设变量，残差自由度<5时不拟合。Treg组成不套用功能程序的深度模型。
3. 固定标签复核：利用Phase211已有患者集中/身份置信度标记，单独评估排除已标记群后的影响。只使用源冻结细胞列表派生子集，原名单和主结果不变；样本减少/不可估计完整登记。不得把全部低表达细胞删除后声称验证了全体接收端。
4. 功能障碍与活化分开：并列患者层两个程序的实际变化；活化调整最多作为一个预设次要解释模型，模型资格不足时停在描述。不得因调整后P更小就改主要定义。

本轮不新增一批基因集、切分点、反复门槛或高低组寻找阳性；冻结时10/20/30%未知身份对照已归档，主规则继续20%。

交付：主结果、完整统计、患者影响及研究异质性表；正文核心保留主要阴性和反向结果，次要技术检查放补充。

## 6. 第四阶段：定向通信与接收端对应

主统计完成后，以既有CellChat/NicheNet及配体受体资源为起点，先核对版本、物种、复合受体组成和目标细胞可测性，再固定SPP1相关候选配体—受体列表。资源与候选表冻结后才计算关系，不按本次免疫P值挑通路。

区分整体成熟巨噬、C4及其他成熟巨噬的配体来源；C4相关比较继承有足够C4细胞的固定子集。对CD8/CD4/NK分别检查受体及必需复合亚基，先判可估计性；若有合格的已知靶基因集合，再作接收端响应对应。避免把用于定义免疫结局的同一套基因重复当成独立机制验证。

这一阶段先交付候选来源/受体/患者覆盖表，再做有限的患者级定向统计；不自动重跑全量CellChat，不把软件通信分数当实验配体结合，不以某条LR阳性替代主要免疫端点阴性。缺受体或配对资格时该条停止，不强行补链。新增LR假设在独立同类家族内校正。

## 7. 外部对应的范围

GSE294109严格按未治疗1例、Axitinib/Pembrolizumab 2例、Ipilimumab/Nivolumab 4例、Lenvatinib/Pembrolizumab 2例分开。只生成固定谱系的患者指标与方向性对照，不作正式重复验证、不对这些层生成推断性P/FDR；尤其不对1例计算相关。分层过小导致某个方向本身不可评估时明确记为不可评估。

GSE264586继续保持未定群和肿瘤取材资格缺口的登记，不在本轮静默加入；若以后补齐原始来源，另建冻结版本。外部旧NK混合群不能当常规NK，外部Macrophage_SPP1不能当冻结C4。GSE207493与发现队列重叠，不作新增独立外部验证。

现有Phase72是旧定义结果，作为历史对照，关键阴性保留；不得把本次重新注释后的不同端点描述成已推翻或已独立验证旧结论。本轮不追加bulk临床疗效、生存或空间模型。

## 8. 执行顺序、交付与停止点

| 顺序 | 工作 | 验收与交付 |
|---|---|---|
| I1 | 冻结输入→文库→患者指标 | 名单、分母、评分覆盖和聚合核验；无关联P |
| I2 | 主要三端点及NK次要端点；C4定位对应 | 主效应/P/FDR、完整12项结果、不可估计原因 |
| I3 | 患者/研究影响及有限解释检查 | 主结论是否被个别样本、研究或深度驱动 |
| I4 | SPP1定向通信资格与有限对应 | 来源—受体—接收端证据表，保留失败条目 |
| I5 | 外部小样本描述性对应与文字整合 | 各治疗层分开，明确不属于正式重复验证 |

先交文字与统计表，审核内容后再安排图件。各阶段使用版本化服务器目录，保存输入身份、脚本、软件、日志及必要独立数值核验；完成后更新四份项目权威记录。无需中间判断的长作业使用原生等待或一个持久监控程序，结束以退出状态、完成标记和必要输出完整性共同判定。

停止/降级规则：身份或测量资格不足、常量或自由度不足时登记不可估计；有效分析不支持时如实报告，不通过改主要端点、聚类或分母追求阳性。免疫分支是否成立由自身数据决定，不以空间/临床/客户湿实验填补缺失结果。

**本次仅计划，I1—I5均未启动。建议下一执行单元为I1，随后按计划进入I2。**
'''
(O/'IMMUNE_RESTART_PLAN_CN.md').write_text(plan)
endpoints=[('CD8_cytotoxicity','negative',16,11),('CD8_dysfunction','positive',16,11),('Treg_CD4','positive',13,10),('NK_cytotoxicity','negative',13,10)]
rows=[]
for pred in ['mature_SPP1_pseudobulk','C4_fraction','C4_SPP1_pseudobulk']:
 for ep,direction,n,nc in endpoints:
  typ='CD8_function' if ep.startswith('CD8') else ('Treg_composition' if ep=='Treg_CD4' else 'NK_function')
  family=('SPP1_' if pred=='mature_SPP1_pseudobulk' else 'secondary_C4_')+typ
  size=(2 if typ=='CD8_function' else 1)*(1 if pred=='mature_SPP1_pseudobulk' else 2)
  rows.append(dict(hypothesis_id=f'I{len(rows)+1:02d}',predictor=pred,endpoint=ep,expected_direction=direction,analysis_tier='primary' if pred=='mature_SPP1_pseudobulk' and ep!='NK_cytotoxicity' else 'secondary',frozen_patient_n=nc if pred=='C4_SPP1_pseudobulk' else n,family=family,family_size=size,adjustment='BH' if size>1 else 'none_prespecified_single_endpoint',test='two_sided_study_blocked_permutation_partial_Spearman',execution_status='NOT_STARTED'))
with (O/'HYPOTHESIS_REGISTER.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
assert len(rows)==12
for family in set(x['family'] for x in rows):
 subset=[x for x in rows if x['family']==family];assert len(subset)==subset[0]['family_size']
source=[F/'FREEZE_SHA256.json',F/'FREEZE_DEFINITION.json',F/'DISCOVERY_PATIENT_ENDPOINT_FREEZE.tsv',F/'DISCOVERY_ENDPOINT_STUDY_COMPOSITION.tsv',F/'EXTERNAL_PATIENT_ENDPOINT_FREEZE.tsv']
(O/'INPUT_PLAN_SOURCE_SHA256.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in source},indent=2))
now=datetime.datetime.now().astimezone().isoformat();(O/'STATUS.json').write_text(json.dumps({'status':'PLAN_ONLY_NOT_STARTED','created_at':now,'hypotheses':12,'family_sizes_consistent':True,'source_freeze':str(F),'analysis_jobs_launched':False,'figures_created':False},indent=2))
backup=O/'root_record_backups';backup.mkdir(exist_ok=False)
for n in ['PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','CODEX_HANDOFF.md']:shutil.copy2(R/n,backup/n)
note='''\n## 2026-09-19 Phase212：免疫分析重启计划，尚未执行\n用户要求“现在重启免疫分析，先列计划”。已保存phase212_immune_restart_plan_20260919_v1/IMMUNE_RESTART_PLAN_CN.md及12项HYPOTHESIS_REGISTER.tsv；源Phase211/freeze_v2不变。依次I1患者指标、I2主统计/C4对应、I3有限稳健性、I4定向通信、I5外部描述性对应与文字整合。主要SPP1→CD8两个功能端点及Treg组成，NK次要；C4比例/群内SPP1为次要定位。CD8 n16（研究2/4/9/1），NK/Treg n13（4/9）；主效应研究调整偏Spearman+研究内双侧置换，原始相关作对照，FDR按预设同类家族2/1/1/4/2/2。外部1/2/4/2例分治疗描述，不称正式验证；不启动空间/临床、重注释或绘图。计划已交付，所有I1-I5尚未启动；下一执行入口I1。\n'''
for n in ['PROJECT_STATUS.md','DECISIONS.md']:
 with (R/n).open('a') as f:f.write(note)
p=R/'CODEX_HANDOFF.md';p.write_text(note+'\n'+p.read_text())
with (R/'RESULTS_MANIFEST.tsv').open('a') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n')
 for p in sorted(O.iterdir()):
  if p.is_file():w.writerow(['phase212_'+p.stem,'phase212','immune_restart_plan',str(p),now,'PLAN_ONLY_NOT_STARTED',str(O/'STATUS.json'),'12 hypotheses; no analysis launched'])
print('PLAN_COMPLETE',len(rows),'hypotheses; no analysis launched')
