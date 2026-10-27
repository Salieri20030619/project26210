from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,shutil,datetime
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=P/'phase161_SLC16A3_cell_origin_20260915_v1'
R=pd.read_csv(O/'PATIENT_WITHIN_CELL_CORRELATIONS.tsv',sep='\t');S=pd.read_csv(O/'PATIENT_LEVEL_CORRELATION_SUMMARY.tsv',sep='\t');SRC=pd.read_csv(O/'COMPARTMENT_SOURCE_SUMMARY.tsv',sep='\t');D=pd.read_csv(O/'C4_OTHER_COMMON_PATIENT_CORRELATION_DIFFERENCE.tsv',sep='\t');X=pd.read_csv(O/'CROSS_COMPARTMENT_COORDINATION.tsv',sep='\t');QC=pd.read_csv(O/'PATIENT_COVERAGE_QC.tsv',sep='\t')
def md(df):
 out=['| '+' | '.join(df.columns)+' |','|'+'|'.join(['---']*len(df.columns))+'|']
 for row in df.itertuples(index=False,name=None):out.append('| '+' | '.join(f'{x:.6g}' if isinstance(x,(float,np.floating)) else str(x) for x in row)+' |')
 return '\n'.join(out)
summary="""# SLC16A3细胞来源核查及“实测乳酸”解释

实测乳酸是公共肾癌研究对患者肿瘤组织进行质谱代谢组检测得到的乳酸相对丰度。上一轮使用RC18/RC20/CPTAC/CPTAC_val的97/52/50/56名肿瘤患者，共255名，与同患者RNA配对。它来自已发布的参考测量矩阵，并非根据RNA预测；是半定量丰度，不是统一的绝对浓度，也不代表恶性细胞向巨噬细胞释放的净通量。来源：https://www.nature.com/articles/s43018-025-00943-0

本轮直接回查冻结原始单细胞counts，保留同16患者共23018细胞，统计仍以患者为单位。

1. SLC16A3在恶性细胞和巨噬细胞均有表达。16/16患者恶性细胞伪bulk表达高于成熟巨噬，中位log2CPM差+1.205，P=3.05e-5、FDR=6.10e-5。两区室表达不等于已发生转运。
2. 巨噬细胞内SLC16A3—SPP1在控制深度与缺氧后，16患者中14正，患者内rho中位0.103，患者级P=0.00269、FDR=0.00269。关联幅度小，但方向一致。
3. C4中符合覆盖标准的9患者全部为正，rho中位0.243；小样本探索性。与相同9患者其他巨噬比较，8/9患者C4相关更高，配对差中位0.131、P=0.0273、FDR=0.0547，为名义/探索性支持。
4. 恶性细胞内TYMP—SLC16A3上游仍弱。深度调整后15患者12正，rho中位0.065、P=0.0833、FDR=0.1665；再加缺氧后11正，rho中位0.047、P=0.1354、FDR=0.2708。保留方向，不宣称已确认TYMP直接调控MCT4。
5. SPP1也在恶性细胞表达，13/16患者其恶性伪bulk表达高于巨噬。外部bulk的SPP1不能全部解释成巨噬来源；之前明确使用巨噬区室SPP1的单细胞结果继续独立保留。

患者内细胞共表达、跨患者亚群平均表达、不同区室配对相关是三个不同层次。C4内新的正向共表达不推翻此前跨患者平均表达的弱关联。维持整体巨噬SPP1功能联系+C4胆固醇处理状态定位，MTA与客户实验背景保留，不预设乳酸为唯一中介。615项计算核查通过，无新图。
"""
(O/'SUMMARY_CN.md').write_text(summary,encoding='utf-8')
report="""# SLC16A3细胞来源与患者内共表达核查

## “实测乳酸”是公共组织质谱测量，不是RNA预测或细胞间通量

这里的“实测”指复用公共研究已经完成的质谱代谢组检测。RC18与RC20来自肾切除患者的新鲜冷冻肿瘤/邻近组织配对RNA和质谱数据；CPTAC及其验证队列也提供配对RNA和质谱峰计数。Phase160仅选肿瘤样本，人数分别为97、52、50、56，共255。我们使用发布的参考测量矩阵，没有运行UnitedMet去预测这些患者的乳酸。

检测值是组织中乳酸的半定量相对丰度，适合在同队列内比较高低，不是统一的mmol/L绝对浓度。组织乳酸水平、局部细胞外乳酸水平、细胞间乳酸转运及净生成通量属于不同测量对象。因此，Phase160中乳酸关联较弱意味着该组织指标未形成清楚桥接，不能据此排除局部转运；也不能反过来用SLC16A3的RNA关联证明乳酸释放。

原始来源说明及处理方法见[UnitedMet原始研究Methods](https://www.nature.com/articles/s43018-025-00943-0)；具体参考文件见[Zenodo 11286535](https://zenodo.org/records/11286535)。研究方法明确区分具有实测数据的参考队列与仅有RNA、需要预测代谢物的目标队列，本项目本轮使用前者。

## 结论：SLC16A3—SPP1巨噬状态联系得到支持，上游恶性TYMP联系仍弱

"""+summary.split('本轮直接回查',1)[1].join(['本轮直接回查',''])+"""

## Methods

使用冻结H5AD原始counts，按Phase154同患者同肿瘤文库匹配，选择16患者23018个恶性细胞与成熟巨噬细胞；C4成员沿用冻结定义。目标基因为TYMP、SLC16A3、LDHA、MSMO1、SPP1、CYP27A1、BSG和SLC16A1。逐细胞表达为log1p(CP10K)，区室表达为原始计数伪bulk log2(CPM+1)。

在各患者各区室内计算细胞表达的秩相关。预定要求每区室至少50细胞，且相关两基因各至少10细胞检出。主要模型控制总计数深度，原始及追加缺氧为敏感性。缺氧分数为冻结Hallmark Hypoxia基因在相应患者区室内标准化表达的均值，排除全部目标基因。

细胞间相关系数作为患者内描述性效应；不将细胞数当作统计样本量，不输出细胞级推断P。随后以独立患者的rho分布进行Wilcoxon符号秩检验，报告中位效应、正负患者数、P及FDR。FDR按区室和调整口径，在同类两条预列关联内计算。C4与其他成熟巨噬的相关差使用相同患者配对，属于后续探索比较；两基因对分别在相应调整家族校正。

恶性与成熟巨噬伪bulk表达差以16患者配对检验，8个基因单独构成来源比较家族。跨区室患者相关单列，区分原始和研究调整，其本轮三项协调家族FDR不替代此前Phase159的候选发现FDR。

## Results

### SLC16A3在两区室均可检出，恶性区室伪bulk表达更高

SLC16A3在16/16患者恶性区室的标准化伪bulk表达高于成熟巨噬，中位差+1.205 log2CPM，P=3.05e-5、FDR=6.10e-5。患者检出率中位数恶性区室70.2%、巨噬区室50.1%。BSG、SLC16A1、LDHA在恶性区室也较高。

TYMP则在16/16患者的成熟巨噬伪bulk表达高于恶性区室。SPP1在恶性细胞中同样广泛表达，13/16患者恶性伪bulk表达更高，中位差+2.009 log2CPM，P=0.000427、FDR=0.000488。该发现改变的是外部bulk来源解释：不能将组织SPP1全部标成巨噬来源；不改变已明确分区的巨噬SPP1结果或C4占巨噬SPP1计数的贡献。

### SLC16A3—SPP1在巨噬细胞内呈患者间一致的正向共表达

控制深度后，16患者中14正、2负，患者内rho中位0.108，患者级P=0.00214、FDR=0.00214；追加缺氧后仍14正、2负，中位rho=0.103，P=0.00269、FDR=0.00269。

LDHA—SPP1在深度与缺氧调整后同样14/16正，中位rho=0.0959，P=0.000305、FDR=0.000610。效应幅度不大，但患者方向一致，支持巨噬细胞内的代谢相关SPP1状态联系。

### C4中的SLC16A3—SPP1共表达提供额外的探索性定位

C4内9患者满足细胞数及检出标准，深度和缺氧调整后9/9正，中位rho=0.243，患者级P=0.00391、FDR=0.00391；按小于10患者规则仍作探索性结果。

与相同9患者其他成熟巨噬比较，8/9患者的C4相关更高，配对rho差中位+0.131，P=0.0273、FDR=0.0547。只控制深度时7/9更高、差+0.149、P=0.0391、FDR=0.0781。两者均为名义/探索性比较，不能称校正后显著的C4特异性。

LDHA—SPP1在C4也为9/9正，但C4相对其他巨噬的差异在深度与缺氧调整后中位仅+0.00583，P=0.570。该阴性比较保留，说明不能把所有糖酵解关联笼统称为C4特异增强。

### 恶性细胞TYMP—SLC16A3具有弱方向性，上游调控尚未确认

15患者满足TYMP与SLC16A3检出要求。原始模型12/15正，中位rho=0.0788，P=0.00153、FDR=0.00305；控制测序深度后中位rho=0.0652，仍12/15正，P=0.0833、FDR=0.1665；进一步控制缺氧后11/15正，中位rho=0.0471，P=0.1354、FDR=0.2708。

这组结果保留了多数患者的微弱正向对应，但表明原始显著性有技术深度和状态共享成分，不能仅挑原始模型作为恶性TYMP直接调控SLC16A3的证据。TYMP—MSMO1在12患者中的原始微弱关联同样在深度调整后减弱（中位rho=0.00815，P=0.301）。

### 跨区室协调与细胞内共表达分别构成证据

16患者恶性SLC16A3与巨噬SLC16A3伪bulk相关rho=0.526、P=0.0362、本轮协调家族FDR=0.0543；研究调整后rho=0.501、P=0.0815。恶性SLC16A3与巨噬SPP1维持此前rho=0.594、P=0.0152；本轮协调家族FDR=0.0457不能用于替代Phase159原候选发现FDR=0.0656，两者均保留。研究调整后rho=0.519、P=0.0693。

这些结果与巨噬细胞内共表达是不同层次，不代表恶性细胞已经向巨噬细胞传递了某种确定信号。跨患者恶性TYMP—恶性SLC16A3仍为rho=0.300、P=0.259，和上一阶段一致。

## Discussion与后续定位

本轮区室来源比较只覆盖已有恶性细胞及成熟巨噬细胞，不是全组织所有细胞类型的来源清单。伪bulk表达高低不等于蛋白分泌量、细胞外浓度或代谢通量。外部bulk中SPP1可有恶性来源，必须避免把组织层SLC16A3—SPP1结果直接写成巨噬特异机制。

患者内共表达、跨患者亚群平均表达和跨区室患者配对相关回答不同问题。本轮C4内9患者的共表达正向，不推翻Phase158中12患者跨患者C4平均SPP1与肿瘤候选缺少关联的结果。C4由表达数据定义，本次是冻结状态内的探索核查，尚不是外部C4独立验证。细胞覆盖QC导致各基因对患者数不同，全部排除与计数保留，未为增强关联调整阈值。

SLC16A3继续作为优先补充候选，其巨噬SPP1关联已有组织层、跨区室和患者内共表达对应；主要缺口仍为恶性TYMP上游与具体跨细胞作用方式。不能把局部转运尚未测到误写成已由实测乳酸证实，也不因组织乳酸相关弱直接排除所有局部作用。

维持用户确认的分工：整体巨噬承载SPP1功能联系，C4提供SPP1信号集中和胆固醇处理状态定位。MTA实验相关假说与客户图背景保留。本轮没有客户实验数据进入统计、没有生成新图。

## 审计及完整附表

原始H5AD重新抽取与Phase154/157缓存的区室表达一致；细胞内调整相关由独立精度矩阵公式核对。615/615计算检查通过。细胞级输入、覆盖QC、患者效应、完整P/FDR、来源文件和脚本全部服务器归档。计算一致性不是外部生物学验证。

"""
report+='\n\n### 患者级共表达全表\n\n'+md(S)
report+='\n\n### 同患者C4—其他巨噬效应差\n\n'+md(D)
report+='\n\n### 来源比较\n\n'+md(SRC)
report+='\n\n### 跨区室患者协调\n\n'+md(X)
report+='\n\n### 完整患者覆盖与排除\n\n'+md(QC)
report+='\n\n### 完整患者内效应\n\n'+md(R)
(O/'SLC16A3_CELL_ORIGIN_AND_LACTATE_CN.md').write_text(report,encoding='utf-8')
status=json.loads((O/'STATUS.json').read_text());assert status['total_checks']==615
status.update(status='COMPLETE',report=str(O/'SLC16A3_CELL_ORIGIN_AND_LACTATE_CN.md'),confirmed_mediator=False,C4_retained=True,MTA_retained=True,completion_time=datetime.datetime.now(datetime.timezone.utc).isoformat())
(O/'STATUS.json').write_text(json.dumps(status,indent=2))
note="""## 2026-09-15 Phase161 SLC16A3细胞来源与实测乳酸解释完成

用户授权下一步并询问实测乳酸。已解释为RC18/RC20/CPTAC/CPTAC_val公共肿瘤组织质谱相对丰度，共255肿瘤患者配对RNA；使用已测参考矩阵，不是RNA预测、绝对浓度或细胞间通量。
同16文库23018细胞直接原始counts核查：SLC16A3恶性区室16/16高于成熟巨噬，差+1.205log2CPM/P3.05e-5/FDR6.10e-5，两区室均表达。SPP1恶性区室13/16更高，外部bulkSPP1不能全部归为巨噬。
患者内共表达以患者为单位：深度+缺氧调整的巨噬SLC16A3—SPP1为14/16正、中位rho0.103/P0.00269/FDR0.00269；C4为9/9正、中位rho0.243，小样本探索。与同9患者其他巨噬配对8/9 C4更高，差0.131/P0.0273/FDR0.0547名义支持；LDHA相关差不支持C4特异增强。
恶性TYMP—SLC16A3深度调整后15患者12正，中位rho0.065/P0.0833/FDR0.1665；再加缺氧11正、0.047/P0.1354/FDR0.2708。保留弱方向，上游未确认。跨患者恶性TYMP—MCT4仍rho0.3/P0.259；恶性MCT4—巨噬SPP1rho0.594/P0.0152与旧结果一致，不能用本轮小家族FDR替换原发现FDR。
报告phase161_SLC16A3_cell_origin_20260915_v1/SLC16A3_CELL_ORIGIN_AND_LACTATE_CN.md；615项检查PASS，无图无客户数据入模。C4与整体巨噬双层分工、MTA客户实验背景保留。主要缺口恶性TYMP上游及具体跨细胞作用方式，暂无确认中介。

"""
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 p=P/name;shutil.copy2(p,O/(name+'.before'));p.write_text(note+p.read_text(encoding='utf-8'),encoding='utf-8')
with (P/'RESULTS_MANIFEST.tsv').open('a',encoding='utf-8') as f:f.write('\t'.join(['phase161_SLC16A3_cell_origin_20260915_v1','Phase161','report',str(O/'SLC16A3_CELL_ORIGIN_AND_LACTATE_CN.md'),status['completion_time'],'PASS','615 computational checks','16 patients; 23018 cells; patient-unit inference; no figures'])+'\n')
p=P/'part2_rebuild_v1/STATUS.json';shutil.copy2(p,O/'PART2_STATUS.before.json');b=json.loads(p.read_text());b.update(latest_phase='Phase161',latest_report=str(O/'SLC16A3_CELL_ORIGIN_AND_LACTATE_CN.md'),stage='SLC16A3_cell_origin_complete',confirmed_mediator=False,C4_retained=True,MTA_retained=True);p.write_text(json.dumps(b,ensure_ascii=False,indent=2))
manifest={x:hashlib.sha256((O/x).read_bytes()).hexdigest() for x in ['SUMMARY_CN.md','SLC16A3_CELL_ORIGIN_AND_LACTATE_CN.md']}
(O/'DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(manifest));print('COMPLETE')

