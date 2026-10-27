from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,datetime,shutil
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
O=P/'phase158_three_candidate_followup_20260915_v1'
S=P/'phase157_ordered_bridge_reassessment_20260915_v1'
R=pd.read_csv(O/'THREE_CANDIDATE_ENDPOINT_RESULTS.tsv',sep='\t')
T=pd.read_csv(O/'HYPOXIA_PROLIFERATION_SENSITIVITY.tsv',sep='\t')
E=pd.read_csv(O/'COMPARTMENT_SPP1_SECONDARY.tsv',sep='\t')
L=pd.read_csv(O/'LEAVE_ONE_RESULTS.tsv',sep='\t')
names={'ligand:IGFBP3':'IGFBP3','program:Cholesterol Homeostasis':'胆固醇稳态','program:Glycolysis':'糖酵解'}
def md(df):
 headers=list(df.columns);out=['| '+' | '.join(headers)+' |','|'+'|'.join(['---']*len(headers))+'|']
 for row in df.itertuples(index=False,name=None):
  out.append('| '+' | '.join(f'{x:.6g}' if isinstance(x,(float,np.floating)) else str(x) for x in row)+' |')
 return '\n'.join(out)
summary="""# 三个候选复核：胆固醇稳态与糖酵解优先用于TYMP—SPP1主线

已完成IGFBP3、胆固醇稳态、糖酵解的同16患者复核。当前优先顺序为胆固醇稳态、糖酵解、IGFBP3；这是主线适配优先级，不是已确认中介的排序。

1. 胆固醇稳态：TYMP—程序rho=0.591，P=0.0159，FDR=0.1259；程序—C4 rho=0.585，P=0.0172，FDR=0.4481。补充SPP1终点后，控制研究和TYMP的条件关联为rho=0.591，P=0.0431，FDR=0.5388；SPP1检出比例对应rho=0.615，P=0.0334，FDR=0.3342。两端有方向性对应，SPP1的名义信号出现在调整模型中；原始SPP1相关P=0.113。
2. 糖酵解：TYMP—程序rho=0.421，P=0.1048；程序—C4 rho=0.565，P=0.0227，FDR=0.4481。控制研究和TYMP后，程序与SPP1检出比例rho=0.580，P=0.0480，FDR=0.3996；SPP1表达条件P=0.106。缺氧调整使TYMP—糖酵解关联减弱。可作为第二条程序线索。
3. IGFBP3：与C4比例rho=0.715，P=0.00186，FDR=0.0348；但与整体巨噬SPP1表达rho=0.062，P=0.820、FDR=0.911，与SPP1检出比例rho=0.124，P=0.649、FDR=0.938。当前更适合作为C4丰度相关候选，不支持将它写成已连接TYMP与SPP1的分子桥。

12例C4细胞数>=20的分层分析中，三个候选与C4内部SPP1表达均未获得明确关联。C4数量变化与C4内部SPP1表达变化必须分开。

外部原始研究为乳酸—SPP1巨噬和氧固醇—SPP1提供了机制先验，但本项目尚未测得这两种代谢物，也未证明TYMP因果调控。MTA继续保留。1046项计算复核通过；未生成新图，未使用客户湿实验数据。
"""
(O/'SUMMARY_CN.md').write_text(summary,encoding='utf-8')
report="""# IGFBP3、胆固醇稳态和糖酵解对TYMP—C4/SPP1主线的复核

## 结论：两个肿瘤代谢程序更贴近TYMP—SPP1主线，IGFBP3主要对应C4丰度

"""+summary.split('\n\n',1)[1]+"""
## Methods

沿用Phase154同一肿瘤文库配对的16名独立患者；恶性细胞TYMP及候选X配对成熟巨噬细胞终点。C4冻结成员不变，C4比例保留零值。X沿用Phase157的恶性细胞IGFBP3伪bulk表达、Hallmark胆固醇稳态及糖酵解程序；程序排除TYMP和SPP1。主要终点为C4/成熟巨噬比例，补充终点为成熟巨噬SPP1伪bulk log2CPM与SPP1阳性细胞比例。

各变量取标准化秩。a为TYMP—X相关，XY为X—终点相关，b为控制TYMP的X—终点偏相关；研究调整使用研究指示变量。偏相关P采用有效残差自由度的t近似。每次仅追加一个缺氧或G2-M程序作C4敏感性分析。留一患者16次、留一研究4次均按原始模型重新计算。C4及其他成熟巨噬的SPP1分层分析采用同文库内各区室至少20细胞标准；C4为12患者、其他成熟巨噬为16患者。

新增SPP1检验对原116候选完整计算，按56配体、50程序、10代谢/转运基因的原家族分别计算BH FDR，并区分终点、边和研究调整口径；没有用筛选后的三个候选重置发现FDR。附加混杂/区室分析的三候选FDR只适用于相应后续探索家族，表中明确标注，不能代替原筛查FDR。

## Results

### 胆固醇稳态同时连接TYMP、C4丰度及调整后的SPP1指标

本地公共数据中，TYMP—胆固醇程序的正向关联在研究调整后仍存在（rho=0.645，P=0.0173；原家族FDR见表）。控制G2-M后保持正向；控制缺氧后效应减小（rho=0.475，P=0.0733；同时研究调整rho=0.524，P=0.0800）。该结果提示它与增殖并非完全重合，但与缺氧存在共享信息。

程序—C4原始边际关联为正，控制TYMP后的b为rho=0.433，P=0.1067，FDR=0.8720。研究调整后b为rho=0.332，P=0.2911。补充SPP1发现，研究和TYMP调整后的SPP1表达及检出比例均有名义正向信号，见完整统计表。原始SPP1表达关联未达名义显著，不能只挑调整模型作为唯一结果。

Phase157中C4间接效应a*b的16次留一患者均为正；去除GSE207493后余6患者的a*b接近零且略负（-0.00238）。这条线适合作为优先程序候选。

### 糖酵解具有SPP1检出比例的调整后正向线索，上游关联与缺氧共享

糖酵解—C4原始边际rho=0.565，P=0.0227，FDR=0.4481；控制TYMP后的b为rho=0.459，P=0.0853，FDR=0.8720。研究与TYMP调整后，糖酵解—SPP1检出比例rho=0.580，P=0.0480，FDR=0.3996；整体SPP1表达条件rho=0.489，P=0.1065。

TYMP—糖酵解原始rho=0.421，P=0.1048；缺氧调整后降至rho=0.159，P=0.5719。因此，现有数据没有把TYMP自身的作用与共同缺氧背景分开。

Phase157中16次留一患者C4间接效应均为正；4次留一研究的乘积也均为正，但去除GSE207493后仅6患者，a和b同时转负。负乘负不能用来声称正向路径跨研究稳定。

### IGFBP3与C4丰度相关，但未连接到SPP1表达

IGFBP3—C4边际关联通过原56配体家族FDR（rho=0.715，P=0.00186，FDR=0.0348）。控制TYMP后的b为rho=0.658，P=0.00769，FDR=0.1077，为名义显著。分别控制缺氧或G2-M后，原始研究混合模型中的条件关联保持；研究调整后b下降至rho=0.351，P=0.2628。

上游TYMP—IGFBP3仅为正向趋势（rho=0.376，P=0.1506），研究调整后rho=0.153，P=0.6176。新增终点中，IGFBP3与成熟巨噬SPP1表达、SPP1检出比例及C4内部SPP1均缺少对应。当前数据不支持把IGFBP3写成SPP1表达上调的中介。

Phase157的TMEM219在24/24名符合巨噬区室覆盖条件的患者C4中达到表达阈值；但C4相对其他成熟巨噬的受体表达差异中位数为-0.108 log2CPM，未显示C4特异富集。接收端可用性支持进一步考察IGFBP3，却不规定作用方向。

### C4内部SPP1表达没有随三个候选呈现一致上升

同文库C4>=20细胞的12患者中，IGFBP3、胆固醇稳态、糖酵解与C4内部SPP1表达的边际rho依次为-0.133、-0.014、-0.098，P依次为0.681、0.966、0.762。其他成熟巨噬中的两个代谢程序则为正向趋势，详情见附表。整体C4比例关联不能自动解释为C4单细胞SPP1表达增强。

## 外部原始研究核查

- 乳酸线：脑肿瘤原始研究报告缺氧适应肿瘤细胞通过乳酸介导的组蛋白乳酰化促进SPP1阳性巨噬状态，为“肿瘤糖酵解相关代谢输出—巨噬SPP1”提供机制先验。[原始研究，PMID 40973181](https://pubmed.ncbi.nlm.nih.gov/40973181/)
- 氧固醇线：宫颈腺癌原始研究报告27-HC处理提高THP-1来源巨噬的SPP1转录及分泌，并考察肿瘤DHCR7干预后的条件培养基效应。[原始研究，Frontiers in Immunology 2026](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2026.1808107/full)
- IGFBP3线：受体鉴定研究在乳腺/前列腺癌体系中支持IGFBP3受体介导caspase-8依赖凋亡；这不构成促进SPP1巨噬的方向证据。[原始受体研究，PMID 20353938](https://pubmed.ncbi.nlm.nih.gov/20353938/)

## Discussion与下一步定位

三个候选均来自同一批患者的筛选后复核，新增SPP1、混杂和区室分析属于后续探索，不是独立患者验证。原始与研究调整模型呈现不同强度，不能按P值择优定义主模型。n<10的留研究结果只作方向性探索；16例模型也不宜不断叠加协变量。

胆固醇稳态和糖酵解是RNA程序，并非胆固醇、27-HC、乳酸浓度或代谢流。外部癌种的SPP1巨噬不等同本项目冻结C4。文献提供的是候选机制，不能补写当前数据未测到的中间物，也不能证明TYMP位于因果上游。IGFBP3受体表达不等于信号激活，受体已知促凋亡功能也不能直接外推至本项目巨噬。

按根本问题的适配度，优先胆固醇稳态，糖酵解作为并行候选，IGFBP3保留在C4丰度关联分支。后续可将两个程序拆为合成、输出、接收和响应环节，先核对本项目表达覆盖与两端方向，再决定具体X；不把27-HC或乳酸预设为唯一答案。MTA维持候选状态，本次结果不推翻此前实测MTA证据。

## 完整统计附表

下表的rho为秩相关或条件秩相关；全部效应、P及相应FDR保留。补充表文件列出各患者/研究的留一结果。
"""
cols=['candidate','endpoint','adjustment','n_patients','a_rho','a_P','a_FDR','XY_rho','XY_P','XY_FDR','b_partial_rho','b_P','b_FDR']
report+='\n\n### 三候选全部终点\n\n'+md(R[cols])+'\n\n### 缺氧与增殖敏感性\n\n'+md(T)+'\n\n### 区室SPP1表达\n\n'+md(E)+'\n\n### 留一稳定性摘要\n\n'
lr=[]
for (candidate,endpoint,kind),g in L.groupby(['candidate','endpoint','kind']):
 lr.append(dict(candidate=candidate,endpoint=endpoint,kind=kind,fits=len(g),a_positive=int((g.a_rho>0).sum()),b_positive=int((g.b_partial_rho>0).sum()),ab_positive=int((g.indirect_ab>0).sum()),XY_min=g.XY_rho.min(),XY_max=g.XY_rho.max()))
report+=md(pd.DataFrame(lr))
report+="""

## 审计与文件

1046/1046计算复核通过，包括全部原始模型X—终点相关与scipy独立直接计算的一致性及路径分解恒等式。此项为计算核查，不代表外部生物学验证。输入SHA256、脚本、内存门槛、完整116候选终点结果、三候选敏感性和患者明细均在服务器本版本目录。未生成新图，未使用客户湿实验；旧结果未覆盖。
"""
(O/'THREE_CANDIDATE_REVIEW_CN.md').write_text(report,encoding='utf-8')
status=json.loads((O/'STATUS.json').read_text());status.update(status='COMPLETE',report=str(O/'THREE_CANDIDATE_REVIEW_CN.md'),priority=['Cholesterol Homeostasis','Glycolysis','IGFBP3_C4_branch'],confirmed_mediator=False,figures_created=False,client_data_used=False)
(O/'STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
note="""## 2026-09-15 Phase158 三候选C4/SPP1复核完成

用户“这三个看一下”按IGFBP3、胆固醇稳态、糖酵解执行。权威报告phase158_three_candidate_followup_20260915_v1/THREE_CANDIDATE_REVIEW_CN.md。
16同文库独立患者，保留C4冻结；新增SPP1终点仍以116原候选家族计算FDR。IGFBP3-C4边际rho0.715/P0.00186/FDR0.0348，但整体SPP1rho0.062/P0.820；保留C4丰度分支。胆固醇程序TYMP关联rho0.591/P0.0159；研究+TYMP调整后SPP1表达rho0.591/P0.0431/FDR0.5388，检出rho0.615/P0.0334/FDR0.3342，原始SPP1P0.113。糖酵解调整后SPP1检出rho0.580/P0.0480/FDR0.3996，但TYMP上游缺氧调整减弱。12患者C4内部SPP1未见三个候选明确关联。主线优先胆固醇稳态、并行糖酵解；无确认中介。外部原始文献只作机制先验，非肾癌C4独立验证。1046检查通过，无图无客户数据。
后续方向：拆分两个程序的合成/输出/接收/响应，保留MTA，不预设特定代谢物。

"""
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 p=P/name
 backup=O/(name+'.before');shutil.copy2(p,backup)
 p.write_text(note+p.read_text(encoding='utf-8'),encoding='utf-8')
with (P/'RESULTS_MANIFEST.tsv').open('a',encoding='utf-8') as f:
 f.write('phase158_three_candidate_followup_20260915_v1\t'+str(O/'THREE_CANDIDATE_REVIEW_CN.md')+'\tCOMPLETE\t16 patients; 1046 checks; no figures; no confirmed mediator\n')
branch=P/'part2_rebuild_v1/STATUS.json'
shutil.copy2(branch,O/'PART2_STATUS.before.json')
b=json.loads(branch.read_text());b.update(latest_phase='Phase158',latest_report=str(O/'THREE_CANDIDATE_REVIEW_CN.md'),stage='three_candidate_followup_complete',confirmed_mediator=False)
branch.write_text(json.dumps(b,ensure_ascii=False,indent=2))
manifest={name:hashlib.sha256((O/name).read_bytes()).hexdigest() for name in ['SUMMARY_CN.md','THREE_CANDIDATE_REVIEW_CN.md']}
(O/'DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(manifest));print('REPORT_COMPLETE')

