from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,datetime,shutil
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=P/'phase160_external_bridge_investigation_20260915_v1'
R=pd.read_csv(O/'COHORT_BRIDGE_RESULTS.tsv',sep='\t');M=pd.read_csv(O/'META_BRIDGE_RESULTS.tsv',sep='\t');OM=pd.read_csv(O/'META_DIRECT_AND_GENE_METABOLITE.tsv',sep='\t');L=pd.read_csv(O/'LEAVE_ONE_COHORT_META.tsv',sep='\t')
client=dict(source='User supplied image, 2026-09-15',local_path='C:/Users/31937/Downloads/湿实验.png',sha256='aabf51942981fbf67cf578259072d4ca9ba84f2efcfde2f3855c087c82fff123',permission='User explicitly permits reference; raw data unavailable',model='786-O shTYMP versus O-NC',reported_directions={'5-prime-MTA':'decreased','Methionine':'decreased','SAM':'decreased','Thymine':'decreased','Thymidine':'increased','Spermidine':'decreased'},use='Experimental background and narrative correspondence only; not reanalyzed or included as public validation, no inferred numerical effects/p-values',not_shown='Direct MTA intervention in SPP1-positive macrophages is not shown in the supplied panels')
(O/'CLIENT_FIGURE_REFERENCE.json').write_text(json.dumps(client,ensure_ascii=False,indent=2))
agreement='\n\n## 2026-09-15 客户湿实验图参考口径再次确认\n用户提供Downloads/湿实验.png，明确可作参考、没有原始数据。允许作为客户实验背景和叙事对应：786-O shTYMP后MTA/甲硫氨酸/SAM/Thymine/Spermidine下降、Thymidine上升。无需原始数据即可继续公共调查；不反推效应量、P值或44倍，不将客户图计入公共验证。当前公共模型均不使用该图数据。MTA继续保留实验相关候选，SLC16A3作为公共数据支持的补充候选，不自动替换客户假说。图中未展示的MTA对SPP1巨噬功能验证不写成已完成。\n'
p=P/'AGENTS.md'
if '客户湿实验图参考口径再次确认' not in p.read_text():p.write_text(p.read_text()+agreement,encoding='utf-8')
def md(df):
 out=['| '+' | '.join(df.columns)+' |','|'+'|'.join(['---']*len(df.columns))+'|']
 for row in df.itertuples(index=False,name=None):out.append('| '+' | '.join(f'{x:.6g}' if isinstance(x,(float,np.floating)) else str(x) for x in row)+' |')
 return '\n'.join(out)
summary="""# 公共队列将SLC16A3确定为优先追踪候选，MTA保留实验主线背景

客户提供的786-O shTYMP图可作参考，不需要原始数据才能继续公共调查。按图，TYMP敲低后MTA、Methionine、SAM、Thymine、Spermidine下降，Thymidine上升。该图作为客户报告的实验背景，不重算效应、P值，不计入公共数据验证；图中未展示直接MTA干预SPP1巨噬细胞的功能验证。

本轮完成RC18、RC20、CPTAC、CPTAC_val四个外部配对RNA–实测代谢组队列，肿瘤患者分别97、52、50、56，共255。使用参考队列实测数据，没有使用UnitedMet预测代谢值。

**SLC16A3（MCT4）成为优先补充候选。** 控制髓系RNA代理和缺氧后，四队列的TYMP—SLC16A3均正向，随机效应合并rho=0.329、P=2.02e-5、FDR=8.06e-5；控制TYMP及上述代理后，SLC16A3—SPP1四队列均正向，合并rho=0.357、P=1.09e-4、FDR=2.18e-4。四次留一队列，两端合并方向均保持正向。原始RC20的TYMP—SLC16A3略负，且完全调整后CPTAC单队列上游很弱，全部保留。不是每个队列都显著，也不是已确认因果中介。

**实测乳酸不支持直接把MCT4解释成乳酸升高。** 原始随机效应：TYMP—乳酸rho=0.090、P=0.477；乳酸—SPP1rho=0.102、P=0.207；MCT4—乳酸rho=0.068、P=0.289。CPTAC存在两端名义正向，但其他队列不一致。胆固醇仅RC18/RC20可测，共149患者，两端同样未形成清楚联系；27-HC未测。

**MTA继续保留。** 同一批公共数据重现已有MTA—组织SPP1正向关系（rho=0.139、P=0.0293，旧FDR=0.0975保留），不是新增独立验证；TYMP—实测MTA合并仍弱。客户实验与公共组织相关回答不同问题，不要求公共数据复刻干预幅度，也不借客户结果弥补公共上游统计。

MSMO1—组织SPP1在外部队列正向，但其与TYMP的上游方向不一致；LDHA—SPP1稳定正向，上游TYMP弱。CYP27A1更适合继续保留在C4胆固醇处理表型与候选响应分支。

已确认的分工保持：整体巨噬承载SPP1代谢联系；C4定位SPP1信号集中的胆固醇处理状态。外部bulk不能直接验证恶性细胞来源或C4特异机制。741项复核通过，无新图；客户图只作背景。
"""
(O/'SUMMARY_CN.md').write_text(summary,encoding='utf-8')
report="""# TYMP—中间环节—SPP1的外部公共队列调查

## 结论：SLC16A3获得跨队列两端支持，MTA保留实验相关机制假说

"""+summary.split('\n\n',1)[1]+"""

## 客户实验背景与公共证据分工

用户于2026-09-15提供湿实验.png并授权作参考。图示786-O细胞TYMP敲低伴随MTA、甲硫氨酸、SAM、胸腺嘧啶、亚精胺下降，胸苷上升，支持继续调查TYMP相关核苷—甲硫氨酸/SAM—多胺代谢变化。上述方向均归属客户报告，不是本项目重新计算的结果。

图中关于后续SPP1阳性巨噬细胞生态位的文字，是机制研究方向；所给图板未展示直接MTA处理巨噬细胞的结果。本轮不索要客户原始数据，不重做客户分析，不从坐标或星号倒算效应量、P值、样本数或44倍变化。客户图不进入候选选择后的公共统计，也不增加公共样本量。

| 联系 | 客户图提供的信息 | 当前公共数据提供的信息 |
|---|---|---|
| TYMP干预—MTA代谢变化 | shTYMP后MTA下降及相关代谢物同步变化 | 组织TYMP—MTA总体相关弱，保留队列异质性 |
| MTA—SPP1巨噬 | 图板中未展示直接功能验证 | 已有组织MTA—SPP1弱正向关系；不是巨噬特异实测 |
| TYMP—SLC16A3—SPP1 | 客户图未直接考察这条路径 | 本轮外部组织层两端相关支持，作为补充候选 |
| C4细胞状态定位 | 客户图不提供C4分群结果 | 单细胞C4的SPP1及胆固醇处理表型保留 |

## Methods

数据来自UnitedMet公开参考队列：RC18、RC20、CPTAC、CPTAC_val，使用其配对原始参考表达和已测代谢物文件，不使用预测代谢矩阵。按MasterMapping_ccRCC.csv匹配RNAID与MetabID，限Tumor，映射为患者标识。相同患者多个肿瘤样本的预定处理为中位数聚合；实际核心分析255条均为不同映射患者，四队列之间未发现重复映射标识。RNA使用log2(TPM+1)，代谢物保留已发布参考丰度，相关分析在各队列内取秩。不插补缺失，按每个检验完整配对取样，至少10患者纳入推断。

预列主要基因候选MSMO1、LDHA、SLC16A3、CYP27A1；HK2和LDHB为背景基因。主要新代谢物为乳酸、胆固醇；MTA为保留的既有候选，另记录丙酮酸、desmosterol、7-hydroxycholesterol（alpha或beta未区分）、胆固醇硫酸酯。未测到的指标保留为未测，7-hydroxycholesterol不作为27-HC替代。

髓系RNA代理由LST1、AIF1、CSF1R、CD68、FCER1G的队列内z表达均值构成；缺氧代理由已冻结Hallmark Hypoxia基因集构成，去除TYMP、SPP1、全部候选与髓系代理基因。分别报告原始、髓系代理调整、髓系与缺氧联合调整。a为TYMP—X，XY为X—SPP1，b为进一步控制TYMP的X—SPP1偏相关。标准化秩回归计算路径系数，a*b仅保留作描述，不计算因果中介比例。偏相关P使用相应有效残差自由度的t近似。

各队列相关经Fisher z合并，方差近似为1/(n-k-3)，k为已控制变量数。固定效应及DerSimonian–Laird随机效应均保存，随机效应为本轮跨队列主要合并口径，采用正态近似区间，同时报告Q、I²与tau²。仅单队列可测的指标作单队列结果标记。四队列均可测的候选补充留一队列分析。

BH FDR分别在主要基因4项、背景基因2项、新主要代谢物2项、既有MTA背景项、其他代谢物以及基因—代谢物联系中计算，并区分边、调整方式、单队列与合并层。预设TYMP—SPP1单一联系不需要多重校正。表中MTA单项FDR等于P；既有MTA边际证据的历史FDR完整保留，不据此重定义其既有证据等级。

## Results

### SLC16A3在外部组织数据中获得两端方向一致的支持

原始随机效应合并TYMP—SLC16A3 rho=0.304、P=0.0166、FDR=0.0665，I²=75.0%，为名义正向且异质性较大。SLC16A3—SPP1边际rho=0.414、P=6.70e-12、FDR=2.68e-11；控制TYMP后rho=0.392、P=5.10e-8、FDR=1.02e-7。

髓系RNA代理调整后，TYMP—SLC16A3四队列均正向，合并rho=0.394、P=1.15e-10、FDR=4.60e-10，I²=0。追加缺氧后rho=0.329、P=2.02e-5、FDR=8.06e-5，I²=31.5%；相应条件b合并rho=0.357、P=1.09e-4、FDR=2.18e-4，I²=51.5%。四次留一队列，完全调整的a范围0.297–0.392、b范围0.277–0.411，全部正向。

关键异质性保留：RC20原始a为-0.055、P=0.698，调整后转正；完全调整后CPTAC的a仅0.068、P=0.645。完全调整的b在RC18与RC20较强，CPTAC与CPTAC_val单队列未达名义显著。因此，结论是跨队列合并和方向一致性支持，不是每个队列各自证明整条路径。

此前同文库单细胞16患者中，恶性SLC16A3与恶性TYMP仅rho=0.300、P=0.259；恶性SLC16A3与巨噬SPP1为名义正向。外部结果增加了组织层两端证据，但未直接替代单细胞恶性来源上游的弱关联。

### 实测乳酸没有重现同样强度的TYMP—X—SPP1联系

乳酸四队列255患者均可测。原始TYMP—乳酸合并rho=0.090、P=0.477、FDR=0.477，I²=73.5%；乳酸—SPP1边际rho=0.102、P=0.207、FDR=0.415，条件b rho=0.083、P=0.197、FDR=0.385。SLC16A3—实测乳酸rho=0.068、P=0.289、FDR=0.404。

CPTAC的TYMP—乳酸rho=0.399、P=0.00411，乳酸—SPP1rho=0.361、P=0.0100，保留为该队列名义正向；RC18、RC20的TYMP—乳酸原始方向为负，CPTAC_val的乳酸—SPP1接近零。联合髓系/缺氧调整后a合并呈名义正向（rho=0.160、P=0.0488、FDR=0.0975），但b仍弱（rho=0.092、P=0.159、FDR=0.318）。

当前数据不支持把SLC16A3表达关联直接翻译为组织乳酸升高介导SPP1变化。LDHB—乳酸合并正向（rho=0.242、P=0.000530、FDR=0.00186）也不能替代TYMP上游缺口，或推断乳酸净生成通量。

### 总胆固醇与27-HC需要区分，C4表型不依赖总胆固醇相关成立

胆固醇仅RC18与RC20可测，共149患者。TYMP—胆固醇两个队列均为负，合并rho=-0.124、P=0.136、FDR=0.273；胆固醇—SPP1rho=0.046、P=0.578、FDR=0.578。联合代理调整仍没有清楚两端联系。

这些参考文件未测27-HC。7-hydroxycholesterol不能代替27-HC，也不能将总胆固醇不相关外推为所有氧固醇机制被排除。C4的CYP27A1、NR1H3、ABCA1等表达特征，以及SPP1原始RNA集中贡献，仍按Phase159结果保留。

MSMO1—desmosterol在RC18可用38患者中rho=0.574、P=0.000111、相应联系家族FDR=0.000777；这是单队列背景联系，不是外部多队列复制，也没有因此确认TYMP—MSMO1—desmosterol—SPP1链条。

### MTA与SPP1的已有弱正向关系保留，客户实验提供另一层上游背景

本轮原始MTA—SPP1合并rho=0.139、P=0.0293，四队列均正向；复现的是Phase98同一批公共参考数据，历史FDR=0.0975保留，不能称新增独立复现。原始条件b合并rho=0.132、P=0.0399。

TYMP—MTA原始随机效应rho=0.058、P=0.585，I²=62.6%；旧固定效应rho=0.0395、P=0.538已复算一致。髓系代理和缺氧调整后a为rho=0.220、P=0.0594，异质性仍高；条件b为rho=0.173、P=0.00798。这些后续调整结果完整保留为探索性支持，不据单个模型抬升MTA为已确定中介。

客户的TYMP敲低—MTA下降属于干预背景，公共组织TYMP—MTA相关属于跨患者观察。两者不要求效应幅度相同，也不能相互替代。MTA继续保留；本轮新增的SLC16A3证据作为并行补充候选。

### 其他基因未获得与SLC16A3同样完整的跨队列两端支持

MSMO1—SPP1四队列原始均正向，条件b合并rho=0.384、P=0.00251、FDR=0.00335；TYMP—MSMO1却方向不一致，合并rho=-0.020、P=0.882，I²=75.9%。这与单细胞恶性MSMO1—TYMP较强的结果构成区室/队列差异，不能隐藏。

LDHA—SPP1四队列原始均正向，条件b合并rho=0.463、P=1.36e-8、FDR=5.45e-8；控制髓系与缺氧后仍为rho=0.414、P=6.31e-6、FDR=2.52e-5，但上游TYMP—LDHA为rho=0.048、P=0.465。

CYP27A1—TYMP原始四队列均正向，合并rho=0.126、P=0.0491、FDR=0.0982；控制髓系代理后降至rho=0.027、P=0.717。CYP27A1—SPP1条件关联接近零。CYP27A1继续用于C4胆固醇处理表型和响应候选定位。

## Discussion与后续调查定位

本轮外部参考数据为整体组织RNA与代谢物测量，不能区分TYMP、SLC16A3、SPP1分别来自恶性、巨噬或其他细胞。髓系RNA代理不是病理纯度或真实巨噬比例；代理调整后的关联不能称完全排除了细胞混合。缺氧分数同样是RNA代理。公开样本标识已完成去重，但数据层的复核不能追踪未公开的跨研究身份。

四队列随机效应合并只有四个研究，区间与异质性估计采用近似；表中保留单队列、固定效应、随机效应和留一结果。单个队列不显著不抹除合并支持，跨队列合并显著也不等于每个队列均成立。样本同数据内的基因关联不是因果中介证据，不能用a*b点估计声称作用份额。

SLC16A3是目前最值得做细胞来源和作用环节核查的补充候选。下一步应优先确认恶性细胞内TYMP与SLC16A3的对应，以及SLC16A3相关巨噬SPP1状态；由于实测乳酸联系弱，不预设乳酸是唯一中间物或唯一作用方式。MTA相关实验背景及公共SPP1关联继续保留，核苷—甲硫氨酸/SAM—多胺方向不因本轮转运基因候选而被替换。

报告采用用户确认的双层分工：整体巨噬细胞承载SPP1代谢联系，C4提供SPP1信号集中的胆固醇处理状态定位。不得为了让机制看起来完整，将整体组织或整体巨噬结果标注成C4特异结果。

## 来源与归档

公共参考文件来自[UnitedMet ccRCC reference datasets，Zenodo 11286535](https://zenodo.org/records/11286535)。原始研究说明这些参考队列具有配对RNA和代谢测量，与模型预测的目标队列分开：[UnitedMet，Nature Cancer 2025](https://www.nature.com/articles/s43018-025-00943-0)。

客户图身份和使用边界保存在CLIENT_FIGURE_REFERENCE.json；图SHA256为aabf51942981fbf67cf578259072d4ca9ba84f2efcfde2f3855c087c82fff123。没有客户原始数据进入公共模型，没有新图。

741项计算核查通过：偏相关由独立精度矩阵公式核对，合并效应与FDR独立复算，并复现历史MTA固定效应结果。原始输入SHA256、患者映射、完整效应/P/FDR、各队列缺测指标和留一结果均服务器归档。

## 统计附表

"""
cols=['candidate','adjustment','edge','cohorts','total_n','positive_cohorts','fixed_rho','fixed_P','fixed_FDR','random_rho','random_low','random_high','random_P','random_FDR','I2']
report+='\n\n### 主要基因与实测代谢物合并结果\n\n'+md(M[M.family.isin(['gene_primary','met_primary','met_context'])][cols])
report+='\n\n### SLC16A3完整队列结果\n\n'+md(R[R.candidate=='SLC16A3'])
report+='\n\n### SLC16A3留一队列结果\n\n'+md(L[L.candidate=='SLC16A3'])
report+='\n\n### 实测乳酸、胆固醇与MTA单队列结果\n\n'+md(R[R.candidate.isin(['MET:lactate','MET:cholesterol','MET:5-methylthioadenosine (MTA)'])])
report+='\n\n### 直接联系与基因—实测代谢物\n\n'+md(OM)
(O/'EXTERNAL_BRIDGE_AND_CLIENT_CONTEXT_CN.md').write_text(report,encoding='utf-8')
status=json.loads((O/'STATUS.json').read_text());assert status['total_checks']==741
status.update(status='COMPLETE',report=str(O/'EXTERNAL_BRIDGE_AND_CLIENT_CONTEXT_CN.md'),priority_candidate='SLC16A3',confirmed_mediator=False,client_reference_used=True,client_data_in_public_models=False,MTA_retained=True,C4_role='SPP1-rich cholesterol-handling state localization',completion_time=datetime.datetime.now(datetime.timezone.utc).isoformat())
(O/'STATUS.json').write_text(json.dumps(status,indent=2))
note="""## 2026-09-15 Phase160 外部公共队列与客户实验背景调查完成

用户确认整体巨噬SPP1联系+C4胆固醇处理状态定位双层分工并授权继续；本轮另提供湿实验.png作参考，无原始数据。客户shTYMP图只作背景：MTA/Methionine/SAM/Thymine/Spermidine下降，Thymidine上升，不重算、不进入公共验证，不索要原始数据。图SHA和使用边界见本阶段CLIENT_FIGURE_REFERENCE.json，MTA继续保留，SLC16A3不替换客户假说。
完成RC18/RC20/CPTAC/CPTAC_val四外部组织RNA–实测代谢组队列，97/52/50/56=255映射患者。报告phase160_external_bridge_investigation_20260915_v1/EXTERNAL_BRIDGE_AND_CLIENT_CONTEXT_CN.md。741项检查PASS，无图、无客户数据进入公共统计。
优先补充候选SLC16A3：髓系RNA代理+缺氧调整后TYMP—SLC16A3四队列均正，随机合并rho0.329/P2.02e-5/FDR8.06e-5；条件SLC16A3—SPP1四队列均正，rho0.357/P1.09e-4/FDR2.18e-4；留一队列两端均正。原始RC20上游略负、完全调整CPTAC上游弱，保留。组织层并非恶性来源/C4特异验证。
实测乳酸n255：原始TYMP—乳酸rho0.090/P0.477，乳酸—SPP1rho0.102/P0.207，SLC16A3—乳酸rho0.068/P0.289；CPTAC名义正向不代表全部队列。总胆固醇n149两端未支持，27-HC未测。
MTA—组织SPP1rho0.139/P0.0293为同旧数据复算，历史FDR0.0975保留，不宣称新增独立复现；TYMP—MTA随机rho0.058/P0.585，旧固定rho0.0395/P0.538一致。调整结果完整保留为探索性。
下一步：优先做SLC16A3的细胞来源/恶性细胞内TYMP对应及巨噬SPP1状态核查，不预设乳酸唯一机制；C4定位与MTA相关背景保留，无确认因果中介。

"""
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 p=P/name;shutil.copy2(p,O/(name+'.before'));p.write_text(note+p.read_text(encoding='utf-8'),encoding='utf-8')
# Preserve manifest schema; repair only our three recent four-column stage rows with backups.
p=P/'RESULTS_MANIFEST.tsv';shutil.copy2(p,O/'RESULTS_MANIFEST.before.tsv');lines=p.read_text(encoding='utf-8').splitlines();new=[lines[0]]
for line in lines[1:]:
 cells=line.split('\t')
 if len(cells)==4 and cells[0].startswith(('phase158_','phase159_')):
  phase,path,validation,notes=cells
  cells=[phase,phase,'report',path,'2026-09-15',validation,'stage STATUS.json',notes]
 new.append('\t'.join(cells))
new.append('\t'.join(['phase160_external_bridge_investigation_20260915_v1','Phase160','report',str(O/'EXTERNAL_BRIDGE_AND_CLIENT_CONTEXT_CN.md'),status['completion_time'],'PASS','741 computational checks','255 patients; observed public metabolomics; client figure background only']))
p.write_text('\n'.join(new)+'\n',encoding='utf-8')
p=P/'part2_rebuild_v1/STATUS.json';shutil.copy2(p,O/'PART2_STATUS.before.json');b=json.loads(p.read_text());b.update(latest_phase='Phase160',latest_report=str(O/'EXTERNAL_BRIDGE_AND_CLIENT_CONTEXT_CN.md'),stage='external_bridge_investigation_complete',confirmed_mediator=False,C4_retained=True,MTA_retained=True,priority_candidate='SLC16A3');p.write_text(json.dumps(b,ensure_ascii=False,indent=2))
manifest={x:hashlib.sha256((O/x).read_bytes()).hexdigest() for x in ['SUMMARY_CN.md','EXTERNAL_BRIDGE_AND_CLIENT_CONTEXT_CN.md']}
(O/'DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(manifest));print('COMPLETE')

