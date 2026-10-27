from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,shutil,datetime
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=P/'phase159_metabolic_parts_C4_value_20260915_v1'
R=pd.read_csv(O/'ALL_COMPONENT_BRIDGE_RESULTS.tsv',sep='\t');PR=pd.read_csv(O/'C4_OTHER_GENE_SUMMARY.tsv',sep='\t');A=pd.read_csv(O/'FAIR_ENDPOINT_COMPARISON.tsv',sep='\t');DI=pd.read_csv(O/'PAIRED_CORRELATION_DIFFERENCE_BOOTSTRAP.tsv',sep='\t');D=pd.read_csv(O/'SAME_LIBRARY_C4_CONTRIBUTION.tsv',sep='\t');RC=pd.read_csv(O/'MACROPHAGE_SIDE_TYMP_X_SPP1.tsv',sep='\t');ST=pd.read_csv(O/'ALL_COMPONENT_LEAVE_ONE.tsv',sep='\t')
def md(df):
 out=['| '+' | '.join(df.columns)+' |','|'+'|'.join(['---']*len(df.columns))+'|']
 for row in df.itertuples(index=False,name=None):out.append('| '+' | '.join(f'{x:.6g}' if isinstance(x,(float,np.floating)) else str(x) for x in row)+' |')
 return '\n'.join(out)
summary="""# C4具有明确贡献，代谢拆解仍未确定完整中介

已完成胆固醇稳态、糖酵解拆解：38项发送端候选（31基因、7环节模块），并检查巨噬侧相关分子。

**C4的正向价值已明确。** 同16患者中，各患者C4细胞占比中位数7.4%，C4贡献的成熟巨噬SPP1原始RNA计数占比中位数54.6%；直接回查冻结H5AD全部16患者一致。该值是所捕获RNA，不是分泌蛋白或组织真实总量。此前约13%的对数表达贡献指标不是原始RNA计数占比；此外旧分析患者集也不同，不能混用。

24名可做肿瘤内C4/其他成熟巨噬配对的患者中，CYP27A1为24/24升高（中位log2CPM差+4.134，P=1.19e-7，FDR=1.79e-6）；NR1H3为23/24升高（+1.605，P=2.38e-7，FDR=2.38e-6）；SPP1为22/24升高（+1.447，P=2.47e-5，FDR=2.47e-5）。C4值得保留为SPP1信号集中的胆固醇处理状态。

**相同口径下C4没有显示普遍更差。** 同16患者，TYMP—C4 rho=0.468、P=0.0678；TYMP—整体巨噬SPP1表达rho=0.106、P=0.696。相关系数差的患者bootstrap区间跨零，未证明C4显著更优；此前17人SPP1名义显著不能直接与此比较。

**具体环节有线索，但两端未完全接上。**
- 肿瘤MSMO1与TYMP rho=0.759，P=0.000654，FDR=0.0118；与C4边际rho=0.591，P=0.0159，FDR=0.143；与整体巨噬SPP1表达rho=0.065，P=0.812。适合保留在TYMP相关甾醇状态—C4分支。
- 肿瘤HK2与TYMP名义正相关，但下游SPP1对应弱。
- 肿瘤LDHA、SLC16A3与巨噬SPP1表达均名义正相关，rho约0.59，P约0.015；各自控制TYMP后的FDR约0.063。但TYMP—LDHA rho=0.053、P=0.846，TYMP—SLC16A3 rho=0.300、P=0.259，上游没有接牢。
- 巨噬CYP27A1与肿瘤TYMP rho=0.582、P=0.0179、FDR=0.143；与SPP1为正向趋势。控制C4比例后TYMP关联减弱，部分联系与组成变化相容。
- 巨噬LDHA/SLC16A3与SPP1的关联在控制C4比例后仍存在，说明有广泛巨噬状态联系，不限于C4。

当前应保留“整体巨噬SPP1功能联系”和“C4的细胞状态定位”两层证据，不能把C4设为TYMP作用的唯一必经终点，也不应删去已经确立的C4正向结果。MTA继续保留；无新图、无客户湿实验使用。431项核查通过。
"""
(O/'SUMMARY_CN.md').write_text(summary,encoding='utf-8')
report="""# 代谢环节拆解与C4实际贡献复核

## 结论：C4具有SPP1与胆固醇处理特征，尚未确立TYMP特异作用路径

"""+summary.split('\n\n',1)[1]+"""

## Methods

沿用Phase154的16名同肿瘤文库配对患者。发送端为恶性细胞伪bulk log2CPM，终点分别为C4占成熟巨噬比例、成熟巨噬SPP1伪bulk表达、SPP1检出比例。按公开路径知识预列38个发送端基因，至少12/16文库有10个原始计数者进入关联；保留31个基因。九类环节中，至少2个且至少一半预列基因达到覆盖的环节形成均值z分数，共7个模块。TYMP、SPP1不参与模块。候选和来源见PLAN_BEFORE_FIT.json，覆盖排除表完整保留。

胆固醇分为甲羟戊酸、甾醇合成、氧固醇生成、后续转化清除、外排/载体；糖酵解分为葡萄糖进入、核心反应、乳酸相互转化及单羧酸转运。LDHA与LDHB联合分数仅称乳酸转化模块，分别单基因结果同时保留。SLC16A1、SLC16A3及BSG分数表示转运相关RNA，不预设净流向。

统计使用标准化秩；a为TYMP—X，XY为X—终点，b为控制TYMP的X—终点偏相关。分别报告原始和研究指示变量调整，偏相关P以有效残差自由度的t近似计算。基因与模块、胆固醇与糖酵解、终点、边及调整口径分别作BH FDR。所有候选做16次留一患者与4次留一研究，并统一补充缺氧调整。FDR用于显著性表述，不作为是否保留方向性候选的硬门槛。

巨噬侧检查预列的脂质摄取、LXR相关、氧固醇代谢和单羧酸转运/感应基因，分别在整体成熟巨噬、C4、其他成熟巨噬计算TYMP—X—SPP1关联。区室至少20细胞，整体及其他成熟巨噬为16患者，C4为12患者；基因至少75%患者有10计数。其FDR按区室、路径、边、调整分别计算。另对CYP27A1、LDHA、SLC16A3追加C4比例协变量，作为组成敏感性分析。

C4状态表型采用全部肿瘤样本患者内聚合，C4与其他成熟巨噬均至少20细胞，共24患者。以配对Wilcoxon检验log2CPM差；标志物、胆固醇基因、糖酵解基因分别校正。所有计数和检出率以原始counts计算。

公平比较在同16患者以及C4至少20细胞的共同12患者分别进行，比较C4比例、整体SPP1表达/检出率及区室SPP1表达。相关差异采用相同患者、研究分层的2000次成对bootstrap，每次重新计算秩。新增终点审阅FDR单独标注，不替换Phase157/158的发现家族FDR。C4原始SPP1计数贡献逐患者计算，并独立从冻结H5AD读取SPP1列核对。

## Results

### C4承载了集中的SPP1转录信号

16患者的C4比例中位数7.406%，C4占成熟巨噬捕获SPP1计数的比例中位数54.590%，占SPP1阳性巨噬细胞数的比例中位数13.083%。RNA计数贡献范围0–93.455%，患者差异明显。三个中位数分别按患者计算，不是将所有细胞混池后的数值。

直接从原始H5AD的CSR counts定位SPP1计数，并以固定cell_id、患者、文库和C4成员回查，16/16患者的C4计数贡献、成熟巨噬SPP1总计数、细胞数全部一致。该结果说明少量SPP1表达较高的C4可贡献较多捕获RNA；“SPP1阳性细胞中仅约13%是C4”不等价于“C4仅贡献约13%的SPP1分子”。

此前Phase153约12.93%是另一患者集的对数表达指数贡献，不能当作原始SPP1计数贡献。两种算法回答不同问题，旧指标仍保留；需要纠正的是用旧指标推断C4总RNA贡献较小的解读。

### C4的胆固醇处理表型获得患者内一致支持

24患者中，CYP27A1全部在C4升高，中位log2CPM差+4.134；NR1H3为23/24升高，+1.605；ABCA1为19/24升高，+0.253；SCARB1为22/24升高，+0.820，均通过各自预列基因家族FDR。SPP1为22/24升高，+1.447。GPNMB与CCL18均24/24升高。

这组结果支持将C4定位为具有SPP1表达及胆固醇处理特征的巨噬状态。CYP27A1在巨噬侧高度富集，提示其代谢所在区室需要重视，不能只假定所有候选代谢物都由肿瘤产生后被C4接收。

与此并存的反向结果包括CH25H在21/24降低、GPR132在23/24降低；HCAR1的C4检出率中位数为0，24患者均未达到10%覆盖标准。SLC16A3及BSG在C4较高，而SLC16A1没有明确升高。当前数据不支持以“C4富集乳酸感应受体”作为首选解释。

### 相同患者下没有证据证明整体SPP1终点普遍优于C4

同16患者，TYMP、IGFBP3、胆固醇稳态、糖酵解与C4比例的原始rho分别为0.468、0.715、0.585、0.565，均高于对应整体SPP1表达的0.106、0.062、0.412、0.385。TYMP与C4和整体SPP1相关之差为0.362，bootstrap 95%区间[-0.357,1.017]；因此不能据此宣布C4统计上更优，也不能用此前异口径的17人SPP1结果宣布C4更差。

12名共同患者中，TYMP—C4比例降至rho=0.105、P=0.746；这是删去C4较少患者后的敏感性结果，不能取代保留低C4值的16人主比例分析。C4内部SPP1与三个上一轮肿瘤候选仍无明确对应。亚群平均表达、亚群比例与整体表达分别反映不同层面的变化。

### 甾醇合成中的MSMO1连接TYMP与C4，尚未连接整体SPP1表达

肿瘤MSMO1—TYMP rho=0.759、P=0.000654、FDR=0.0118；MSMO1—C4边际rho=0.591、P=0.0159、FDR=0.1429。控制TYMP后的C4关联rho=0.410、P=0.1286、FDR=0.5584；研究调整后进一步减弱。MSMO1—整体SPP1表达rho=0.065、P=0.812。缺氧调整后TYMP关联仍为rho=0.700、P=0.00368、FDR=0.0663，说明上游对应没有被缺氧分数完全解释。

FDPS与C4的条件关联较强（rho=0.750、P=0.00128、FDR=0.0230），但TYMP—FDPS rho=0.144、P=0.594。IDI1与TYMP正相关，控制TYMP后下游关联转负。胆固醇外排/载体模块的下游方向为负。不能把MSMO1的上游与FDPS的下游合并为同一个分子中介。

肿瘤CYP27A1—TYMP为负向（rho=-0.221、P=0.412），与C4同样为负向。DHCR7仅10/16文库达到预定计数阈值，未进入关联；其排除表示当前覆盖不足，不是阴性检验。CH25H仅1/16达到阈值。当前公共数据不能复刻一条已确定的“肿瘤TYMP—CYP27A1/27-HC—C4”路径。

### 乳酸转化和转运环节对应整体SPP1，上游TYMP仍缺少一致连接

肿瘤HK2—TYMP rho=0.659、P=0.00551、FDR=0.0716，为名义正向；HK2—整体SPP1表达rho=0.200、P=0.458，下游较弱。

LDHA/LDHB转化模块与整体SPP1表达的边际rho=0.647、P=0.00674、FDR=0.0270；控制TYMP后rho=0.655、P=0.00801、FDR=0.0320。研究调整后条件rho=0.666、P=0.0181、FDR=0.0723。但TYMP—该模块rho=-0.038、P=0.888；这条下游关联成立，并未同时构成TYMP的正向中介。

单基因LDHA与整体巨噬SPP1表达rho=0.597、P=0.0146、FDR=0.0656；肿瘤SLC16A3对应rho=0.594、P=0.0152、FDR=0.0656。控制TYMP后两者均为名义正向、FDR约0.0634。缺氧调整后下游效应仍约0.59。上游TYMP—LDHA rho=0.053、P=0.846；TYMP—SLC16A3 rho=0.300、P=0.259。当前可保留LDHA/SLC16A3作为SPP1相关代谢环节，不能声称TYMP已驱动乳酸升高。

### 巨噬CYP27A1提供另一条位置明确的候选，广泛巨噬糖酵解联系不依赖C4比例

整体巨噬CYP27A1—肿瘤TYMP rho=0.582、P=0.0179、FDR=0.1435；控制TYMP的CYP27A1—SPP1 rho=0.429、P=0.110、FDR=0.2937。控制C4比例后上游rho降为0.393、P=0.147，下游rho=0.335、P=0.241。研究与C4调整结果完整保留。这一方向性线索与细胞组成贡献相容，尚未显示独立于组成的中介联系。

整体巨噬LDHA和SLC16A3与SPP1分别rho=0.788和0.815；控制TYMP、C4比例后仍约0.818和0.794，敏感性家族FDR均约0.00103；同时研究调整也保持正向。但它们与恶性TYMP的上游相关接近零。C4内LDHA—SPP1在12患者同样有正向关联，而其他巨噬也有对应。糖酵解相关SPP1状态具有广泛巨噬分布，不能因C4发现而忽略整体巨噬的证据。

## Discussion与报告定位

目前应使用两层叙事：整体巨噬层呈现SPP1相关代谢状态及跨区室候选联系；C4层提供SPP1信号集中、CYP27A1/LXR相关胆固醇处理表型和细胞定位。C4值得保留，但“TYMP主要/特异通过C4作用”尚不是已确立结果。并列证据无需强行按一种亚群终点统一。

C4由同一表达数据发现，标志物富集是冻结亚群后的内部表型核对，不是外部独立验证。RNA计数占比受捕获效率、细胞总RNA及状态影响，不是组织真实分子总量，也不是分泌蛋白占比；患者内标准化表达差提供了另一层支持。个别研究样本多，研究调整和留研究结果需要同时保留；n<10留研究结果只作探索。后续筛选、模块拆解、接收端分析均属同一发现数据的探索，不能把新增家族FDR当作独立复现。

LDHA/LDHB和单羧酸转运RNA不能确定乳酸净生成或流向；CYP27A1、CH25H作用于不同氧固醇分支。当前没有同患者代谢物实测，不能把程序提升写成乳酸/27-HC提升。MTA未被本次分析排除。现有横断面数据不确定先后或因果；也不应因为未得到完整中介而否定已确立的C4表型和SPP1联系。

进一步优先核查的具体环节可收敛为：肿瘤MSMO1与C4表型的对应；肿瘤SLC16A3及巨噬LDHA/SLC16A3与SPP1的对应；巨噬CYP27A1与TYMP的联系。三者各有缺口，保留完整负向/不支持结果，不再仅凭一个整体Hallmark名称推定特定代谢机制。

## 公共机制来源

胆固醇合成基因及环节参考[Reactome胆固醇合成路径](https://reactome.org/content/detail/R-HSA-191273)。CYP27A1所在代谢路径参考[Reactome氧固醇路径](https://reactome.org/content/detail/R-HSA-193807)。两者支持环节划分，不意味着本项目已测得代谢流。

SLC16A3为单羧酸转运蛋白，其功能注释支持考察转运环节，RNA表达本身不确定净输出。[UniProt SLC16A3](https://www.uniprot.org/uniprotkb/O15427/entry)

GPR132的乳酸相关巨噬作用来自乳腺癌实验体系。本项目C4中GPR132反而较低，故保留这一接收端矛盾，不将该机制直接移植为C4解释。[GPR132原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC5255630/)

## 完整统计与审计附表

"""
report+='\n\n### C4表型基因全表\n\n'+md(PR)
report+='\n\n### 同患者原始RNA贡献及区室SPP1\n\n'+md(D)
report+='\n\n### 同患者终点比较\n\n'+md(A)
report+='\n\n### 成对相关差异bootstrap\n\n'+md(DI)
key=['gene:MSMO1','gene:FDPS','gene:HK2','gene:LDHA','gene:SLC16A3','gene:CYP27A1','module:glycolysis_lactate_conversion']
cols=['candidate','endpoint','adjustment','n_patients','a_rho','a_P','a_FDR','XY_rho','XY_P','XY_FDR','b_partial_rho','b_P','b_FDR']
report+='\n\n### 重点发送端候选全部终点\n\n'+md(R[R.candidate.isin(key)][cols])
report+='\n\n### 巨噬侧重点候选\n\n'+md(RC[RC.gene.isin(['CYP27A1','NR1H3','LDHA','SLC16A3'])])
sl=[]
for (candidate,endpoint,kind),g in ST[ST.candidate.isin(key)].groupby(['candidate','endpoint','kind']):
 sl.append(dict(candidate=candidate,endpoint=endpoint,kind=kind,n_fits=len(g),a_positive=int((g.a_rho>0).sum()),b_positive=int((g.b_partial_rho>0).sum()),ab_positive=int((g.indirect_ab>0).sum()),ab_min=g.indirect_ab.min(),ab_max=g.indirect_ab.max()))
report+='\n\n### 留一结果摘要\n\n'+md(pd.DataFrame(sl))
report+="""

全部38项发送端、其缺氧敏感性、接收端、覆盖排除与逐患者留一结果均以TSV保留在本目录。关联结果表列出未挑选的所有效应及P/FDR。仅作同数据计算复核，不能称外部验证。

431项核查通过：发送端与比例检查344项，原始H5AD计数/细胞核查48项，巨噬侧直接Spearman复核39项。脚本与输入SHA256、独立原始计数核查表、内存门槛日志全部服务器归档。无新图、无客户湿实验使用。
"""
(O/'METABOLIC_COMPONENTS_AND_C4_VALUE_CN.md').write_text(report,encoding='utf-8')
status=json.loads((O/'STATUS.json').read_text());assert status['total_checks']==431
status.update(status='COMPLETE',report=str(O/'METABOLIC_COMPONENTS_AND_C4_VALUE_CN.md'),confirmed_mediator=False,C4_retained=True,completion_time=datetime.datetime.now(datetime.timezone.utc).isoformat())
(O/'STATUS.json').write_text(json.dumps(status,indent=2))
note="""## 2026-09-15 Phase159 代谢拆解与C4实际贡献复核完成

用户授权继续拆解，同时质疑C4价值。完成38发送端候选（31基因+7模块）、巨噬侧对应及同患者公平比较；431项核查PASS，无图无客户数据。报告phase159_metabolic_parts_C4_value_20260915_v1/METABOLIC_COMPONENTS_AND_C4_VALUE_CN.md。
关键更正：同16患者C4比例中位7.406%，其所捕获成熟巨噬SPP1原始RNA计数占比中位54.590%，直接回查原H5AD16患者全部一致。此前约13%为不同表达/阳性细胞指标和患者集，不能当原始RNA量占比。24配对患者CYP27A1在C4全部升高、NR1H3在23人、SPP1在22人升高，FDR通过。C4保留为SPP1信号集中的胆固醇处理状态；不能写TYMP特异作用已由C4解释。
同16患者TYMP-C4 rho0.468/P0.0678高于TYMP-整体SPP1 rho0.106/P0.696，但相关差区间跨零，不宣布显著优劣。原17患者不同汇总口径不能混比。
拆解：肿瘤MSMO1-TYMP rho0.759/P0.000654/FDR0.0118，MSMO1-C4名义正向而SPP1无对应；HK2贴TYMP，LDHA/SLC16A3贴整体SPP1但上游弱。巨噬CYP27A1-TYMP rho0.582/P0.0179/FDR0.143，控制C4比例减弱。巨噬LDHA/SLC16A3-SPP1关系控制C4仍在，属广泛巨噬联系。暂无确认中介，MTA保留。
叙事采用整体巨噬功能联系+C4细胞状态定位双层，不强制所有SPP1关系经过C4。后续具体线索MSMO1、SLC16A3/LDHA、巨噬CYP27A1，各保留对应缺口。

"""
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 p=P/name;shutil.copy2(p,O/(name+'.before'));p.write_text(note+p.read_text(encoding='utf-8'),encoding='utf-8')
with (P/'RESULTS_MANIFEST.tsv').open('a',encoding='utf-8') as f:f.write('phase159_metabolic_parts_C4_value_20260915_v1\t'+str(O/'METABOLIC_COMPONENTS_AND_C4_VALUE_CN.md')+'\tCOMPLETE\t431 checks; frozen C4; 38 sender candidates; original RNA verified\n')
p=P/'part2_rebuild_v1/STATUS.json';shutil.copy2(p,O/'PART2_STATUS.before.json');b=json.loads(p.read_text());b.update(latest_phase='Phase159',latest_report=str(O/'METABOLIC_COMPONENTS_AND_C4_VALUE_CN.md'),stage='metabolic_components_and_C4_value_complete',confirmed_mediator=False,C4_retained=True);p.write_text(json.dumps(b,ensure_ascii=False,indent=2))
manifest={x:hashlib.sha256((O/x).read_bytes()).hexdigest() for x in ['SUMMARY_CN.md','METABOLIC_COMPONENTS_AND_C4_VALUE_CN.md']}
(O/'DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(manifest));print('COMPLETE')

