from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,datetime,shutil
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=P/'phase163_polyamine_MTA_serial_receiver_20260915_v1'
M=pd.read_csv(O/'META_CONDITIONAL_EDGES.tsv',sep='\t');R=pd.read_csv(O/'COHORT_SERIAL_DECOMPOSITION.tsv',sep='\t');E=pd.read_csv(O/'COHORT_CONDITIONAL_EDGES.tsv',sep='\t');RP=pd.read_csv(O/'MTA_RECEIVER_PANEL_SUMMARY.tsv',sep='\t')
def md(df):
 out=['| '+' | '.join(df.columns)+' |','|'+'|'.join(['---']*len(df.columns))+'|']
 for row in df.itertuples(index=False,name=None):out.append('| '+' | '.join(f'{x:.6g}' if isinstance(x,(float,np.floating)) else str(x) for x in row)+' |')
 return '\n'.join(out)
def st(x,e,adj='MTAP_myeloid'):
 r=M[(M.candidate==x)&(M.edge==e)&(M.adjustment==adj)].iloc[0]
 return f"rho={r.random_rho:.3f}、P={r.random_P:.3g}、FDR={r.random_FDR:.3g}"
refs=[
dict(source='Hansen et al., Scientific Reports 2022',url='https://www.nature.com/articles/s41598-022-07697-0',context='GBM and mouse macrophage models',evidence='MTA can promote alternative macrophage polarization involving ADORA2B and STAT3; Arg1/IL10/VEGFA endpoints',scope='Mechanistic prior, not proof of human ccRCC SPP1 or frozen C4; search-accessible primary abstract/method excerpts; full-text access was intermittently blocked'),
dict(source='Keyel et al., 2014',url='https://pmc.ncbi.nlm.nih.gov/articles/PMC4130577/',context='TLR-stimulated macrophage models',evidence='MTA alters inflammatory activation through adenosine receptor signaling',scope='Context-dependent activation, not a specific C4 signature'),
dict(source='GSE249082',url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249082',context='MTA-treated fibroblasts and secreted purine derivatives',evidence='Public study describes fibroblast metabolic reprogramming and indirect macrophage polarization',scope='Already inventoried in Phase77; fibroblast RNA signature cannot be relabeled direct macrophage MTA response')]
(O/'MECHANISTIC_SOURCE_REVIEW.json').write_text(json.dumps(refs,ensure_ascii=False,indent=2))
retained="""# 已确认保留的MTA–多胺轴叙事

第一部分显示，TYMP阳性恶性细胞伴随甲硫氨酸/多胺相关代谢状态，SRM在14名配对患者中11名升高（P=0.00671、FDR=0.0403）。外部公共组织队列中，SMS与实测MTA及SPP1呈正向关联，MTAP与实测MTA总体呈负向关联。上述结果支持保留多胺合成与MTA清除共同参与的代谢背景，并与客户TYMP干预后的MTA变化形成分层对应。

当前可写为分段关联支持，不写成已确认的TYMP—SRM/SMS—MTA—SPP1串联中介。MTA作为实验相关机制候选继续保留，MTAP单列清除和遗传背景，SLC16A3作为并行补充线索。整体巨噬承载SPP1功能联系，C4保留SPP1信号集中和胆固醇处理状态定位。
"""
(O/'RETAINED_NARRATIVE_CN.md').write_text(retained,encoding='utf-8')
summary=f"""# 保留多胺–MTA部分，串联检验定位到下游缺口

Phase162的第一部分SRM上调、公共SMS—MTA/SPP1以及MTAP清除端结果继续保留，已写入项目规则和保留叙事。

本轮在同四队列255患者中，以SRM、SMS、ODC1/SRM/SMS生产端RNA模块分别作为X，检查TYMP—X—实测MTA—组织SPP1。完成36个队列模型，每个2000次患者bootstrap；所有36个串联乘积区间均跨零。

主模型控制MTAP和髓系RNA代理：
- SMS—MTA保留正向，{st('SMS','d')}，四队列均正。
- 进一步控制SMS后，MTA—SPP1为{st('SMS','b')}，接近零。
- 控制MTA后，SMS—SPP1仍为{st('SMS','X_direct')}，四队列均正。
- SRM的TYMP上游得到支持，但SRM—MTA条件相关接近零；生产端模块的完整串联同样没有确认。
这表明当前组织数据不支持用实测MTA解释完整串联，不能据此删除此前分段联系或排除所有MTA作用。

接收端补查24患者：ADORA2B/ADORA2A在C4和其他巨噬中均未达到任何患者10%细胞检出的覆盖标准；C4的STAT3、SOCS3、IL10、VEGFA转录量多低于其他巨噬。不能据RNA判断受体蛋白缺失或STAT3磷酸化，更不能直接套用其他癌种的MTA–ADORA2B–STAT3机制。

文献保留直接巨噬作用与成纤维细胞中转两类机制先验；GSE249082为成纤维细胞背景，不能当作肾癌巨噬直接MTA响应验证。

本轮600项独立/一致性检查通过，无图、无客户数据入模。主线保留，当前缺口集中在MTA到SPP1的具体作用方式与细胞位置。
"""
(O/'SUMMARY_CN.md').write_text(summary,encoding='utf-8')
report="""# 多胺—MTA—SPP1同队列串联检验与接收端核查

## 结论：分段证据保留，当前组织数据没有确认MTA串联中介

"""+summary.split('\n\n',1)[1]+"""

## 已确认保留的内容

"""+retained.split('\n\n',1)[1]+"""

## Methods

沿用Phase162四个公共配对参考队列的255名肿瘤患者，人数97/52/50/56。每个队列内分别以SRM、SMS、生产端RNA模块为X；模块由ODC1/SRM/SMS既定队列标准化表达构成，保持其定义。A为TYMP，Z为实测MTA，Y为组织SPP1。所有变量取标准化秩。

拟合X~A+C、Z~A+X+C、Y~A+X+Z+C，其中C分别为无额外协变量、MTAP、MTAP加髓系RNA代理。MTAP加髓系代理为预列主要模型，原始与仅MTAP模型为敏感性。计算a（A到X）、d（X到Z）、b（Z到Y），以及X与Y、A与Z、A与Y的剩余条件系数。串联乘积为a*d*b，使用总效应等于各路径之和的恒等式核对。总效应不要求预先显著，不计算中介百分比。

各队列每模型进行2000次患者有放回重抽样，每次重新取秩和标准化，计算串联乘积的百分位95%区间。生产端RNA分数为预定义观测特征，不在bootstrap中按结局重选基因。三候选、三调整、四队列共36模型。

各条件边采用偏秩相关，P按有效自由度计算；以Fisher z和协变量修正方差进行随机效应合并，报告I²。BH FDR在每个队列或合并层、每条边和调整方式下的三候选家族内计算。新增三候选FDR不代替Phase162六基因家族FDR。留一队列结果完整保存。不将分别合并后的a、d、b相乘来制造跨队列串联效应。

接收端依据公开机制先验，预列ADORA2B、ADORA2A、STAT3、SOCS3、IL10、VEGFA、ARG1。在冻结肿瘤C4与其他成熟巨噬中按患者聚合；两区室均至少20细胞，共24患者。报告检出率和患者内log2CPM差，七基因构成单独表达比较家族。10%检出阈值仅用于覆盖描述。

## Results

### SMS的MTA关联保留，但MTA未解释SMS—SPP1联系

控制TYMP、MTAP和髓系代理后，SMS—MTA的条件相关仍为正，四队列一致，合并rho=0.250、P=0.000578、本轮三候选FDR=0.00173。此系数复现Phase162的相应关联，其原六基因层面FDR=0.00289仍保留。

在SPP1模型中再控制SMS后，MTA—SPP1合并rho=0.0112、P=0.866、FDR=0.866，仅两个队列方向为正；而SMS与SPP1的剩余条件相关四队列均为正，rho=0.312，完整P/FDR见表。当前组织数据不支持把SMS—SPP1联系归结为已测组织MTA的串联解释。

原始模型敏感性同样显示：仅控制TYMP、不加SMS时，MTA—SPP1合并rho=0.132、P=0.0399；加入SMS后rho=0.0396、P=0.544，而SMS与SPP1的条件rho=0.353、P=2.15e-7。加入MTAP和髓系代理本身也使不含X的MTA—SPP1基线变弱（rho=0.0870、P=0.185），因此不能把全部变化只归因于SMS。

### SRM与生产端模块的完整串联没有确认

主模型TYMP—SRM合并rho=0.230、P=0.00213、FDR=0.00577；但SRM—MTA在控制TYMP、MTAP和髓系代理后rho=0.00292、P=0.965。SRM对应第一部分上游状态的事实保留，而其作为单个串联X的中间一段仍弱。

生产端模块的TYMP上游rho=0.228、P=0.00384、FDR=0.00577；模块—MTA条件rho=0.139、P=0.0865、FDR=0.130；控制模块后的MTA—SPP1rho=0.0573、P=0.388、FDR=0.582。模块与SPP1的剩余条件相关仍为正，rho=0.257、P=7.73e-5、FDR=0.000116。

主模型12个队列/候选组合中，只有5个组合的a、d、b同时为正；三个调整口径全部36个串联乘积的bootstrap区间均跨零。没有将某个队列或单个乘积正向当作整条链的确认。

### 经典腺苷受体机制没有在当前RNA中获得接收端支持

24名配对患者中，C4的ADORA2B检出率中位0，其他成熟巨噬中位0.239%；两区室没有患者达到10%检出标准。ADORA2A中位均0，也没有达到该覆盖标准。ARG1在本面板内未检出。

STAT3在两个区室有表达，但C4相对其他巨噬23/24患者更低；SOCS3、VEGFA同为23/24更低，IL10为21/24更低。此结果不支持直接将这组转录标志物命名为当前C4的MTA响应程序。

## 外部机制资料核查

GBM与鼠巨噬模型的原始研究报道MTA可经ADORA2B/STAT3相关机制促进替代性巨噬激活，考察Arg1、IL10、Vegfa等终点。它提供机制先验，不能直接作为人肾癌SPP1/C4已获验证的依据。[Scientific Reports原始研究](https://www.nature.com/articles/s41598-022-07697-0)

另一原始研究报道MTA通过腺苷受体改变TLR刺激下的巨噬炎症反应，提示作用依赖细胞及刺激背景。[MTA与巨噬激活原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC4130577/)

GSE249082公开研究描述MTA被成纤维细胞摄取代谢，并通过其分泌的嘌呤衍生物影响巨噬极化。这是可继续核查的间质中转路线。该数据此前已在Phase77登记，其成纤维细胞RNA签名不能重命名为直接巨噬MTA响应。[GSE249082](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249082)

本阶段没有取得可直接重算“人肾癌MTA处理巨噬—SPP1”的独立转录组输入，以上机制文献不计作本项目患者验证。部分全文访问受站点访问限制，来源核查范围记录于MECHANISTIC_SOURCE_REVIEW.json；未以受阻页面的搜索失败推断某项实验不存在。

## Discussion与后续定位

本轮结果保留第一部分SRM上调、公共SMS—MTA/SPP1、MTAP—MTA等分段事实，同时限定其串联解释。条件相关减弱不是因果排除：组织代谢物水平可能无法代表局部细胞外暴露，bulk表达也混合多个区室；横断面协变量可能位于共同状态或作用路径上。当前能说的是同批组织数据没有确认特定串联模型，不能说MTA不参与任何巨噬作用。

每个队列没有时间先后与干预，a*d*b只是标准化秩的关联分解。bootstrap没有消除未测混杂、测量误差或方向不确定性；四队列合并近似与中等异质性保留。没有用FDR阴性作为停止或删除线索的门槛。

ADORA2A/B低RNA检出不能证明受体蛋白不存在；STAT3转录量也不代表磷酸化活性。不同癌种、鼠巨噬与人肾癌之间的差异不能通过一个小标志物列表消除。当前接收端证据不足以给C4指定经典MTA—ADORA2B—STAT3机制。

因此，MTA部分继续保留，正文定位为代谢背景和分段联系；下一步重点应是下游具体作用方式与细胞位置，可核查直接巨噬作用和间质中转两类路线。MTAP仍作清除背景，SLC16A3并行，整体巨噬与C4双层分工不变。客户图仅为实验背景，不进入公共模型或替代缺失的中介证据。

## 审计与完整统计

600/600检查通过，包括路径恒等式、独立精度矩阵偏相关、bootstrap抽样重建、合并效应和FDR复核。所有36模型、72000个串联bootstrap数值、单队列条件边、留一队列和接收端患者明细均服务器保存。无新图、无客户数据入模。

"""
for title,df in [('所有条件边合并结果',M),('主模型串联分解',R[R.adjustment=='MTAP_myeloid']),('全部串联区间',R[['cohort','candidate','adjustment','n','serial_adb','serial_low','serial_high','bootstrap_positive_fraction']]),('接收端表达',RP)]:
 report+='\n\n### '+title+'\n\n'+md(df)
(O/'SERIAL_MTA_AND_RECEIVER_REVIEW_CN.md').write_text(report,encoding='utf-8')
status=json.loads((O/'STATUS.json').read_text());assert status['total_checks']==600
status.update(status='COMPLETE',report=str(O/'SERIAL_MTA_AND_RECEIVER_REVIEW_CN.md'),phase162_retained=True,MTA_retained=True,confirmed_serial_mediator=False,completion_time=datetime.datetime.now(datetime.timezone.utc).isoformat())
(O/'STATUS.json').write_text(json.dumps(status,indent=2))
note="""## 2026-09-15 Phase163 保留多胺轴并完成同队列串联/接收端核查

用户确认Phase162部分保留并继续。已写入AGENTS及RETAINED_NARRATIVE_CN.md，原第一部分SRM/公共SMS—MTA/SPP1等分段事实不删除。报告phase163_polyamine_MTA_serial_receiver_20260915_v1/SERIAL_MTA_AND_RECEIVER_REVIEW_CN.md。
同4队列255患者，SRM/SMS/生产端模块三候选，raw/MTAP/MTAP+髓系三口径，共36模型各2000患者bootstrap。全部36个a*d*b区间跨零，主模型12组合仅5个a/d/b同时正。无确认串联中介。
主模型SMS—MTA条件rho0.250/P0.000578/本轮FDR0.00173，原Phase162六基因FDR0.00289保留；加入SMS后MTA—SPP1rho0.011/P0.866，而SMS—SPP1控制MTA仍rho0.312正。SRM上游支持，但SRM—MTA条件rho0.0029/P0.965。MTA背景不含X但已调TYMP/MTAP/髓系时rho0.087/P0.185，不把全部衰减归因SMS。
24患者接收端：C4/其他巨噬ADORA2A/B均无患者达到10%RNA检出，C4 STAT3/SOCS3/VEGFA多数降低；不能判定蛋白缺失或pSTAT3活性，也不能直接套用GBM鼠巨噬ADORA2B机制。
机制先验保留直接巨噬作用（SciRep2022 s41598-022-07697-0）及成纤维细胞中转（既有GSE249082），未得到直接人肾癌MTA处理巨噬SPP1独立输入。600检查PASS，无图无客户数据入模。
推荐保留分段证据，继续定位MTA下游作用方式与细胞位置；MTAP清除背景、SLC16A3并行、整体巨噬+C4双层分工不变。

"""
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 p=P/name;shutil.copy2(p,O/(name+'.before'));p.write_text(note+p.read_text(encoding='utf-8'),encoding='utf-8')
with (P/'RESULTS_MANIFEST.tsv').open('a',encoding='utf-8') as f:f.write('\t'.join(['phase163_polyamine_MTA_serial_receiver_20260915_v1','Phase163','report',str(O/'SERIAL_MTA_AND_RECEIVER_REVIEW_CN.md'),status['completion_time'],'PASS','600 computational checks','255 patients; 36 models; 72000 bootstrap draws; Phase162 retained; no figures'])+'\n')
p=P/'part2_rebuild_v1/STATUS.json';shutil.copy2(p,O/'PART2_STATUS.before.json');b=json.loads(p.read_text());b.update(latest_phase='Phase163',latest_report=str(O/'SERIAL_MTA_AND_RECEIVER_REVIEW_CN.md'),stage='polyamine_MTA_serial_receiver_complete',MTA_retained=True,C4_retained=True,confirmed_mediator=False);p.write_text(json.dumps(b,ensure_ascii=False,indent=2))
manifest={x:hashlib.sha256((O/x).read_bytes()).hexdigest() for x in ['SUMMARY_CN.md','SERIAL_MTA_AND_RECEIVER_REVIEW_CN.md','RETAINED_NARRATIVE_CN.md']};(O/'DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(manifest));print('COMPLETE')

