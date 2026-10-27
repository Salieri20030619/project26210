from pathlib import Path
import json,hashlib,datetime
import pandas as pd
w=Path(__file__).resolve().parents[1];q=w/'qc'
b=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/09_metabolism_trajectory')
reg=json.loads((q/'FIGURE_REGISTRY.json').read_text());validation=json.loads((q/'INDEPENDENT_VALIDATION.json').read_text())
assert validation['status']=='PASS' and 'PASS' in (q/'R_NEW_STATISTICS_VALIDATION.txt').read_text()
read=lambda p:pd.read_csv(p,sep='\t')
def fulltable(p):return read(p).to_csv(sep='\t',index=False,float_format='%.10g')
notes={}
notes['2_6']='''2.6 恶性上皮代谢通路与酶/转运体图件 v1

本包共9套图，每套提供300 dpi PNG及PDF，加本说明共19个文件。建议先看01–03；04–05为指定酶和转运体；06–08为scFEA预测代理；09为全部85条通路的探索性补充。

结果：三条代谢通路与TYMP/F1呈一致正向关联
14位合格配对患者中，TYMP阳性相对阴性的嘧啶、半胱氨酸/蛋氨酸、精氨酸/脯氨酸AUCell中位差分别为0.009482、0.010186、0.014007；正向患者分别13/14、14/14、13/14。P分别0.001696、0.001097、0.001367，四通路家族FDR均为0.002262。
F1高主克隆与其他主克隆的患者配对比较同样支持这三条通路升高。当前数据不支持甘油磷脂AUCell的配对差异：TYMP比较中位差=-0.0001500，7/14正向，P=FDR=0.801732；F1克隆比较中位差=0.0010186，9/14正向，P=FDR=0.379471。
20位合格患者的F1连续关联中，三条通路中位rho分别0.104698、0.167113、0.199657，FDR分别0.019298、0.002482、0.002482；甘油磷脂rho=0.050023，FDR=0.489783。
四个scFEA轴与F1均为正向患者关联并通过各自四轴家族FDR；其数值是原模块标准化后的预测代理，详见完整统计。
12个原方案指定基因的F1相关检验中，TYMP、TK1、ODC1、SRM、PCYT2通过12基因家族FDR。其余指定基因的完整结果全部保留；当前数据不支持MTAP随F1稳定降低，也不支持全部12个酶/转运体共同激活。
全85通路的TYMP配对比较为新增探索性补充，49/85通过该85项家族的FDR<0.05，包含正向和负向变化。不得把85项FDR与四条预设通路的FDR互换。

方法
复用正式11,950个恶性上皮细胞、26位患者、5队列；UMAP、cNMF、CNV克隆及原代谢模型保持冻结。AUCell使用原85个KEGG基因集及原始UMI计分，本次不重新运行AUCell、scMetabolism、Compass或scFEA。
TYMP阳性定义为原始UMI>0，阴性为UMI=0。配对比较要求同一患者两组各至少20细胞，以各组中位数之差作为患者效应。F1高克隆沿用原正式主克隆定义与0.25阈值，不改为F1主导模块或全局中位数分组。
连续关联在每位患者内部计算Spearman rho，要求至少30细胞且两变量有变异，再用患者rho进行双侧Wilcoxon signed-rank检验。原代谢检验使用正态近似及连续性校正；BH分别在每个算法、每个比较或每个预测变量的4通路/轴家族内计算。
新增85通路比较使用相同14位合格配对患者、相同组内中位数和双侧近似Wilcoxon，85项独立BH。新增12基因F1关联以原方案全部12基因为固定家族独立BH。基因的可估计患者数因表达变异而不同，完整n见图05和下表。
12基因直接从正式H5AD的counts层按精确细胞键提取，按全部保留基因的文库总UMI计算log1p(CP10K)，TYMP原始计数逐值一致，标准化表达在数值容差内一致。
图04为患者内F1秩四分位，保留同分值；只纳入四组各至少5细胞的17位患者。点面积为先患者内计算再跨患者等权平均的检出比例，颜色为同口径平均log1p表达，不按基因重新缩放。
scFEA轴沿用原模块集合，每模块在全细胞中作样本标准差z标准化再取轴内均值；精确恒定模块按原R safe_z置零。本次独立复算发现NumPy对恒定值的浮点标准差残差，已在复核脚本中按常量范围明确识别；原结果未改动，四轴最大误差<5e-14。胆碱代理M83为恒定模块，变化来自M21。

讨论与解读范围
AUCell反映转录组通路活性，scFEA为转录组推断的相对通量代理；两者不等于实际代谢物浓度或实测通量。MTA不是本模型显式输出，胆碱M21/M83只覆盖局部反应，不能写成完整肯尼迪通路通量。本输入没有匹配非靶代谢组，不能验证原方案所述44倍MTA变化。
TYMP阴性指单细胞UMI未检出，不能等同于蛋白缺失或基因敲除。酶RNA不等于酶活性，MTAP表达图不能证明内源性缺失；F1与代谢的相关不建立TYMP驱动的因果机制。全85通路比较为探索性扩展，队列结构、转录本检测和患者间异质性应结合解释。

图件索引
'''
notes['2_7']='''2.7 恶性上皮拟时序与代谢/RNA动态图件 v1

本包共10套图，每套提供300 dpi PNG及PDF，加本说明共21个文件。建议先看01、04、09；02–03展示CellRank命运概率与转移方向；05为scFEA代理动态；06–08为12个指定酶/转运体RNA动态；10为状态构成。

结果：代谢轨迹存在方法差异，当前数据不支持统一的恶性进展时间轴
Monocle3、Slingshot、scTour均提供正式的原始细胞排序。四通路中，Slingshot的嘧啶通路患者相关为正（中位rho=0.116058，14/20正向，P=0.00365448，原12终点家族FDR=0.04385376）。该方向未在Monocle3或scTour获得校正后支持；其余代谢通路的方向和完整统计均保留。
Slingshot和scTour的F1中位相关分别0.287632和0.288109，但FDR分别0.116258、0.113228；不能称为校正后显著。Monocle3的F1中位rho=-0.034036，FDR=0.313319。
CytoTRACE2干性/分化潜能与F1的患者关联为负（中位rho=-0.203566，5/20正向，P=0.0120799，FDR=0.048317），与甘油磷脂AUCell也为负（rho=-0.131633，3/20正向，P=0.00198555，FDR=0.0119133）。它衡量潜能，不作为第四种拟时序。
CellRank给出Cycle、Hypoxia、TYMP-state三类命运概率及原转移向量，所有细胞三项概率之和为1。12基因RNA曲线与四通路和scFEA代理均按三种主拟时序完整展示，没有为了符合预期而反转方向或筛掉方法。

方法
复用11,950细胞、26位患者、5队列的冻结UMAP、表达状态和既有模型。原Monocle3主图含1个partition和6个cluster；原Slingshot有3条lineage，正式主展示沿用已选lineage1。scTour沿用原200 epoch模型。本轮不重跑轨迹模型或改动root。
原root/方向锚点为患者内标准化的CNV burden、Hypoxia、Cycle和负向Base四项均值，未使用TYMP、F1或代谢分数；候选root细胞共922个。CytoTRACE2按原5队列分别运行，显示原始Score。CellRank沿用0.8 pseudotime kernel + 0.2 connectivity kernel及原三终态。
01三种主拟时序色阶统一0–1，为原患者内缩放数值；Monocle3叠加原真实principal graph。CytoTRACE2与CellRank分别标注潜能与最大命运概率。02三类命运概率共用0–1色阶。03将原CellRank转移向量按28×22个UMAP空间格汇总，至少10细胞的格显示平均向量，箭头共用显示尺度。
04–08和10采用20位至少30细胞患者。每位患者先按原拟时序秩分20箱，箱内求均值；固定Gaussian核带宽0.075秩单位，按箱细胞数加权拟合。04–08线为患者拟合曲线中位数，带为患者曲线IQR，同一特征跨方法共用纵轴。10为各患者状态比例拟合后跨患者等权均值；每个秩点四状态之和为1。所有分箱与拟合值保留服务器。
04–05及09的统计使用原未平滑的患者Spearman rho，与显示曲线的峰值不作互换。原轨迹汇总采用双侧signed-rank检验（原Python实现默认exact/auto口径），每方法在12个原预定义关联终点内BH。连续相关通常20位可估计患者，CNV仅19位；原48项汇总及956条患者相关复算通过。06–08仅作RNA动态描述，不新增基因轨迹推断性检验。

讨论与解读范围
横截面拟时序是计算状态排序，不是随访时间、临床早晚期或谱系追踪；不能将原方案“低恶性到高恶性”的预期写成已经证实的事实。锚点用于root或方向选择，锚点相关属于方法诊断，不是独立的恶性进展验证。
CytoTRACE2 Score代表潜能。CellRank使用了拟时序与既定表达终态，没有spliced/unspliced层，箭头不是RNA velocity；TYMP-state命运概率不能独立验证TYMP导致终末状态。不同方法结果不一致的主要证据全部保留。
RNA曲线不代表蛋白活性、代谢通量或代谢物含量。scFEA的MTA缺失和胆碱局部代理范围与2.6一致。患者IQR不是置信区间，细胞和秩箱不作为独立样本量。

图件索引
'''
for sec,note in notes.items():
 for r in reg:
  if r['section']==sec:note+=r['name']+' | '+r['title']+'\n'
 note+='\n完整统计（原精度数值；q/FDR为对应家族校正值）\n'
 if sec=='2_6':
  for prefix in ['four_pathway_AUCell','scFEA_four_axis']:
   for suffix in ['paired_statistics_v1.tsv','correlation_summary_v1.tsv']:
    p=b/'metabolism_v1/tables'/(prefix+'_'+suffix);note+='\n'+p.name+'\n'+fulltable(p)
  for name in ['ENZYME_F1_STATISTICS.tsv','ALL85_TYMP_PATIENT_STATISTICS.tsv']:note+='\n'+name+'\n'+fulltable(q/name)
 else:note+=fulltable(b/'trajectory_v1/tables/trajectory_all_methods_correlation_summary_v1.tsv')
 note+='\n唯一分析与进度目录：'+str(w)+'\n所有脚本、逐细胞数据、完整患者表、输入SHA256与审计留存服务器。用户审阅后继续。\n'
 (w/('figures_'+sec)/('README_'+sec.replace('_','.')+'_阅读说明.txt')).write_text(note,encoding='utf8')
# Hash exact original source files and the extraction source; stream large files.
paths=[Path(p) for p in json.loads((q/'SOURCE_PATHS.json').read_text())]
paths.append(Path(json.loads((q/'RNA_EXTRACTION_VALIDATION.json').read_text())['source']))
paths += [b/'trajectory_v1/source_data/monocle3_principal_graph_segments_v1.tsv',b/'trajectory_v1/source_data/cellrank_pseudotime_transition_vectors_v1.tsv.gz']
paths += list((w/'scripts').glob('*'))
paths += [w.parent/'local_records_migration_20260908_150455/project/phase24_part1_analysis_completion_audit_local_v1/extracted_v1/checklist_structure_v1.txt']
manifest=[]
for p in sorted(set(paths)):
 before=p.stat();h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
 after=p.stat();assert before.st_size==after.st_size and before.st_mtime_ns==after.st_mtime_ns
 manifest.append({'path':str(p),'bytes':before.st_size,'mtime_ns':before.st_mtime_ns,'sha256':h.hexdigest()})
 print('HASHED',p.name,before.st_size,flush=True)
(q/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
rows=[]
for sec in ['2_6','2_7']:
 for p in sorted((w/('figures_'+sec)).iterdir()):
  if p.is_file():rows.append({'section':sec,'name':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
assert len(rows)==40
(q/'DELIVERY_MANIFEST.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2))
print('DELIVERY_PREPARED',len(rows),sum(r['bytes'] for r in rows),flush=True)
