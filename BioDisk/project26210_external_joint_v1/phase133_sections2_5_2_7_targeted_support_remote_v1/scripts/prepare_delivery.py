from pathlib import Path
import json,hashlib,subprocess,xml.etree.ElementTree as ET,datetime
import pandas as pd
w=Path(__file__).resolve().parents[1];q=w/'qc';root=w.parent
assert (q/'R_INDEPENDENT_VALIDATION.txt').read_text().startswith('PASS')
assert json.loads((q/'PRIOR_AXIS_RECHECK.json').read_text())['status']=='PASS'
a=pd.read_csv(q/'TYMP_ACTIVITY_SUMMARY.tsv',sep='\t');b=pd.read_csv(q/'TERMINAL_ENZYME_SUMMARY.tsv',sep='\t');old=pd.read_csv(q/'PRIOR_ENZYME_AXIS_COMPLETE_SUMMARY.tsv',sep='\t')
note='''2.5 / 2.7 支持性补充及正文调整 v1

本轮按用户最新定位：湿实验承担主要验证，生信提供方向、患者一致性和名义显著的对应线索；不要求生信独立完成因果机制闭环。2.3 scDRS从正文主线和主图清单移出，其原始结果、图件和阴性记录完整保留。正文不再声称获得scDRS遗传支持。

本包4套图，均有360 dpi PNG/PDF，加本说明共9文件。
01：TYMP RNA与四个MR靶点推断活性的患者内关联，以及技术协变量校正；RRM2B作为单独RNA端点。
02：13个指定酶基因×三种拟时序的全部39项早段—末段pseudobulk RNA比较。
03：02中全部4项名义显著端点的患者配对图，同时包含正向和反向。
04：恢复展示已有Phase30严格/广义酶功能代理的完整轴结果；旧P/FDR不变，名义线索不再被原“高等级支持”标签掩盖。

可以放入结果主线的支持性线索
一、TYMP表达与部分MR靶点活性正向对应。
20位患者的未调整患者内关联中，CASP9中位rho=0.076646，18/20正向，P=0.000636、FDR=0.001668；CHMP1A rho=0.117349，16/20正向，P=0.001825、FDR=0.002434；PRSS8 rho=0.110424，17/20正向，P=0.000834、FDR=0.001668。这比跨4个蛋白比较MR效应大小更直接对准TYMP相关状态。
二、既有TYMP严格功能代理在F1高组中名义升高。
F1高主克隆相对其他主克隆：n=14，中位delta=0.057178，11正/2负/1零，P=0.027708、FDR=0.055416。患者内F1 Q4相对Q1：n=20，delta=0.042944，13正/5负/2零，P=0.031104、FDR=0.062208。均作为名义显著支持，而非校正后显著。
三、TYMP广义功能代理在两种拟时序末段名义升高。
Slingshot：n=20，delta=0.162888，15/20正向，P=0.044054、FDR=0.088108。
scTour：n=20，delta=0.300574，14/20正向，P=0.044054、FDR=0.088108。
这两项是原表已有的结果，本轮恢复其展示，不声称新独立验证。它们支持部分功能程序随状态排序变化。
四、新增RNA端点提供方法特异性对应。
scTour下TYMP：n=20，末段减早段的中位delta=+0.312972 log2(CPM+1)，13/20正向，P=0.034919、39项FDR=0.340461。
Slingshot下SRM：n=20，delta=+0.635434，13正/6负/1零，P=0.025520、FDR=0.331758。
两项均为探索性、名义显著，不能合并成同一条轨迹或一致共激活。逐一留出患者后，中位效应范围分别为0.253763–0.372180、0.540414–0.730454，方向均保持正向。

主要不一致和敏感性结果
技术调整后，CASP9和CHMP1A关联减弱；PRSS8保留名义正向（rho=0.054238，P=0.045793、FDR=0.091586）；PXN从未调整的近零负向变为小幅正向（rho=0.046029，16/20正向，P=0.004823、FDR=0.019292）。两套模型均展示，不能只挑有利模型。
TYMP与RRM2B RNA的未调整相关很小（rho=0.024532，13/20正向，单RNA终点P=0.045793），技术调整后接近零（P=0.695063）。这只提供转录层候选线索，不是RRM2B蛋白活性证据。
新增39项RNA比较中，Slingshot的PCYT1A下降（delta=-0.689216，P=0.002091、FDR=0.081553）、scTour的SMS下降（delta=-0.226539，P=0.016042、FDR=0.312820）同样达到名义显著并完整展示。
原TYMP广义代理在Monocle3方向为负且未名义显著；RRM2B广义代理与F1为负向。RRM2B严格代理虽在Slingshot/scTour有名义P，但均只有5位患者正向、14或15位为零，总体中位效应为0；不称为普遍共同激活。

方法
本轮补充方案在新端点计算前固定，属于既有数据上的探索性分析。使用原20位患者、11,921个恶性细胞、4个有效metaVIPER靶点和固定13个酶基因；没有按本轮结果增删基因或反转拟时序。
2.5以患者内Spearman rho为效应，每患者至少30细胞且变量有变异。技术敏感性把TYMP、目标和原文库量/检出基因数分别转秩后，去除两技术协变量的线性部分，再相关。推断检验以患者效应为单位，双侧Wilcoxon正态近似、连续性校正、零差值剔除。四个蛋白在主模型、技术敏感性中分别作为预定义4项BH家族；单RNA主端点不额外作FDR。
2.7保留Monocle3/Slingshot/scTour原患者内拟时序。早段<=20%分位，末段>=80%分位，阈值不同且每臂至少5细胞、患者总细胞>=30。每患者每臂汇总原始基因UMI和全保留基因总UMI，计算log2(1+1e6*gene_UMI/total_UMI)，末段减早段为效应。这个新pseudobulk端点补充原细胞中位数端点，不替换它。13基因×3方法共39项统一BH，所有结果保留。
中位效应区间使用5,000次患者bootstrap；与Wilcoxon检验不完全相同，区间和P原样展示。逐患者留出只检查效应稳定性，没有据此剔除患者。
原Phase30功能代理沿用其冻结定义、患者内z标准化、组内中位数、原P及“终点×方法×feature_class”BH。严格TYMP代理=TYMP RNA、M45、M164的患者内z均值；广义TYMP代理=TYMP RNA与Pyrimidine_flux的z均值。其余轴及组件在服务器原登记保留。
历史字段F1_dominant_group在Phase30代码实际比较F1最高主要克隆与其他主要克隆，本轮图例已按真实定义写为F1-high clone；不与2.5原9患者的F1主导细胞比较混淆。

讨论与对应范围
生信结果按方向、效应、患者一致性和名义证据支持既定实验主线；“名义显著”不改写为FDR显著。同一数据的多种分析不增加独立样本量。
metaVIPER是推断活性，功能轴是RNA/scFEA组合代理；含TYMP RNA的组合不能当作TYMP独立蛋白测量。原MR效应—活性负向描述保留补充，正文不再要求或声称“MR效应越强活性越高”。
拟时序为横截面排序，具体RNA和代理变量的方向可能不同；不把TYMP广义代理的升高改写成所有方法TYMP RNA均升高。
用户湿实验已验证作为本项目叙事依据；本轮尚未收到具体干预和方向清单，因此上述先作为待精确对应的支持线索，不新增未告知的湿实验结果。

正文组织建议
保留2.1、2.2、2.4、2.5、2.6、2.7原章节编号便于追踪；2.3不进入正文主图序列。
2.5优先呈现TYMP与可推断蛋白活性关联、TYMP严格/广义功能代理；技术敏感性和原跨蛋白MR关系保留说明及补充。
2.7优先说明TYMP相关功能代理在Slingshot/scTour中的方向一致和名义证据，再展示scTour TYMP与Slingshot SRM的分开RNA例证。完整三方法39项比较及主要反向结果随补充保留。
后续若需要分支分析，应先恢复真实Slingshot分支权重；现有all_lineages文件使用na=FALSE，不能把它直接当作分支归属。未据此开展无依据的分支检验。

完整统计
'''
for title,frame in [('本轮TYMP关联',a),('本轮39项末端RNA',b),('原功能轴完整汇总',old)]:note+='\n'+title+'\n'+frame.to_csv(sep='\t',index=False,float_format='%.10g')
note+='\n复核：200条患者相关、780条计数推导的末段差值由R独立复算通过；相关最大误差2.51e-15，P最大误差2.23e-16，全部本轮BH通过。原72项汇总的效应、正负零计数及全部BH家族复核通过。\n唯一记录目录：'+str(w)+'\n所有逐细胞/患者数据、脚本、统计和审计仅在服务器。\n'
(w/'figures/README_支持性补充与正文调整.txt').write_text(note,encoding='utf8')
(q/'MAIN_TEXT_SECTION_STATUS.tsv').write_text('section\tmain_text_status\tdecision\n2.1\tretain\tclone heterogeneity and nominal TYMP association\n2.2\tretain\tTYMP-associated transcriptional program\n2.3\tomit\tuser decision; negative scDRS evidence retained in archive\n2.4\tretain\tcandidate regulon correspondence\n2.5\tretain_with_supplement\tdirect TYMP/activity links and clearly labeled enzyme-function proxies\n2.6\tretain\tmetabolic pathway correspondence\n2.7\tretain_with_supplement\tmethod-specific directional and nominal dynamics\n')
checks=[]
for p in sorted((w/'figures').glob('*.pdf')):
 info=subprocess.run(['pdfinfo',str(p)],capture_output=True,text=True,check=True).stdout;assert any(x.startswith('Pages:') and x.split(':')[1].strip()=='1' for x in info.splitlines())
 xml=subprocess.run(['pdftotext','-bbox',str(p),'-'],capture_output=True,text=True,check=True).stdout;tree=ET.fromstring(xml);bad=[]
 for page in tree.iter('{http://www.w3.org/1999/xhtml}page'):
  width=float(page.attrib['width']);height=float(page.attrib['height'])
  for word in page.iter('{http://www.w3.org/1999/xhtml}word'):
   d=word.attrib
   if float(d['xMin'])<-.5 or float(d['yMin'])<-.5 or float(d['xMax'])>width+.5 or float(d['yMax'])>height+.5:bad.append(word.text)
 checks.append({'name':p.name,'pages':1,'outside_page_text':bad})
assert len(checks)==4 and all(not r['outside_page_text'] for r in checks),checks
(q/'PDF_VALIDATION.json').write_text(json.dumps(checks,indent=2))
(q/'VISUAL_REVIEW.json').write_text(json.dumps({'status':'PASS','renderings_inspected':[r['name'] for r in checks],'rendering':'Poppler 125 dpi','checks':['All patient directions and uncertainty visible','Nominal P and FDR separated','Primary and technical sensitivity both shown','RNA and protein-activity types separated','All 39 endpoints and opposing nominal effects retained','Original strict/broad axes and ties explicitly displayed'],'user_review':'PENDING'},indent=2))
paths=[Path(p) for p in json.loads((q/'SOURCE_PATHS.json').read_text())]+list((w/'scripts').glob('*'))+[w/'SUPPLEMENT_ANALYSIS_PLAN.md']
prior=root/'local_records_migration_20260908_150455/project/phase30_part1_enzyme_function_axis_local_v1'
paths += [prior/'results_local_v1/enzyme_axis_patient_effects_v1.tsv',prior/'01_analyze_enzyme_function_axes_v1.py',prior/'PHASE30_ENZYME_AXIS_REGISTRY_v1.md']
manifest=[]
for p in sorted(set(paths)):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
 manifest.append({'path':str(p),'bytes':p.stat().st_size,'sha256':h.hexdigest()})
(q/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
rows=[{'name':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted((w/'figures').iterdir()) if p.is_file()];assert len(rows)==9
(q/'DELIVERY_MANIFEST.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2))
print('DELIVERY_PREPARED',len(rows),flush=True)
