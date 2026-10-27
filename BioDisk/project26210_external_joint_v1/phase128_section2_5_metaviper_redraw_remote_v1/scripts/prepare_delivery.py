from pathlib import Path
import json,hashlib,subprocess
import pandas as pd
w=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase128_section2_5_metaviper_redraw_remote_v1');q=w/'qc'
note="2.5 MR–metaVIPER 与 F1 / 拟时序图件（v1）\n\n本包包含10套图件，每套提供360 dpi PNG及PDF，加本说明共21个文件。建议先看01、02、03、09；04为原始活性UMAP，05–08为拟时序趋势展示，10为明确分开的RNA补充图。\n\n结果\nCASP9、CHMP1A和PRSS8的推断活性与F1呈正向患者层关联；PXN呈负向或当前数据不支持关联。连续F1及患者内四分位比较均以20位患者为统计单位；主分组比较为9位配对患者，按探索性结果呈现。\n四蛋白MR效应与患者等权平均活性的散点呈负向描述性关系（Pearson r=-0.6483），当前数据不支持“MR效应越强，单细胞推断活性越高”的主张。固定这4个蛋白后，在20位患者内重复计算的Pearson相关均为负，中位r=-0.6637，P=1.90735e-6，BH FDR=3.81470e-6。该敏感性结果不增加MR蛋白数量。\nCASP9、CHMP1A、PRSS8的末端共同升高未在Monocle3、Slingshot、scTour三种方法间重复。Slingshot中的部分正向结果与其他方法不同，完整保留。PXN的正向末端变化仅在Monocle3获得原FDR支持。\n8个最终MR蛋白中，CASP9、CHMP1A、PRSS8、PXN有合格肾癌调控网络；KDELC2（映射POGLUT3）、RRM2B、TMEM25、TYMP缺少合格网络，不能生成其metaVIPER活性。因此TYMP/RRM2B只提供原注册的RNA敏感性图。\nHIF1A是单独的背景靶点，不属于这4个最终MR靶点。其metaVIPER活性与F1为负向；这与2.4中pySCENIC调控子未检出是不同方法的结果，不互相替代。\n\n图件说明\n01 A：横轴为最终MR表报告OR的自然对数，纵轴为患者等权平均metaVIPER NES。四个点对应真实蛋白；不制作四点密度等高线，不标注四蛋白推断性P值。右侧每点为一位患者在相同四蛋白上的Pearson相关，黑菱形为中位数。\n02 B：主分组是cNMF中F1 usage最大的细胞与其他主导模块细胞，不是简单按全局F1中位数切分。小提琴仅使用9位合格配对患者的7,021个细胞（F1主导4,069，其他2,952）；两组各患者均至少10细胞。连线连接同一患者的组中位数，避免用细胞数替代样本量。\n03 C：并列展示患者内F1连续相关、F1主导组差值及患者内上/下四分位差值。不同终点使用各自的效应单位。n=20面板保留效应、P和FDR；n=9面板优先展示效应和样本一致性，完整统计在本说明及服务器表中。\n04 D：11,921个正式可评分细胞的原始metaVIPER NES。保留原11,950细胞UMAP作为坐标背景，另29个无正式评分细胞仅以灰色显示，不补值。所有色阶以NES=0为中性颜色；正负两侧各按完整数据范围映射。\n05 E：四个MR蛋白在三种拟时序中的患者平衡趋势曲线；同一蛋白的三个面板共享纵轴。线为各患者固定核拟合后的等权均值，带为患者曲线IQR，不能称为置信区间；小点为原患者平衡分箱均值。\n06–08 F：分别为Monocle3、Slingshot、scTour。先在每位患者内根据其拟时序秩拟合活性，再把该患者每个细胞的拟合值映射到冻结UMAP。属于显示用的拟合投射，不是原始单细胞活性。每个蛋白跨三方法共用数值上下限且NES=0居中；不对UMAP空间做人工平滑或添加轨迹箭头。\n09 G：左图为患者内末端20%减初始20%的活性中位数差及正向人数，右图为三个蛋白对同时正向变化的患者比例。效应量、P和FDR保留。右图仅描述共同方向，不进行细胞层检验。\n10 S：TYMP、RRM2B原始RNA的log1p CP10K轨迹敏感性；采用同一显示规则。不能将这些曲线称为蛋白活性、酶活或救援后的metaVIPER结果。\n\n方法\n复用通过终末验证的5,566调控因子×11,921细胞metaVIPER矩阵，包含20位患者。此次直接从正式RDS提取4个MR蛋白及HIF1A，与历史导出逐值一致。cNMF标签、MR表坐标、原轨迹及全部统计均按精确细胞/患者键核对。\n正式统计使用原注册的患者单位：组差每臂至少10细胞；连续关联至少20细胞并有变异；四分位由患者自己的F1分布计算。拟时序末端/初始取患者内80%/20%分位，每臂至少5细胞。\n显示平滑使用原20秩箱的实际患者均值与箱细胞数，Gaussian核带宽固定为0.075个患者拟时序秩单位（1.5箱）。患者内部按箱细胞数加权，再跨患者等权汇总；参数对全部蛋白、RNA及三种方法一致。原始患者效应和统计不由平滑曲线重新计算。\nBH家族沿用原注册定义：F1每个终点在4个最终MR蛋白内校正；HIF1A单独计算；拟时序每种方法及证据类型分别校正；RNA与MR活性分开；MR患者敏感性在Pearson/Spearman两种指标内校正。CytoTRACE2与CellRank只作为原有潜能/命运概率敏感性，未替换为主拟时序。\n\n讨论与解读范围\nmetaVIPER NES是表达与调控网络推断的相对调控活性，不是直接测量的蛋白含量或生化酶活；F1关联不建立TYMP的因果驱动作用。F1主导组的9位患者按探索性定位，其完整P/FDR保留，不因图形美观改变统计口径。\n跨蛋白关系只涵盖4个可推断靶点，不代表完整MR蛋白空间。患者内重复的敏感性分析仍以同一4个蛋白为条件。\n拟时序是计算排序，不等同于随访时间或临床早晚期。拟合UMAP保留患者基线差异，其区域颜色不能单独证明状态演进或终末共激活。当前数据不支持稳定的跨方法末端共同升高，也不能由RNA替代缺失蛋白活性得出该结论。\n\n复核\n300条F1效应、700条轨迹相关、420条末端差值复算通过；原分箱均值最大误差约8e-15；原BH检验家族复算通过。10张最终PDF均完成单页结构与渲染视觉检查。旧排版初稿保存在服务器。\n"
m12=pd.read_csv(q/'m12_f1_activity_patient_summary_v1.tsv',sep='\t')
note+='\nF1终点完整统计\n'
for r in m12.itertuples():note+=f'{r.protein} | {r.endpoint} | n={r.n_eligible_patients} | effect={r.median_patient_effect:.7g} | P={r.wilcoxon_p_two_sided:.7g} | FDR={r.BH_FDR:.7g} | 正/负/零={r.n_positive}/{r.n_negative}/{r.n_zero}\n'
ter=pd.read_csv(q/'m13_terminal_shift_summary_v1.tsv',sep='\t')
note+='\n末端差值完整统计（RNA与活性分别标注）\n'
for r in ter.itertuples():note+=f'{r.feature} | {r.feature_type} | {r.method} | n={r.n_eligible_patients} | delta={r.median_patient_effect:.7g} | P={r.wilcoxon_p_two_sided:.7g} | FDR={r.BH_FDR:.7g} | 正/负/零={r.n_positive}/{r.n_negative}/{r.n_zero}\n'
note+='\n唯一记录目录：'+str(w)+'\n脚本、完整统计、逐细胞原值/拟合值、来源SHA256及进度只写服务器。本地交付目录：C:/Users/31937/Downloads/26210/06_2.5_metaVIPER_v1。\n'
(w/'figures/README_2.5_metaVIPER.txt').write_text(note,encoding='utf8')
checks=[]
for f in sorted((w/'figures').glob('*.pdf')):
 info=subprocess.run(['pdfinfo',str(f)],capture_output=True,text=True,check=True).stdout
 assert any(l.startswith('Pages:') and l.split(':')[1].strip()=='1' for l in info.splitlines())
 checks.append({'name':f.name,'pages':1,'status':'PASS'})
assert len(checks)==10
(q/'PDF_VALIDATION.json').write_text(json.dumps(checks,indent=2))
(q/'VISUAL_REVIEW.json').write_text(json.dumps({'status':'PASS','final_pdf_renderings_inspected':10,'render_scale':1700,'refinements':['Violin title spacing','RNA footer spacing','Shared within-protein curve y-axis','NES zero-centred colour scales','Sparse fitted colourbar ticks'],'labels':'No clipping or overlap in final reviewed renderings','user_review':'PENDING'},indent=2))
sources=json.loads((q/'SOURCE_MANIFEST.json').read_text());known={r['path'] for r in sources}
for f in (w/'scripts').iterdir():
 if f.is_file() and str(f) not in known:sources.append({'path':str(f),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
(q/'SOURCE_MANIFEST.json').write_text(json.dumps(sources,indent=2))
rows=[{'name':f.name,'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in sorted((w/'figures').iterdir()) if f.is_file()]
assert len(rows)==21
(q/'DELIVERY_MANIFEST.json').write_text(json.dumps(rows,indent=2))
print(json.dumps({'files':len(rows),'bytes':sum(r['bytes'] for r in rows),'PDF_checks':10,'visual_review':'PASS'}))
