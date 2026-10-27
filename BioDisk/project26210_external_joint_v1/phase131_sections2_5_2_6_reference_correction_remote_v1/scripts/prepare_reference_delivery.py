from pathlib import Path
import json,hashlib,subprocess,xml.etree.ElementTree as ET
import pandas as pd
w=Path(__file__).resolve().parents[1];q=w/'qc';root=w.parent
v25=json.loads((q/'SECTION25_VALIDATION.json').read_text());v26=json.loads((q/'SECTION26_VALIDATION.json').read_text());assert v25['status']==v26['status']=='PASS'
old25=root/'phase128_section2_5_metaviper_redraw_remote_v1';old26=root/'phase129_sections2_6_2_7_redraw_remote_v1'
note25='''2.5 密度等高线相关图修正版 v2

按用户要求，两种关系均绘制，并严格区分：
01：MR效应 × metaVIPER活性，4个蛋白点及其二维Gaussian KDE等高线。
02：F1 × 各蛋白metaVIPER活性，CASP9/CHMP1A/PRSS8/PXN四个独立面板，密度等高线和患者层统计分开标注。
每张提供360 dpi PNG和PDF；加本说明共5文件。

参考对应
原方案2.5第1项后接image13与image14，两者都是二维核密度填色加等高线，非环形弦图。此次重绘使用这一图形结构，替换旧01_A普通散点作为本要求的主图；原其他蛋白活性分组/拟时序图件仍保留。

结果
MR效应与患者等权平均活性的Pearson r=-0.648335，四个实际点为CASP9、CHMP1A、PRSS8、PXN。当前方向为负，不能写成MR效应越强活性越高。
F1与CASP9、CHMP1A、PRSS8的患者内活性相关为正并通过原四MR蛋白家族FDR；PXN的患者内中位相关为负，FDR=0.053169。每蛋白对应的完整效应、P/FDR保留如下。

方法
01横轴使用原MR报告OR的自然对数，不误标为log HR。纵轴为先患者内平均活性、再跨20患者等权平均的metaVIPER NES。每个蛋白一个点，4点等权Gaussian KDE，Scott带宽。图中保留实际点和蛋白名，展示完整等高线；带宽0.7/1/1.4倍敏感性图留服务器QC。
02使用原11,921细胞及20患者的正式活性矩阵，F1与精确细胞ID逐一对齐。Gaussian KDE中每位患者的总权重相同，患者内每细胞等权；四蛋白使用同一带宽规则。色阶为各蛋白内部的相对密度，等高线水平固定。Weighted cell r使用与KDE相同权重计算Pearson相关，作为分布描述；右侧独立标注原患者内Spearman相关汇总、P和FDR。
观察到的F1范围是0–1；图内虚线标出该区间。为展示完整密度等高线，绘图区延伸至-0.2至1.2。区间外仅为Gaussian核尾部，不新增负F1或大于1的细胞观测。

讨论与解读范围
01只有4个可评估MR蛋白，等高线是这4点的平滑显示，不能当成大量独立蛋白观测；不进行4蛋白推断性相关检验。20位患者不增加MR蛋白数量。
02的细胞描述性相关与患者内相关是不同统计量；患者等权密度仍保留患者间基线差异，因此PXN的加权细胞r为正而患者内中位rho为负，两项结果均展示。不得用细胞r取代患者层检验。
密度色阶不表示P值或显著性。TYMP/RRM2B/KDELC2/POGLUT3/TMEM25仍缺合格网络，本次不补造活性、不以RNA冒充蛋白活性，也不增加临床或因果结论。

完整统计
'''
note25+=pd.read_csv(q/'F1_DENSITY_STATISTICS.tsv',sep='\t').to_csv(sep='\t',index=False,float_format='%.10g')
note26='''2.6 原示例布局修正版 v2

本包7套图，每套360 dpi PNG+PDF，加本说明共15文件。本轮修正图形与参考不对应的问题，正式评分、细胞身份、UMAP和原推断性统计均保留。

逐图对应
01 对应image16：横排四张UMAP、带坐标和浅网格、magma紫橙黄色阶及侧边小色条，展示四条原AUCell通路。
02 对应image15：TYMP-/TYMP+双分面气泡矩阵，纵轴四条代谢通路，横轴F1–F7 cNMF程序，灰色分面标题、网格和颜色/大小独立图例。
03 同image15布局，完整展示85条通路，未按结果方向选择保留通路。
04 对应正文Corrplot要求：F1与四通路的圆点相关矩阵，圆面积=|rho|、红蓝颜色=方向和强度，旁列原4项F1关联统计。
05 对应image15气泡图结构及正文12酶/转运体要求：基因为纵轴，F1患者内四分位为横轴，气泡颜色为RNA均值，大小为检出比例。
06 采用image16四联UMAP形式，补充展示四条scFEA代理评分。
07 对应image15显著性气泡的统计编码：四通路×两种患者配对比较，颜色为效应、面积为-log10(P)、黑色轮廓为FDR<0.05。

结果
既有患者层结论保持：TYMP阳性和F1高主克隆均支持嘧啶、Cys/Met、Arg/Pro三通路评分升高；当前数据不支持甘油磷脂的这两项配对差异。20患者的F1连续关联也支持前三条通路，甘油磷脂未获支持。
双分面程序内气泡图为探索性描述，展示程序之间和患者间的真实差异；不从这些小样本分面作新的显著性判断。F3没有符合两臂细胞阈值的患者配对，明确以灰叉标记，不填造数值。

方法
全部11,950恶性细胞、26患者、5队列的原UMAP及评分冻结。01与06只改变绘图布局和色阶，完整显示每个指标的原始取值范围。
02–03要求同一患者、同一cNMF程序的TYMP-/TYMP+两臂各至少10细胞；每臂先求患者均值/比例，再跨相同患者等权平均。F1–F7配对患者数为6/9/0/2/3/2/2。气泡颜色=患者等权平均AUCell，大小=患者等权AUCell>0细胞比例。这里没有使用细胞P值；图例明确写比例，不能套用参考图的配受体P值含义。
04在每位至少30细胞且变量有变异的患者内求Spearman，再取患者rho中位数；共20位。F1与四通路的原P/FDR不变，四项为原预定义同类家族。新增通路之间矩阵项仅描述，不新增检验。
05复用原17位患者的四分位展示表；F1同分值保留，各四分位至少5细胞。RNA均值不按基因重新缩放。
07使用原14位合格患者的配对结果，各组每患者至少20细胞；双侧近似Wilcoxon及原四通路BH家族不变，两种比较分别校正。大小依据有效患者层P，非细胞层P。完整P/FDR和效应如下。

讨论与解读范围
气泡分面的cNMF程序是恶性细胞内程序，不能直接改称参考中的正常/炎症细胞类型。AUCell>0表示该算法评分大于零，不等同于通路显著激活；因为四条重点通路几乎均为正分，02中气泡面积接近是实际结果。07另用有效患者层统计表现差异证据。
程序内n<10按探索性结果呈现，不给伪重复P值。颜色显示转录组通路评分或RNA，scFEA是预测代理；MTA不是显式模型输出，不能由这些图验证实测代谢物44倍变化或TYMP因果作用。

完整患者层统计
'''
b=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/09_metabolism_trajectory/metabolism_v1')
for n in ['four_pathway_AUCell_paired_statistics_v1.tsv','four_pathway_AUCell_correlation_summary_v1.tsv']:
 note26+='\n'+n+'\n'+pd.read_csv(b/'tables'/n,sep='\t').to_csv(sep='\t',index=False,float_format='%.10g')
for sec,note in [('2_5',note25),('2_6',note26)]:
 note+='\n唯一分析记录目录：'+str(w)+'\n所有脚本、患者/细胞表、密度网格、统计、输入SHA256和审计只保存在服务器。\n'
 (w/('figures_'+sec)/('README_'+sec.replace('_','.')+'_参考图修正版.txt')).write_text(note,encoding='utf8')
checks=[]
for sec in ['2_5','2_6']:
 for p in sorted((w/('figures_'+sec)).glob('*.pdf')):
  info=subprocess.run(['pdfinfo',str(p)],capture_output=True,text=True,check=True).stdout
  assert any(x.startswith('Pages:') and x.split(':')[1].strip()=='1' for x in info.splitlines())
  xml=subprocess.run(['pdftotext','-bbox',str(p),'-'],capture_output=True,text=True,check=True).stdout;tree=ET.fromstring(xml);bad=[]
  for page in tree.iter('{http://www.w3.org/1999/xhtml}page'):
   width=float(page.attrib['width']);height=float(page.attrib['height'])
   for word in page.iter('{http://www.w3.org/1999/xhtml}word'):
    d=word.attrib
    if float(d['xMin'])<-.5 or float(d['yMin'])<-.5 or float(d['xMax'])>width+.5 or float(d['yMax'])>height+.5:bad.append(word.text)
  checks.append({'section':sec,'name':p.name,'pages':1,'outside_page_text':bad})
(q/'PDF_VALIDATION.json').write_text(json.dumps(checks,indent=2));assert len(checks)==9 and all(not c['outside_page_text'] for c in checks),checks
paths=list((w/'scripts').glob('*'))+[old25/'data/SELECTED_METAVIPER_ACTIVITY.tsv',old25/'qc/m11_mr_effect_activity_points_v1.tsv',old25/'qc/m12_f1_activity_patient_summary_v1.tsv',old26/'qc/ENZYME_DOTPLOT_PATIENT_PROFILES.tsv',b/'source_data/four_pathway_AUCell_cell_scores_v1.tsv.gz',b/'source_data/KEGG_AUCell_all_pathways_v1.tsv.gz',b/'source_data/scFEA_four_axis_cell_scores_v1.tsv.gz',b/'tables/four_pathway_AUCell_paired_statistics_v1.tsv',b/'tables/four_pathway_AUCell_correlation_summary_v1.tsv',b.parent/'trajectory_v1/input/scvi_umap_v1.tsv.gz',b.parent.parent/'07_pyscenic_v1/input/cell_metadata_v1.tsv.gz']
paths += [root/'local_records_migration_20260908_150455/references'/f'image{i}.png' for i in [13,14,15,16]]
paths += [root/'local_records_migration_20260908_150455/project/phase24_part1_analysis_completion_audit_local_v1/extracted_v1/checklist_structure_v1.txt']
manifest=[]
for p in sorted(set(paths)):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
 manifest.append({'path':str(p),'bytes':p.stat().st_size,'sha256':h.hexdigest()})
(q/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
rows=[]
for sec in ['2_5','2_6']:
 for p in sorted((w/('figures_'+sec)).iterdir()):
  if p.is_file():rows.append({'section':sec,'name':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
assert len(rows)==20
(q/'DELIVERY_MANIFEST.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2))
print('DELIVERY_PREPARED',len(rows),'PDF',len(checks),flush=True)
