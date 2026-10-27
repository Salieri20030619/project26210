from pathlib import Path
import json,hashlib,subprocess,xml.etree.ElementTree as ET
w=Path(__file__).resolve().parents[1];q=w/'qc';root=w.parent
v=json.loads((q/'VALIDATION.json').read_text());assert v['status']=='PASS'
registry=json.loads((q/'FIGURE_REGISTRY.json').read_text());assert len(registry)==3
note='''2.7 拟时序渐变热图补图 v1

本包补齐原方案image19所示的平滑渐变热图。建议先看11_Monocle3_gradient_heatmap_v1；12为Slingshot对照，13为scTour对照。每张提供360 dpi PNG及PDF，加本说明共7文件。

图形结构
主体为沿拟时序秩排列的表达趋势热图，采用蓝紫—紫色—橙色—黄色渐变。左侧为基因聚类树、C1–C4色条及基因数；上方为真实患者等权细胞状态构成和早→晚拟时序色带。右侧给出完整图例，并以引线标出12个指定酶/转运体在实际行中的位置。
黄色表示该基因的拟合表达高于其自身沿轨迹的均值，紫蓝表示相对低值；不能用颜色直接比较不同基因的绝对表达量。

数据与方法
使用原正式轨迹输入的2,500个HVG，补入其中未含的11个原方案指定酶/转运体：TYMP、SLC29A1、SLC29A2、ODC1、SRM、SMS、MTAP、CHKA、PCYT1A、PCYT2、CEPT1。TK1已在原HVG中。总计2,511基因全部保留，没有按图形效果挑选动态基因，也未筛选细胞层P值。
HVG原始计数来自正式trajectory_v1输入，按正式全保留基因的细胞文库UMI归一化为log1p(CP10K)；补充基因来自Phase129经核对的同一源表达矩阵。
原对象11,950细胞、26患者；趋势汇总要求每患者至少30恶性细胞，因此使用20位合格患者的11,921细胞。另6位小样本患者共29细胞未进入趋势拟合，源对象保留。
每种方法分别保留原患者内拟时序，不重新拟合Monocle3/Slingshot/scTour或改动方向。每位患者内按拟时序秩分20箱，箱内求基因均值，以箱细胞数作权重，采用固定Gaussian带宽0.075秩单位拟合；随后跨20位患者等权平均。
热图原始横轴为0.025–0.975的200个拟合秩位置，各基因在该方法内沿秩作z标准化。渲染时仅在水平方向插值至1,400位置以形成连续渐变，不在基因行之间平滑。三个方法共用[-2.5,2.5]色阶，超界值按端点饱和，完整原值仍保留。
Monocle3的行z曲线使用Ward/Euclidean层次聚类，固定切为4簇。仅旋转原树分支以便阅读早晚峰值，未改动树距离或基因值。C1–C4分别含408/741/340/1,022基因。三方法完全共用Monocle3的基因行顺序和聚类标签，便于逐行比较。
顶部状态带来自原Base、Hypoxia、Cycle、TYMP_state标签，采用相同患者内分箱/平滑后跨患者等权平均，每个秩点四状态之和为1。它显示构成比例，没有把一个共同秩位置虚构成单一细胞类型。

结果与解读
热图显示不同基因的拟合表达峰值分布在早、中、晚不同秩区间；指定酶基因的位置均由真实行索引引出。Slingshot和scTour对照保留自身动态，不按Monocle3结果重新排序或反转方向。
本图用于描述和比较表达动态，不增加推断性检验。2.7之前的患者层相关性、P/FDR和方法间不一致结果保持原结论。

讨论范围
z标准化突出相对变化形状，不代表绝对效应很大，也不代表这些基因已被证明为显著动态基因。C1–C4是基因曲线显示簇，不是新的恶性细胞亚型。列是拟合秩位置，不是独立生物学样本。
拟时序是横截面细胞状态排序，不能直接解释为临床早晚期或真实纵向时间；这些渐变不证明TYMP因果驱动或跨方法一致的恶性进展。RNA表达不等于蛋白活性或代谢通量。

核查
细胞/基因身份、归一化文库量、患者纳入和12基因覆盖核对通过。独立按Gaussian加权和公式复算拟合，最大误差2.7e-15。三个方法基因行顺序相同，状态比例之和为1。PDF单页及视觉检查完成后交付。
唯一分析记录目录：'''+str(w)+'\n原值、全部患者曲线、拟合矩阵、行簇清单、输入哈希、软件与脚本只保存在服务器。\n'
(w/'figures/README_2.7_渐变热图.txt').write_text(note,encoding='utf8')
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
(q/'PDF_VALIDATION.json').write_text(json.dumps(checks,indent=2));assert len(checks)==3 and all(not r['outside_page_text'] for r in checks),checks
b=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/09_metabolism_trajectory/trajectory_v1');old=root/'phase129_sections2_6_2_7_redraw_remote_v1'
paths=[b/'input/genes_v1.tsv',b/'input/cells_v1.tsv',b/'input/raw_counts_genes_by_cells_v1.mtx.gz',b/'input/cell_metadata_v1.tsv.gz',b/'source_data/trajectory_all_methods_cell_scores_v1.tsv.gz',old/'qc/RNA_EXTRACTION_ROWS.tsv',old/'data/ENZYME_LOG1P_CP10K.tsv.gz',root/'local_records_migration_20260908_150455/references/image19.png']+list((w/'scripts').glob('*'))
manifest=[]
for p in sorted(set(paths)):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
 manifest.append({'path':str(p),'bytes':p.stat().st_size,'sha256':h.hexdigest()})
(q/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
rows=[{'name':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted((w/'figures').iterdir()) if p.is_file()]
assert len(rows)==7
(q/'DELIVERY_MANIFEST.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2))
print('DELIVERY_PREPARED',len(rows),flush=True)
