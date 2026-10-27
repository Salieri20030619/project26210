from pathlib import Path
import json,hashlib,datetime
w=Path(__file__).resolve().parents[1];q=w/'qc'
b=json.loads((q/'VALIDATION.json').read_text());e=json.loads((q/'E_NETWORK_VALIDATION.json').read_text());pdf=json.loads((q/'PDF_VALIDATION.json').read_text())
assert b['status']==e['status']=='PASS' and len(pdf)==5
note='''2.4 B/E 修正版 v2

本包5套图：B热图、E基因连接总图、E1 ILF2、E2 HOXD3、E3 HNF4A。每套提供360 dpi PNG和PDF，加本说明共11文件。用于替换旧2.4包中的B和E；其余面板不变。

B：补齐缺失的标记和图例
左侧新增原行聚类树的4簇显示标记、色条和每簇调控子数：C1(9)、C2(8)、C3(4)、C4(3)。编号按热图从上到下排列，由原average-linkage/Euclidean树按k=4切分，仅用于图形注释。
顶部新增明确的Cell type、cNMF program、Algorithm三行：细胞类型统一为Malignant epithelial；F1–F7为cNMF程序；G3/GB2对应GENIE3/GRNBoost2。右侧给出所有颜色、算法、聚类数目和加粗候选基因的图例。
原24调控子、14列、行列顺序、患者平衡AUCell汇总、方法内行z标准化均保持一致，原值复算最大误差<1e-15。21位患者贡献37份患者×程序谱；F1–F7独立患者数仍为9/12/2/4/4/2/4。同一患者可贡献多个程序，不能相加当作独立样本量。未新增推断性检验。

E：纠正旧图将方法画成连线端点的错误
旧E的连线端点为GENIE3/GRNBoost2和靶基因，虽然表格记录的是TF→靶基因边，但图形表达成了方法→基因，不能作为基因调控网络使用。当前版本将每条线的起点落实到TF基因，终点落实到靶基因；网络节点全部为基因。
总图包含ILF2、HOXD3、HNF4A三位TF和68个原展示靶基因，共71个基因节点。保留原68条唯一TF→靶基因连接，并从正式完整边表补回同一节点集合内真实存在的3条交叉边：ILF2→PTER、HNF4A→PTN、HNF4A→TRPV4。总计71条唯一有向基因边，对应89条原方法支持记录。
边颜色表示方法支持：紫色为双算法共有（18条），棕色为GENIE3特有（20条），蓝色为GRNBoost2特有（33条）。算法名称只在图例中作为边证据类型，不作为网络节点。两算法权重尺度不同，未将原始权重混作跨方法可比的线宽。
菱形为TF基因；普通圆点为靶基因；金色圆点PTER、PTN、TRPV4表示被多个展示TF共同连接的靶基因。单独网络分别含ILF2的15个靶基因、HOXD3的24个靶基因、HNF4A的32个靶基因；总图中的共享靶基因只保留一个节点。
展示节点沿用旧图的固定选择：每个TF按方法支持类别，各取方法内权重百分位均值最高的最多10个靶基因；修正版在这个既定节点集合上补齐所有真实存在的边。因此图为可读性选择的子网络，不是全部调控网络。全部非自身靶点数原为GENIE3/GRNBoost2：ILF2 4/1133、HOXD3 14/53、HNF4A 87/518。未为视觉效果添加无依据的靶基因间连线，也未重新推断网络。
这些是pySCENIC预测的调控子边，不能直接称为实验确认的物理相互作用或因果调控。现有ILF2→TYMP仅有GRNBoost2预测，未进入既定展示节点；HOXD3/HNF4A没有TYMP边，不新增直接调控TYMP的主张。

核查
B矩阵、标签和患者数复核通过；E全部71条唯一基因边逐一对应原边表，原展示68边全部保留，算法节点数为0。5张PDF均为单页，页面边界检查及逐图视觉审阅通过。
唯一记录目录：'''+str(w)+'\n脚本、完整矩阵、节点/边清单、来源SHA256、核查及进度仅保存在服务器。\n'
(w/'figures/README_2.4_B_E_修正版.txt').write_text(note,encoding='utf8')
(q/'VISUAL_REVIEW.json').write_text(json.dumps({'status':'PASS','reviewed_at':datetime.datetime.now().isoformat(),'PDF_renderings_inspected':[r['name'] for r in pdf],'rendering':'Poppler 125 dpi','B_refinement':'Moved Cluster heading clear of Algorithm annotation label','E_review':'Every endpoint is a gene; arrows point from TF to target; all gene labels and legends readable','user_review':'PENDING'},indent=2))
manifest=[]
for f in sorted((w/'figures').iterdir()):
 if f.is_file():manifest.append({'name':f.name,'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
assert len(manifest)==11
(q/'DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
print('DELIVERY_READY',len(manifest),sum(r['bytes'] for r in manifest),flush=True)
