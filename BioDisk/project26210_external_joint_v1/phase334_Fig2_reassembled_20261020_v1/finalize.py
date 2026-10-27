from pathlib import Path
import json,hashlib
import pandas as pd
import pymupdf as fitz
w=Path(__file__).resolve().parent;o=w/'delivery';q=w/'qc'
pdf=o/'Fig2_reassembled_A4_v1.pdf';d=fitz.open(pdf);p=d[0];text=p.get_text()
assert len(d)==1 and abs(p.rect.width*25.4/72-210)<.01 and abs(p.rect.height*25.4/72-297)<.01
nodes=pd.read_csv(w/'inputs/GENE_NETWORK_NODE_POSITIONS.tsv',sep='\t').gene
missing=[g for g in nodes if g not in text];assert not missing,missing
v=json.loads((q/'VALIDATION.json').read_text(encoding='utf-8'))
v['status']='PASS';v['visual_review']='Final A4 render reviewed; horizontal network labels; no clipped colour scale or residual CNV annotation; new A and A-F panel order confirmed.'
v['network_all_71_node_labels_present']=True
v['source_hashes']={x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in (w/'inputs').iterdir() if x.is_file()}
v['pdf_sha256']=hashlib.sha256(pdf.read_bytes()).hexdigest()
(q/'VALIDATION.json').write_text(json.dumps(v,ensure_ascii=False,indent=2),encoding='utf-8')
caption='''图2 TYMP相关恶性上皮程序（A4重组版v1）

A：上皮恶性身份与恶性上皮转录程序。左图在冻结UMAP中展示38,670个肿瘤组织上皮细胞，颜色对应原CNV分类：恶性共同支持10,844个、仅inferCNV支持恶性1,106个、非恶性共同支持498个、方法不一致22,216个、未定4,006个。前两类构成11,950个恶性细胞分析集。右侧三图分别显示该恶性细胞集的七个cNMF主导程序、F1使用量及归一化TYMP RNA；三图使用相同原始UMAP坐标，颜色分别表示程序身份、使用量和log1p表达。保留原99百分位显示上限及全部细胞，不重算嵌入或分群。

B：GSE156632 RCC2的原染色体CNV示例。横轴按染色体1–22排列，标示1、5、10、15、20号；纵轴上方200个参考细胞，下方701个肿瘤上皮查询细胞。颜色为原log2相对信号。完整原矩阵图像和色标保留，文字在最终尺寸下重新排布。

C：31个患者×程序表达谱间的Spearman相关，来自19名独立患者，每个表达谱至少20个细胞。横纵轴均为原表达谱；P01–P31为表达谱短索引，不是31名独立患者。C1–C7为原展示簇。保留原矩阵、展示排序、树及程序/研究色条；顶部色条为F1–F7，左侧研究色条沿原颜色：GSE156632蓝、GSE159115橙、GSE178481绿、GSE207493紫、GSE242299粉。完整身份对应已保存在服务器qc/PROFILE_LABEL_KEY.tsv；原附表S12仍为既有身份索引。

D：原24个调控子、七个程序、两种算法的患者平衡AUCell矩阵。横轴F1–F7，每个程序相邻两列依次为GENIE3和GRNBoost2；纵轴为调控子。各算法内按行Z标准化。21名独立患者形成37个患者×程序表达谱；F1–F7有效患者数分别为9、12、2、4、4、2、4。粗体标示ILF2、HOXD3和HNF4A，原行树和顺序保持。

E：三个TF与68个靶基因的预测调控网络，保留全部71个节点及71条有向连接。菱形为TF，灰蓝点为靶基因，金色为共享靶点；紫线为两算法共同支持、棕线仅GENIE3、蓝线仅GRNBoost2。全部基因标签改为水平排列，左右分列仅为排版，不代表新的分组或效应方向。连线为预测调控证据，不是因果验证。

F：F1使用量与四种蛋白推断活性的原密度网格。横轴F1使用量，纵轴metaVIPER NES；依次展示CASP9、CHMP1A、PRSS8和PXN，颜色为相对密度。原网格、等高线阈值和蛋白集合保留；统计单位为患者，原效应量、n、P、FDR仍见既有附表S13及S3。

面板对应：新A=候选①+原C；B保留；原D→C、E→D、F→E、G→F。当前独立图已实际拼好，未替换Phase332完整稿。未来并入全文时须同步修改旧图2C–G正文引用及附表S12标题，不能只换PDF图片。
'''
(o/'Fig2_新版图注与编号对应.txt').write_text(caption,encoding='utf-8')
placeholder={'status':'PLACEHOLDER_ONLY','candidate':'Candidate_3_cNMF_representative_genes','existing_preview':'phase333_Fig2A_candidates_20261019_v1/delivery/Candidate_3_cNMF_representative_genes.pdf','supplement_number':None,'supplement_figure_generated':False,'inserted_into_manuscript':False,'note':'仅占位，沿用既有候选预览；不新建空白补充图PDF，不插入全文。'}
(o/'补充图占位说明.txt').write_text('补充图占位：cNMF程序代表基因图（候选③）。\n当前只登记占位；正式编号待全文整合时确定，已有S1/S2保持。\n未新生成补充图、未插入全文；既有Phase333候选③预览仅作为后续来源。\n',encoding='utf-8')
(q/'SUPPLEMENT_PLACEHOLDER.json').write_text(json.dumps(placeholder,ensure_ascii=False,indent=2),encoding='utf-8')
manifest=[{'name':x.name,'bytes':x.stat().st_size,'sha256':hashlib.sha256(x.read_bytes()).hexdigest()} for x in sorted(o.iterdir()) if x.is_file()]
(q/'DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print('PASS: A4, 71 network labels, final render, immutable input hashes; supplement placeholder only')
