from pathlib import Path
import json,hashlib,shutil,datetime
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase194_mainfigure_assembly_20260917_v1'
txt="""# 2.1—2.12主图外部拼版 v1

本次将既有图件中的代表性面板重排为3张独立主图。正文PDF未修改；分析数据、阈值、颜色数值、原始轨迹和统计结果均未重算。主图3使用最新2.11等高线修订图。各图为420×470 mm矢量拼版母版，便于独立放大审阅；PDF保留原图矢量要素及原始嵌入图像，另有SVG和600 dpi PNG。JPG用于快速预览。若用于期刊固定版芯，需要按该版芯再作字号与布局调整。

## 主图1｜TYMP相关恶性上皮转录程序及其调控背景

A，整合单细胞图谱及主要细胞组成，258,120个细胞、75个文库、40名个体。B，原始代表性患者GSE156632 RCC2的染色体CNV图，保留参考细胞及肿瘤取材上皮细胞区域。C，恶性上皮cNMF程序、F1使用度及TYMP RNA在同一嵌入中的分布。D，程序基因贡献与按F1使用度排序的表达热图，保留F1从高到低的原顺序。E，两种算法的患者平衡调控子活性。F，四个候选蛋白推断活性与F1使用度的既有密度图，同时保留原患者级统计及PXN的非显著结果。

此主图承接2.1—2.5。F1主导细胞在原分析中高度集中于GSE207493（99.66%），原队列组成结果仍保留在原Figure3及报告材料中；面板C/D的相似表达分布不等于跨队列普遍验证。E中的C1—C4是原调控子显示聚类，与后续巨噬细胞C4没有身份对应关系。蛋白活性为推断值，F图保留原始加权细胞密度和独立患者统计口径。

## 主图2｜代谢线索与SPP1富集巨噬状态的发现及组织联系

A，恶性上皮细胞四个代谢通路的原AUCell图。B，从蛋氨酸/SAM与多胺反应衔接MTA的生化示意。C/D，保留Monocle3与scTour两种方法的原恶性上皮表达动态热图，并列呈现方法间差别。E，冻结成熟巨噬细胞图谱及C4定位。F，C4和其他候选群的多基因表达热图。G，4个队列255例实测MTA与组织表达图谱；先展示独立发现的状态E/F，再展示组织联系G。

此主图承接2.6—2.10。C/D使用恶性上皮轨迹，不能与主图3的髓系轨迹拼成同一条细胞转化路径。B为已知生化反应，不表示TYMP直接催化MTA。G保留全部患者及按MTA排序的原图，MTA—组织SPP1为名义显著候选联系（原P=0.02926、FDR=0.09755），不能据此认定MTA诱导C4。

## 主图3｜C4功能身份与髓系状态排序

A，成熟巨噬候选群的功能程序热图，保留抗原呈递、吞噬体等并非统一增强的程序。B，2.11最新分层等高线版，四小面板依次为SPP1、胆固醇代谢、溶酶体及抗原呈递；灰线为等表达线，紫色虚线为C4轮廓，实际数值色标保留。C，5研究原始分区轨迹，保留无有效伪时间的灰色区域。D，患者平衡的主分析及敏感性表达动态，两种结果并列。

此主图承接2.11—2.12。C4不是SPP1的唯一来源；嵌入图不代表组织空间。主分析获得有效伪时间的C4为442/2674，敏感性分析为1743/2674；D分别为14与24名患者。主分析SPP1动态原FDR=0.07177，敏感性分析原FDR=0.03069。拼版不改变原证据层级。

## 取舍与原图保留

为保持3张主图的可读性，本次选择代表性面板，不把15张独立/复合图整页缩成缩略图。未纳入的面板仍保存在原图与报告中，包括全细胞标志点图及多基因分布、CNV分类与亚克隆小提琴图、cNMF患者结构及队列构成、调控子排名和完整网络、四蛋白MR效应密度图、85通路全图和scFEA代理、恶性上皮原始轨迹总览。重要证据边界和阴性结果仍在原报告、原图及上述说明中保留。

新增A—G等大写字母表示主图面板；所含原图的小面板字母在当前试拼版中保留，用于复核。详细裁切坐标和源PDF校验值在服务器归档。
"""
(d/'MAIN_FIGURE_GUIDE_CN.md').write_text(txt,encoding='utf8')
rows=json.loads((d/'tables/PANEL_PROVENANCE.json').read_text())
text='# 主图面板来源表\n\n| 主图 | 面板 | 内容 | 源图 |\n|---|---|---|---|\n'
for z in rows:text+=f"| {z['main_figure'].split('_')[2]} | {z['panel']} | {z['label']} | {Path(z['source']).name} |\n"
(d/'PANEL_INDEX_CN.md').write_text(text,encoding='utf8')
inputs=set(z['source'] for z in rows)
report=r/'phase191_discovery_order_report_20260917_v1/report/26210_TYMP_MTA_SPP1_macrophage_integrated_report_v3.pdf'
inputs.add(str(report))
(d/'tables/SOURCE_SHA256.json').write_text(json.dumps({s:hashlib.sha256(Path(s).read_bytes()).hexdigest() for s in sorted(inputs)},indent=2))
print('Guide and provenance saved')
