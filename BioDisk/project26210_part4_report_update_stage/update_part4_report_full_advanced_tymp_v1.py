#!/usr/bin/env python3
"""Add all post-advanced-analysis TJ-RCC and TYMP state results to the current report."""

from __future__ import annotations

import json
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from ooxml_report_utils_v1 import (
    A, REL, W, WP, any_paragraph, caption_figure, direct_paragraph, image_paragraph,
    insert_after, node_text, page_break, paragraph, png_size, q, set_paragraph_text,
    sha256, table, table_row,
)

ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
SOURCE = ROOT / "docs/client_reports/part4_spatial/project26210_第四部分空间转录组阶段性报告_客户样式综合图版_POSTN_APOE论文对标增强版_v2.docx"
STAGE = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_report_full_update_stage")
OUTPUT = STAGE / "project26210_part4_full_advanced_tymp_v2_candidate.docx"
METADATA = STAGE / "project26210_part4_full_advanced_tymp_v2_metadata.json"
EXPECTED_SOURCE_SHA256 = "f3b61f3334c8bdb5dbbd6cb5b674f5ea13acbb3bb7fa0ebc863c4d3b2cb86bf4"

ADV = ROOT / "external/ccRCC_spatial/TJ_RCC_figshare_24599295/advanced_full_v2"
COHORT = ADV / "cohort_synthesis_v1"
ADV_VERIFY = ADV / "verification_v2/verification_summary_v1.json"
ADV_MANUAL = ADV / "verification_v2/manual_figure_review_v2.tsv"
STATE = ROOT / "results/part4_spatial/tymp_spp1_state_v1"

ADV_IMAGES = [
    COHORT / "bayesspace_patient_overview_v1.png",
    COHORT / "bayesspace_program_effects_v1.png",
    COHORT / "spatial_program_top_recurrent_genes_v1.png",
    COHORT / "spatial_program_candidate_replication_v1.png",
    COHORT / "misty_patient_performance_v1.png",
    COHORT / "misty_patient_view_contributions_v1.png",
    COHORT / "commot_patient_transport_rank_v1.png",
    COHORT / "commot_patient_hypothesis_proxy_v1.png",
]
STATE_IMAGES = [
    STATE / "figures/cohort_primary_state_patient_effects_v1.png",
    STATE / "figures/cohort_route_concordance_v1.png",
    STATE / "figures/cohort_exploratory_mediation_heatmap_v1.png",
    STATE / "figures/cohort_interface_downstream_summary_v1.png",
]
IMAGES = ADV_IMAGES + STATE_IMAGES

TABLE_PATHS = {
    "bayes": COHORT / "bayesspace_patient_summary_v1.tsv",
    "bayes_effect": COHORT / "bayesspace_patient_program_effects_v1.tsv",
    "candidate": COHORT / "candidate_program_gene_replication_v1.tsv",
    "spatial_genes": COHORT / "spatial_program_gene_replication_v1.tsv",
    "misty": COHORT / "misty_patient_performance_summary_v1.tsv",
    "commot": COHORT / "commot_patient_hypothesis_proxy_summary_v1.tsv",
    "state": STATE / "tables/cohort_state_direction_summary_v1.tsv",
    "mediation": STATE / "tables/cohort_exploratory_mediation_summary_v1.tsv",
    "downstream": STATE / "tables/cohort_interface_downstream_summary_v1.tsv",
}


for path in [SOURCE, ADV_VERIFY, ADV_MANUAL, STATE / "verification_metadata_v1.json", STATE / "manual_figure_review_v1.tsv", *IMAGES, *TABLE_PATHS.values()]:
    if not path.is_file():
        raise FileNotFoundError(path)
if sha256(SOURCE) != EXPECTED_SOURCE_SHA256:
    raise RuntimeError("Current report hash changed; refusing to update an unreviewed source")
if OUTPUT.exists():
    raise RuntimeError(f"Refusing to overwrite staged candidate: {OUTPUT}")

adv_verify = json.loads(ADV_VERIFY.read_text())
state_verify = json.loads((STATE / "verification_metadata_v1.json").read_text())
assert adv_verify["status"] == "passed" and adv_verify["n_checks"] == 733 and adv_verify["n_failed"] == 0
assert adv_verify["n_samples"] == 12 and adv_verify["n_spots"] == 34609 and adv_verify["file_manifest_rows"] == 1716
assert state_verify == {"status": "passed", "checks": 211, "passed": 211, "failed": 0}
adv_manual = pd.read_csv(ADV_MANUAL, sep="\t")
assert len(adv_manual) == 12
assert int(adv_manual["n_figures"].sum()) == 56
assert set(adv_manual["status"].str.lower()) == {"passed"}
state_manual = pd.read_csv(STATE / "manual_figure_review_v1.tsv", sep="\t")
assert len(state_manual) == 16 and state_manual["overall"].eq("PASS").all()

bayes = pd.read_csv(TABLE_PATHS["bayes"], sep="\t")
effects = pd.read_csv(TABLE_PATHS["bayes_effect"], sep="\t")
candidate = pd.read_csv(TABLE_PATHS["candidate"], sep="\t")
misty = pd.read_csv(TABLE_PATHS["misty"], sep="\t")
commot = pd.read_csv(TABLE_PATHS["commot"], sep="\t")
state = pd.read_csv(TABLE_PATHS["state"], sep="\t")
mediation = pd.read_csv(TABLE_PATHS["mediation"], sep="\t")
downstream = pd.read_csv(TABLE_PATHS["downstream"], sep="\t")
assert len(bayes) == 12 and int(bayes.n_spots.sum()) == 34609 and set(bayes.selected_q) == {4, 5, 6}
assert len(candidate) == 37 and len(misty) == 120 and len(commot) == 19
assert len(state) == 48 and len(mediation) == 4 and len(downstream) == 10

tumor_effect = effects[effects.feature.eq("Tumor-epithelial")].tumor_minus_non_tumor_median
postn_effect = effects[effects.feature.eq("POSTN-ECM")].tumor_minus_non_tumor_median
assert len(tumor_effect) == len(postn_effect) == 12 and int((tumor_effect > 0).sum()) == 12 and int((postn_effect > 0).sum()) == 3
assert np.isclose(tumor_effect.median(), 0.39996737501423696) and np.isclose(wilcoxon(tumor_effect).pvalue, 0.00048828125)
assert np.isclose(postn_effect.median(), -0.22408572251509476) and np.isclose(wilcoxon(postn_effect).pvalue, 0.05224609375)

program_rows = [["候选程序", "基因数", "12/12双方法严格基因", "患者支持范围", "队列解读"]]
program_interpretation = {
    "POSTN-ECM": "7/7基因均稳定，但不证明POSTN直接通讯",
    "APOE-TAM": "APOE/C1Q稳定；TREM2/LPL为9/12",
    "TYMP-context": "SPP1/MIF/CD74/CXCR4/APP为12/12；TYMP为11/12",
    "CD8-effector": "无基因达到12/12双方法严格支持",
    "TGFb-response": "完整程序不稳定；SERPINE1为12/12",
    "Tumor-epithelial": "KRT8/KRT18为12/12；其余存在异质性",
}
for program in ["Tumor-epithelial", "POSTN-ECM", "APOE-TAM", "TYMP-context", "CD8-effector", "TGFb-response"]:
    frame = candidate[candidate.program.eq(program)]
    program_rows.append([
        program, str(len(frame)), f"{int((frame.n_both_strict == 12).sum())}/{len(frame)}",
        f"{int(frame.n_both_strict.min())}–{int(frame.n_both_strict.max())}/12", program_interpretation[program],
    ])
assert program_rows[2][2] == "7/7" and program_rows[3][2] == "4/6" and program_rows[4][2] == "5/6" and program_rows[5][2] == "0/6"

misty_targets = [
    ("program_POSTN_ECM", "POSTN-ECM"), ("program_APOE_TAM", "APOE-TAM"),
    ("program_CD8_effector", "CD8-effector"), ("cell_TYMP_High_Malignant", "TYMP-High Malignant"),
    ("cell_Endothelial", "Endothelial"),
]
misty_rows = [["目标", "中位ΔR²", "正向患者", "ΔR² BH FDR", "中位RMSE改善%", "解释"]]
for target, label in misty_targets:
    r2 = misty[(misty.target == target) & (misty.measure == "gain.R2")].iloc[0]
    rmse = misty[(misty.target == target) & (misty.measure == "gain.RMSE")].iloc[0]
    interpretation = "患者级稳定多尺度预测增益" if r2.wilcoxon_bh_fdr_within_measure < 0.05 else "未形成稳定多尺度增益"
    misty_rows.append([label, f"{r2['median']:.3f}", f"{int(r2.positive_n)}/12", f"{r2.wilcoxon_bh_fdr_within_measure:.4g}", f"{rmse['median']:.3f}", interpretation])
assert int((misty[(misty.measure == "gain.R2")].wilcoxon_bh_fdr_within_measure < 0.05).sum()) == 18

commot_selected = commot[((commot.axis_id.isin(["APOE_TREM2_TYROBP", "CCL2_CCR2", "POSTN_ITGAV_ITGB5", "SPP1_CD44", "MIF_CD74_CXCR4"])) & (commot.scale.isin(["secreted_200um", "secreted_500um"])))].copy()
commot_rows = [["候选轴", "尺度", "可估计患者", "正富集患者", "中位log2富集", "判断"]]
for row in commot_selected.itertuples(index=False):
    judgement = "稳定正向组成代理" if row.positive_enrichment_n == row.n_patients else ("方向趋势" if row.positive_enrichment_n >= 8 else "未获队列支持")
    commot_rows.append([row.axis_id, row.scale.replace("secreted_", "").replace("um", " µm"), str(row.n_patients), f"{row.positive_enrichment_n}/{row.n_patients}", f"{row.median_proxy_log2_enrichment:.3f}", judgement])

primary_state = state[(state.adjustment == "composition_qc") & (state.subset == "macrophage_enriched")]
state_rows = [["路线/指标", "尺度", "正向患者", "中位ρ", "方向P", "BH FDR", "解释"]]
state_selection = [
    ("reference_dependent", "SPP1_state_no_SPP1", "参考依赖：无SPP1极化评分"),
    ("reference_dependent", "SPP1_expression", "参考依赖：SPP1表达"),
    ("reference_light", "SPP1_state_no_SPP1", "参考轻依赖：无SPP1极化评分"),
    ("reference_light", "SPP1_vs_C1QC_contrast", "参考轻依赖：SPP1相对C1QC"),
]
for route, outcome, label in state_selection:
    for row in primary_state[(primary_state.route == route) & (primary_state.outcome == outcome)].sort_values("grid").itertuples(index=False):
        state_rows.append([label, row.grid, f"{row.positive_patients}/12", f"{row.median_rho:.3f}", f"{row.sign_test_p:.4g}", f"{row.sign_test_fdr_bh:.4g}", "趋势性正向" if row.median_rho > 0 else "SPP1特异性不支持"])

mechanism_rows = [["分析层", "尺度", "方向一致性", "中位效应", "BH FDR", "结论"]]
for row in mediation.sort_values(["route", "grid"]).itertuples(index=False):
    mechanism_rows.append([f"MTA间接效应（{row.route}）", row.grid, f"{row.indirect_a_times_b_positive_patients}/12正向", f"{row.indirect_a_times_b_median:.3f}", f"{row.indirect_a_times_b_sign_fdr_bh:.3f}", "未支持转录代理中介"])
for row in downstream[downstream.outcome.eq("Invasion-EMT context")].sort_values("grid").itertuples(index=False):
    mechanism_rows.append(["TYMP–SPP1界面→侵袭/EMT", row.grid, f"{row.positive_patients}/12正向", f"ρ={row.median_rho:.3f}", f"{row.sign_test_fdr_bh:.3f}", "探索性表型趋势"])

advanced_summary_rows = [
    ["方法", "12例主要结果", "患者级证据", "证据边界"],
    ["BayesSpace", "q=4–6；肿瘤样spot中位37.4%；边界spot中位9.4%", "Tumor-epithelial 12/12升高，P=0.000488", "计算域，不是病理边界"],
    ["SPARK-X/Hotspot", "345基因12/12双方法严格支持", "POSTN-ECM 7/7；TYMP-context 5/6", "空间变异不等于通讯或因果"],
    ["MISTy", "20目标中18个ΔR²经患者级BH FDR<0.05", "POSTN-ECM、APOE-TAM和TYMP-high均12/12正增益", "预测重要性无方向且非因果"],
    ["靶向COMMOT", "9轴12/12可估计；CCL2–CCR2 10/12；ITGB3两轴0/12", "APOE–TREM2/TYROBP两尺度12/12正代理富集", "组成加权最优传输，不证明直接结合"],
]

with ZipFile(SOURCE) as source_zip:
    document_root = ET.fromstring(source_zip.read("word/document.xml"))
    rels_root = ET.fromstring(source_zip.read("word/_rels/document.xml.rels"))
    source_names = source_zip.namelist()
body = document_root.find(q(W, "body"))
if body is None:
    raise RuntimeError("DOCX body missing")
all_source_text = "\n".join(node_text(p) for p in document_root.iter(q(W, "p")))
if "4.12 【高级分析补丁】" not in all_source_text or "4.13 【TJ-RCC 12例全套高级分析】" in all_source_text:
    raise RuntimeError("Unexpected report section state")

set_paragraph_text(direct_paragraph(body, "正式报告｜"), "正式报告｜GLMF1/GLMF2与TJ-RCC 12例全套空间分析｜高级空间分析与TYMP–MTA–SPP1机制验证整合版", bold=True)
set_paragraph_text(any_paragraph(body, "2026年8月2日"), "2026年8月3日（12例高级分析与TYMP状态更新）", bold=True)
set_paragraph_text(
    direct_paragraph(body, "【增强版新增说明】"),
    "【本次整合新增说明】保留原4.12节GLMF1/GLMF2高级补丁，新增4.13节TJ-RCC 12例BayesSpace、SPARK-X/Hotspot、MISTy与靶向COMMOT完整队列结果，以及4.14节TYMP–MTA-context–SPP1巨噬细胞状态分析。湿实验结论作为因果锚点；空间结果用于人组织相关性、空间组织和患者异质性验证。",
    bold=True, shading="FFF2CC",
)
insert_after(body, direct_paragraph(body, "高级分析补丁显示："), [
    paragraph("12例全套高级分析进一步将两切片发现扩展到34,609个spots和12名独立患者：BayesSpace肿瘤样程序12/12升高；SPARK-X/Hotspot得到345个12/12双方法严格空间基因；MISTy 20个锁定目标中18个患者级ΔR²经BH校正通过；靶向COMMOT最稳定的是APOE–TREM2/TYROBP组成代理，而POSTN–APOE和SPP1–CD44未获队列支持。"),
    paragraph("TYMP状态分析不再只比较SPP1+ TAM数量，而是在巨噬细胞富集区域校正巨噬细胞量和技术协变量后检验极化状态。参考轻依赖的无SPP1极化评分在两个尺度均10/12正向，但全局BH FDR=0.154；SPP1相对C1QC对比方向为负，探索性MTA间接效应接近零。故空间数据提供趋势性人组织支持，但不能独立验证5′-MTA中介或SPP1特异性。"),
])

existing_rids = [int(x.get("Id")[3:]) for x in rels_root.findall(q(REL, "Relationship")) if (x.get("Id") or "").startswith("rId") and x.get("Id")[3:].isdigit()]
next_rid = max(existing_rids, default=0) + 1
docpr_ids = [int(n.get("id")) for n in document_root.iter(q(WP, "docPr")) if (n.get("id") or "").isdigit()]
next_docpr = max(docpr_ids, default=0) + 1
image_specs = []
for index, image in enumerate(IMAGES, start=1):
    rid = f"rId{next_rid}"; next_rid += 1
    media_name = f"full_advanced_tymp_v1_{index:02d}_{image.name}"
    ET.SubElement(rels_root, q(REL, "Relationship"), {"Id": rid, "Type": "http://schemas.openxmlformats.org/officeDocument/2006/relationships/image", "Target": f"media/{media_name}"})
    width, height = png_size(image)
    image_specs.append((image, rid, media_name, width, height, next_docpr)); next_docpr += 1

section = [
    page_break(),
    paragraph("4.13 【TJ-RCC 12例全套高级分析】BayesSpace、SPARK-X/Hotspot、MISTy与靶向COMMOT", style="Heading2"),
    paragraph("在GLMF1/GLMF2两切片发现层之后，本节以12名独立TJ-RCC患者为生物学重复，使用与发现层一致的四类高级方法完成队列验证。12例共34,609个spots、每例11,702个共享基因；所有方法均12/12完成，独立验收733/733项通过。患者级汇总优先于spot级显著性。"),
    table(advanced_summary_rows, 14),
    paragraph("4.13.1 BayesSpace计算域及程序分布", style="Heading3"),
    paragraph("12例选择q=4–6；计算肿瘤样spot比例患者中位数为0.374，计算边界spot比例中位数为0.094。Tumor-epithelial程序在计算肿瘤样域相对非肿瘤样域12/12升高，中位差0.400，患者级Wilcoxon双侧P=0.000488。POSTN-ECM仅3/12升高，中位差−0.224，P=0.0522；APOE-TAM和CD8-effector没有队列一致方向。所有域及边界均为表达驱动计算结果，未经过病理医师边界标注。"),
    *caption_figure(image_specs[0], "图4.38  TJ-RCC 12例BayesSpace患者级计算域概览。q为每例独立选择的计算域数；肿瘤样域依据Tumor-epithelial程序与恶性细胞丰度注释，不能替代病理边界。"),
    *caption_figure(image_specs[1], "图4.39  预设空间程序在BayesSpace计算肿瘤样域与非肿瘤样域之间的患者级差异。"),
    page_break(),
    paragraph("4.13.2 SPARK-X/Hotspot无偏空间程序的患者级重复", style="Heading3"),
    paragraph("全基因分析得到345个基因在12/12患者中同时满足SPARK-X严格阈值与Hotspot正空间自相关严格阈值。锁定候选中，POSTN-ECM的7/7基因均为12/12双方法严格支持；APOE-TAM为4/6，TYMP-context为5/6，CD8-effector为0/6。TYMP-context中的SPP1、MIF、CD74、CXCR4和APP为12/12，TYMP为11/12。该结果验证的是基因空间异质性，不意味着这些基因属于同一模块或构成因果通讯链。"),
    table(program_rows, 14),
    *caption_figure(image_specs[2], "图4.40  TJ-RCC队列中SPARK-X/Hotspot双方法最稳定的重复空间基因。"),
    *caption_figure(image_specs[3], "图4.41  37个锁定候选程序基因的患者级空间重复次数。"),
    page_break(),
    paragraph("4.13.3 MISTy多尺度生态预测", style="Heading3"),
    paragraph("固定20个目标、intra/juxta/para三视图和相同模型参数后，18/20目标的患者级ΔR²经目标内12患者Wilcoxon检验并在同一measure内BH校正后FDR<0.05。POSTN-ECM中位ΔR²=1.587且12/12为正，APOE-TAM为0.269且12/12为正，TYMP-High Malignant为1.508且12/12为正；CD8-effector中位仅0.046、7/12为正、FDR=0.519。MISTy表示多尺度上下文的预测贡献，不提供正负作用方向、直接配体–受体关系或因果解释。"),
    table(misty_rows, 14),
    *caption_figure(image_specs[4], "图4.42  TJ-RCC 12例固定目标的MISTy患者级预测增益。"),
    *caption_figure(image_specs[5], "图4.43  intra、juxta与para视图对20个目标的患者级贡献构成。"),
    page_break(),
    paragraph("4.13.4 靶向COMMOT及可估计性边界", style="Heading3"),
    paragraph("12条预设轴中9条在12/12患者可估计，CCL2–CCR2在10/12可估计；依赖ITGB3的POSTN和SPP1两轴因12例ITGB3均为0而0/12可估计。APOE–TREM2/TYROBP在APOE+C1QC+ TAM自代理中200/500 µm均12/12正富集，中位log2富集0.357/0.274。CCL2–CCR2为8/10和9/10正向。POSTN–ITGAV/ITGB5的SMC/Pericyte→APOE+C1QC+ TAM代理仅4/12正向且中位为负；SPP1–CD44亦未获队列支持。COMMOT结果为混合spot组成加权的集体最优传输代理，不证明直接结合、细胞分辨通讯或因果方向。"),
    table(commot_rows, 14),
    *caption_figure(image_specs[6], "图4.44  TJ-RCC靶向COMMOT候选轴的患者级传输质量排名。"),
    *caption_figure(image_specs[7], "图4.45  预设来源–靶代理的患者级COMMOT组成富集重复。", width_emu=5_250_000),
    paragraph("4.13.5 完成性、QA与患者图", style="Heading3"),
    paragraph("BayesSpace、SPARK-X/Hotspot、MISTy和COMMOT均12/12完成。独立验证733/733项通过，登记1,716个文件；8张队列图和BayesSpace、空间程序、MISTy、COMMOT四类各12张患者关键图共56张正式PNG均已人工打开检查。正文展示8张队列图，逐病例图保留在正式结果目录作为补充审计材料。"),
    page_break(),
    paragraph("4.14 【TYMP机制状态分析】TYMP–MTA-context–SPP1巨噬细胞状态与下游表型", style="Heading2"),
    paragraph("客户湿实验已提出“TYMP介导5′-MTA蓄积，通过重塑SPP1⁺巨噬细胞生态位促进ccRCC恶性演进”。本节不重复建立湿实验因果，而检验该机制在人肿瘤组织中的空间相关性。主问题由‘SPP1+ TAM数量是否升高’改为‘在巨噬细胞数量相当时，邻近TYMP状态是否对应更强SPP1样极化’。"),
    paragraph("主评分以GPNMB、LGALS3、FABP5、CD9和LPL构成，并明确排除SPP1本身；参考依赖路线使用cell2location的TYMP-high恶性状态与巨噬细胞总量，参考轻依赖路线使用TYMP表达残差和泛巨噬细胞转录评分。两条路线均校正UMI、检测基因数及肿瘤/免疫组成，并在巨噬细胞富集spots中做4×4、5×5分块与9,999次置换。"),
    table(state_rows, 14),
    *caption_figure(image_specs[8], "图4.46  12例患者中校正后SPP1状态与邻域TYMP-high恶性状态的空间分块效应。"),
    *caption_figure(image_specs[9], "图4.47  参考依赖与参考轻依赖TYMP–SPP1状态效应的一致性。"),
    paragraph("结果解读：参考依赖主分析的无SPP1极化评分在4×4为9/12正向、中位ρ=0.179、FDR=0.305；5×5为7/12正向、中位ρ=0.214、FDR=0.864。参考轻依赖路线两尺度均10/12正向，中位ρ=0.278/0.230，名义P=0.0386但全局BH FDR=0.154，属于趋势性人组织支持。相反，SPP1相对C1QC对比在参考轻依赖巨噬细胞富集分析中4×4为11/12负向、5×5为12/12负向，5×5 FDR=0.00781，提示更接近广义巨噬细胞状态重塑，而非稳定的SPP1特异性替代。"),
    page_break(),
    paragraph("4.14.1 MTA转录代理中介与下游表型", style="Heading3"),
    paragraph("MTA-context由MAT2A、MAT2B、AMD1、SRM和SMS的转录评分构成，仅代表多胺/MTA相关转录背景。探索性路径分别估计TYMP→MTA-context、MTA-context→SPP1状态及间接效应a×b。两条路线、两个尺度的间接效应仅4–5/12为正，中位值接近0，BH FDR=0.775–0.826，未支持转录代理中介。"),
    table(mechanism_rows, 14),
    *caption_figure(image_specs[10], "图4.48  TYMP→MTA-context→SPP1状态的探索性患者级间接效应。数值为标准化路径系数乘积，不是代谢物实测或因果中介证明。"),
    *caption_figure(image_specs[11], "图4.49  TYMP–SPP1共激活界面与下游免疫及恶性程序的患者级方向汇总。"),
    paragraph("TYMP–SPP1共激活界面与侵袭/EMT背景在4×4和5×5均10/12正向，中位ρ=0.269/0.304，名义P=0.0386，BH FDR=0.193；免疫抑制、T细胞功能障碍、细胞毒性和缺氧/血管生成均未形成校正后的队列一致性。因此侵袭/EMT只能作为值得实验对照的探索性表型趋势，不能替代纵向恶性进展、预后或干预证据。"),
    paragraph("4.14.2 QA与机制证据边界", style="Heading3"),
    paragraph("12例34,609个spots均完成；独立验证211/211项通过，4张队列图及12张患者图共16/16人工图检通过。Visium spot混合多种细胞；MTA-context不是5′-MTA浓度；探索性中介不是生化中介；患者方向检验才是生物学重复层。客户湿实验应作为因果锚点，本节只提供人组织空间相关性、患者异质性及可转化验证。", bold=True, shading="FFF2CC"),
]

heading5 = direct_paragraph(body, "5. 主要结论", style="Heading1")
insert_index = list(body).index(heading5)
for node in section:
    body.insert(insert_index, node); insert_index += 1

insert_after(body, direct_paragraph(body, "【高级分析补丁结论】"), [
    paragraph("【12例全套高级分析结论】12例患者级结果确认了可重复的空间异质性和多尺度可预测结构：Tumor-epithelial程序在计算肿瘤样域12/12升高，345个基因获12/12 SPARK-X/Hotspot双方法严格支持，18/20 MISTy目标的患者级ΔR²通过BH校正。靶向COMMOT最稳定的是APOE–TREM2/TYROBP组成代理；POSTN–APOE和SPP1–CD44未获队列支持。"),
    paragraph("【TYMP机制状态结论】校正巨噬细胞数量后，参考轻依赖的无SPP1极化评分在两个尺度均10/12正向，但全局FDR=0.154；参考依赖结果较弱且尺度敏感。SPP1相对C1QC状态主要为负，MTA转录代理中介未获支持。结合客户湿实验，可表述为TYMP相关恶性区域与巨噬细胞状态重塑在人ccRCC中具有趋势性空间相关，但空间数据不能独立证明5′-MTA中介、SPP1特异性或促恶性因果。"),
])
insert_after(body, direct_paragraph(body, "【高级分析补丁限制】"), [
    paragraph("【12例高级分析限制】BayesSpace域不是病理边界；SPARK-X/Hotspot高比例空间信号可受大尺度组织结构影响；MISTy重要性无效应方向；COMMOT含对角运输和混合spot组成代理。患者数为12，病例内spots和空间块不能当作独立患者。"),
    paragraph("【TYMP状态分析限制】SPP1状态仍由混合spot转录评分近似；参考依赖路线沿用既有cell2location参考，参考轻依赖路线也不能获得单细胞分辨率。MTA-context没有测量5′-MTA、MTAP功能或代谢通量；缺少纵向进展、预后、干预和治疗应答标签。"),
])

status_tables = [node for node in body.iter(q(W, "tbl")) if "TJ-RCC 12例cell2location及全套下游分析" in node_text(node)]
if len(status_tables) != 1:
    raise RuntimeError(f"Expected one chapter-5 status table, found {len(status_tables)}")
status_tables[0].append(table_row(["TJ-RCC 12例四类全套高级分析", "已完成；733/733独立检查与56/56图检通过", "已纳入第4.13节和第5章"], font_half_points=16))
status_tables[0].append(table_row(["TYMP–MTA-context–SPP1状态分析", "已完成；211/211独立检查与16/16图检通过", "已纳入第4.14节和第5章"], font_half_points=16))

set_paragraph_text(
    any_paragraph(body, "最终总述"),
    "最终总述　第四部分现已形成三层证据：GLMF1/GLMF2用于空间发现和方法深化，TJ-RCC 12名独立患者用于全套高级方法重复，客户湿实验作为TYMP–5′-MTA–SPP1机制的因果锚点。12例分析稳定支持肿瘤空间程序、POSTN-ECM基因空间异质性、广泛MISTy多尺度预测增益和APOE–TREM2/TYROBP组成代理；TYMP相关巨噬细胞状态获得趋势性空间支持，但SPP1特异性和MTA转录代理中介未获独立验证。故报告可以支持该湿实验机制在人ccRCC组织中的空间相关性与患者异质性，不能把空间共现、预测或最优传输改写为已独立证明的代谢中介、直接通讯或促恶性因果。",
    bold=True,
)

document_bytes = ET.tostring(document_root, encoding="utf-8", xml_declaration=True)
rels_bytes = ET.tostring(rels_root, encoding="utf-8", xml_declaration=True)
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
with ZipFile(SOURCE) as source_zip, ZipFile(OUTPUT, "w", ZIP_DEFLATED, allowZip64=True) as output_zip:
    for name in source_names:
        if name == "word/document.xml":
            output_zip.writestr(name, document_bytes)
        elif name == "word/_rels/document.xml.rels":
            output_zip.writestr(name, rels_bytes)
        else:
            output_zip.writestr(source_zip.getinfo(name), source_zip.read(name))
    for image, _, media_name, _, _, _ in image_specs:
        output_zip.writestr(f"word/media/{media_name}", image.read_bytes())

with ZipFile(OUTPUT) as archive:
    if archive.testzip() is not None:
        raise RuntimeError("Candidate DOCX failed ZIP CRC")
    final_root = ET.fromstring(archive.read("word/document.xml"))
    final_text = "\n".join(node_text(p) for p in final_root.iter(q(W, "p")))
    for marker in ["4.13 【TJ-RCC 12例全套高级分析】", "4.14 【TYMP机制状态分析】", "图4.38", "图4.49", "733/733", "211/211", "【TYMP机制状态结论】", "【TYMP状态分析限制】"]:
        if marker not in final_text:
            raise RuntimeError(f"Required marker missing: {marker}")
    for _, _, media_name, _, _, _ in image_specs:
        if f"word/media/{media_name}" not in archive.namelist():
            raise RuntimeError(f"Embedded image missing: {media_name}")

metadata = {
    "status": "staged_and_structurally_verified",
    "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    "source": str(SOURCE), "source_size": SOURCE.stat().st_size, "source_sha256": sha256(SOURCE),
    "candidate": str(OUTPUT), "candidate_size": OUTPUT.stat().st_size, "candidate_sha256": sha256(OUTPUT),
    "added_sections": ["4.13", "4.14"], "added_figures": [f"4.{number}" for number in range(38, 50)],
    "added_tables": 6, "advanced_verification": "733/733", "state_verification": "211/211",
    "embedded_images": [{"source": str(path), "embedded_name": media_name, "size": path.stat().st_size, "sha256": sha256(path)} for path, _, media_name, _, _, _ in image_specs],
    "source_tables": [{"path": str(path), "size": path.stat().st_size, "sha256": sha256(path)} for path in TABLE_PATHS.values()],
    "report_date": "2026-08-03",
    "inference_boundary": "Wet experiments remain the causal anchor; spatial analyses support human-tissue relevance and heterogeneity only.",
}
METADATA.write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps(metadata, ensure_ascii=False, indent=2))
