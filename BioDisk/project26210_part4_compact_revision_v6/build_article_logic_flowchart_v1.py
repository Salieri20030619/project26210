#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


OUTDIR = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_compact_revision_v6")
PNG = OUTDIR / "article_logic_flowchart_v2.png"
PDF = OUTDIR / "article_logic_flowchart_v2.pdf"
META = OUTDIR / "article_logic_flowchart_metadata_v2.json"
FONT_PATH = Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")

INK = "#24344D"
MUTED = "#66758A"
BLUE = "#2E75B6"
TEAL = "#169C78"
PURPLE = "#7656A8"
ORANGE = "#D98B18"
RED = "#C64F5A"
LIGHT_BLUE = "#EAF3FA"
LIGHT_TEAL = "#E9F6F1"
LIGHT_PURPLE = "#F1EDF7"
LIGHT_ORANGE = "#FFF3DE"
LIGHT_GRAY = "#F3F5F8"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def box(ax, x, y, w, h, title, body, face, edge, title_size=13, body_size=10.2):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        facecolor=face, edgecolor=edge, linewidth=1.7,
    )
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h * .68, title, ha="center", va="center",
            fontsize=title_size, fontweight="bold", color=edge)
    ax.text(x + w / 2, y + h * .30, body, ha="center", va="center",
            fontsize=body_size, color=INK, linespacing=1.28)


def arrow(ax, x1, y1, x2, y2, color="#8A96A6", lw=2.0):
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=14,
        linewidth=lw, color=color, connectionstyle="arc3,rad=0.0",
    ))


def main() -> None:
    for path in (PNG, PDF, META):
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    font_manager.fontManager.addfont(str(FONT_PATH))
    font_name = font_manager.FontProperties(fname=str(FONT_PATH)).get_name()
    plt.rcParams.update({
        "font.family":font_name,
        "axes.unicode_minus":False,
        "pdf.fonttype":42,
        "savefig.facecolor":"white",
    })
    fig, ax = plt.subplots(figsize=(20, 8.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.subplots_adjust(left=.02, right=.98, bottom=.04, top=.90)
    fig.suptitle("报告逻辑主线：从空间发现到患者验证与证据收敛",
                 fontsize=22, fontweight="bold", color=INK, y=.97)
    ax.text(.5, .94, "核心问题：TYMP相关恶性状态是否对应可重复的髓系/血管空间生态，并能否在独立患者中形成可信的人组织证据？",
            ha="center", va="center", fontsize=12.8, color=MUTED)

    # Column 1: data foundation.
    ax.text(.105, .875, "① 数据基础与统计单位", ha="center", fontsize=14.5,
            fontweight="bold", color=BLUE)
    box(ax, .025, .64, .16, .17, "GLMF1 / GLMF2", "2张发现切片\nspot与切片内空间结构", LIGHT_BLUE, BLUE)
    box(ax, .025, .405, .16, .17, "TJ-RCC", "12名独立患者\n患者作为生物学重复", LIGHT_BLUE, BLUE)
    box(ax, .025, .17, .16, .17, "单细胞参考", "53,794个髓系细胞\n63个合格样本进入DPT检验", LIGHT_BLUE, BLUE)

    # Column 2: analysis strategy.
    ax.text(.305, .875, "② 分析与验证策略", ha="center", fontsize=14.5,
            fontweight="bold", color=PURPLE)
    box(ax, .225, .64, .16, .17, "空间定位", "统一参考映射、GSVA\n同spot与邻域结构", LIGHT_PURPLE, PURPLE)
    box(ax, .225, .405, .16, .17, "稳健性审计", "4×4 / 5×5分块置换\n组成/QC校正与尺度复核", LIGHT_PURPLE, PURPLE)
    box(ax, .225, .17, .16, .17, "跨患者收敛", "患者效应、方向检验\nBH FDR与路线敏感性", LIGHT_PURPLE, PURPLE)
    for y in (.725, .49, .255):
        arrow(ax, .185, y, .225, y)

    # Column 3: three parallel biological themes.
    ax.text(.52, .875, "③ 三条生物学主线", ha="center", fontsize=14.5,
            fontweight="bold", color=TEAL)
    box(ax, .425, .64, .19, .17, "主线A｜巨噬细胞状态", "TYMP–MTA-context–SPP1\n空间相容性、状态偏移与中介审计", LIGHT_TEAL, TEAL)
    box(ax, .425, .405, .19, .17, "主线B｜血管生态位", "SPP1–vascular niche\n与T-cell balance的患者重复关系", LIGHT_TEAL, TEAL)
    box(ax, .425, .17, .19, .17, "主线C｜髓系连续谱", "横断面DPT状态排序\n比较C1QC增强与SPP1相对优势", LIGHT_TEAL, TEAL)
    for y in (.725, .49, .255):
        arrow(ax, .385, y, .425, y)

    # Column 4: evidence grading.
    ax.text(.725, .875, "④ 证据分级与边界", ha="center", fontsize=14.5,
            fontweight="bold", color=ORANGE)
    box(ax, .645, .64, .16, .17, "获得支持", "患者重复的空间关系\nC1QC偏移的横断面连续谱", LIGHT_ORANGE, ORANGE)
    box(ax, .645, .405, .16, .17, "条件性支持", "部分关系受路线、尺度\n或细胞组成影响", LIGHT_ORANGE, ORANGE)
    box(ax, .645, .17, .16, .17, "未被证明", "5′-MTA生化中介、直接通讯\n谱系转化与治疗耐药", "#FCECEE", RED)
    for y in (.725, .49, .255):
        arrow(ax, .605, y, .645, y)

    # Column 5: final conclusion.
    ax.text(.902, .875, "⑤ 结论表达", ha="center", fontsize=14.5,
            fontweight="bold", color=INK)
    box(ax, .835, .52, .14, .26, "可信的人组织框架",
        "TYMP相关恶性区域\n—髓系/血管生态位—\nT细胞平衡\n\n强调患者级关联，\n不外推为闭合因果链",
        LIGHT_GRAY, INK, title_size=13.5, body_size=10.4)
    arrow(ax, .805, .605, .835, .65, color=INK, lw=2.4)
    ax.text(.902, .41, "正文按“发现 → 验证 → 收敛 → 分级”展开",
            ha="center", fontsize=10.8, fontweight="bold", color=MUTED)

    ax.text(.5, .055,
            "阅读提示：图4.1聚焦巨噬细胞状态，图4.2聚焦血管生态位与T细胞平衡，图4.3聚焦髓系横断面连续谱。",
            ha="center", fontsize=11.3, color=INK,
            bbox={"facecolor":"#F7F9FC", "edgecolor":"#D6DDE7", "boxstyle":"round,pad=.45"})

    fig.savefig(PNG, dpi=360, bbox_inches="tight", facecolor="white")
    fig.savefig(PDF, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    metadata = {
        "status":"generated_pending_report_integration_qc",
        "generated_at":datetime.now().astimezone().isoformat(timespec="seconds"),
        "png":str(PNG), "png_size":PNG.stat().st_size, "png_sha256":sha256(PNG),
        "pdf":str(PDF), "pdf_size":PDF.stat().st_size, "pdf_sha256":sha256(PDF),
        "logic":"data foundation -> validation strategy -> three biological themes -> evidence grading -> bounded conclusion",
    }
    META.write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
