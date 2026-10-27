from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from docx import Document
from docx.text.paragraph import Paragraph


ROOT = Path(__file__).resolve().parents[1]
PHASE = Path(__file__).resolve().parent
SOURCE = (
    ROOT
    / "phase46_part1_clone_tree_and_fig14_revision_local_v1"
    / "PART1_EXTERNAL_SINGLE_CELL_CLIENT_REPORT_WITH_FIGURES_REVIEW_v10.docx"
)
OUTPUT = PHASE / "PART1_EXTERNAL_SINGLE_CELL_CLIENT_REPORT_WITH_FIGURES_REVIEW_v11.docx"
QC = PHASE / "qc_local_v1"


def find_paragraph(document: Document, exact_text: str) -> Paragraph:
    matches = [paragraph for paragraph in document.paragraphs if paragraph.text == exact_text]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one paragraph matching {exact_text!r}; found {len(matches)}")
    return matches[0]


def find_paragraph_starting(document: Document, prefix: str) -> Paragraph:
    matches = [paragraph for paragraph in document.paragraphs if paragraph.text.startswith(prefix)]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one paragraph starting with {prefix!r}; found {len(matches)}")
    return matches[0]


def replace_text(paragraph: Paragraph, new_text: str) -> None:
    if not paragraph.runs:
        paragraph.add_run(new_text)
        return
    paragraph.runs[0].text = new_text
    for run in paragraph.runs[1:]:
        run.text = ""


def insert_before(anchor: Paragraph, text: str, style: str = "Normal") -> Paragraph:
    paragraph = anchor.insert_paragraph_before(text, style=style)
    paragraph.paragraph_format.keep_with_next = style.startswith("Heading")
    return paragraph


def add_labeled_paragraph_before(anchor: Paragraph, label: str, body: str) -> Paragraph:
    paragraph = anchor.insert_paragraph_before(style="Normal")
    lead = paragraph.add_run(label)
    lead.bold = True
    paragraph.add_run(body)
    return paragraph


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {OUTPUT}")
    if not SOURCE.exists():
        raise FileNotFoundError(SOURCE)

    QC.mkdir(parents=True, exist_ok=True)
    document = Document(SOURCE)

    version_note = find_paragraph_starting(document, "版本说明  客户图文审阅稿")
    replace_text(
        version_note,
        "版本说明  客户图文审阅稿：保留全部26项分析与15组正文结果图，并纳入新增患者三层恶性口径复核及参考文献桥接筛查；正文继续采用患者层统计和分级判定。",
    )

    summary_anchor = find_paragraph(document, "1  数据范围与分析设计")
    add_labeled_paragraph_before(
        summary_anchor,
        "增量复核  ",
        "新增患者复核不改变原有CNV阈值，而是在严格保守CNV、任一CNV方法阳性的扩展口径和无阳性CNV要求的暂定肿瘤样口径下并列评估。GSE178481仅在暂定口径支持TYMP—SPP1正相关；参考桥接筛查支持SPP1髓系免疫抑制状态及局部POSTN⁺ CAF—APOE⁺ TAM连接，但尚不支持完整的F1/TYMP→CAF→TAM→CD8因果链。",
    )

    statistics = find_paragraph_starting(document, "患者是独立统计单位。")
    replace_text(
        statistics,
        "患者是独立统计单位。细胞用于形成患者内效应，不将同一患者的细胞当作独立重复。患者数不少于10时，同时报告效应和BH校正后的FDR；患者数3–9时，采用双侧精确置换检验，P<0.05可报告为小样本显著结果，否则只描述方向。新增患者复核沿用冻结阈值，并按严格保守CNV、任一CNV方法阳性的扩展口径和暂定肿瘤样口径分层报告。",
    )

    section7 = find_paragraph(document, "7  综合解释")
    insert_before(section7, "C06  新增患者三层复核与参考桥接筛查", "Heading 2")
    insert_before(
        section7,
        "为避免把更宽的上皮口径误写成严格恶性证据，新增患者结果在三套既有定义下并列复核：严格层要求冻结的保守CNV支持；扩展层接受inferCNV高负荷或CopyKAT非整倍体任一阳性；暂定层仅要求肿瘤样上皮注释。各层阈值均未调整，结果按患者汇总。参考文献提出的CAF—TAM—CD8方向则拆成肿瘤→CAF、CAF→TAM、TAM状态和TAM→T细胞四段逐段筛查。",
    )
    add_labeled_paragraph_before(
        section7,
        "三层复核结果  ",
        "GSE178481全队列TYMP—SPP1在严格层仅2名患者，无法估计；扩展层5名患者呈正向但未显著（rho=0.700，P=0.188）；暂定层8名患者呈显著正相关（rho=0.929，精确P=0.00223，层内BH校正q=0.0179）。因此该结果属于依赖肿瘤样口径的正向支持，不能作为严格恶性细胞确认。GSE207493在严格、扩展和暂定三层均未形成正向队列证据。",
    )
    add_labeled_paragraph_before(
        section7,
        "参考桥接结果  ",
        "GSE264586中POSTN⁺ CAF与APOE⁺ TAM在5名可评价患者中rho=1.000，精确P=0.0167；GSE178481和GSE202374均支持无基因重叠的SPP1髓系状态与髓系免疫抑制程序同向。APOE⁺ TAM与免疫抑制评分虽在多队列相关，但两套评分共享APOE、C1QC、CD163和MRC1，只能解释为状态一致性，不能当作独立机制复现。",
    )
    add_labeled_paragraph_before(
        section7,
        "证据边界  ",
        "现有数据未支持F1/TYMP→POSTN⁺ CAF，也未支持TAM→T细胞功能障碍、细胞毒性下降或血管重塑。GSE207493新增患者中的三项免疫关系P=0.02345、q=0.09380，属于名义显著的小样本次级结果，不应写作FDR显著。当前最稳妥的结论是：TAM免疫抑制状态和局部CAF—TAM连接值得继续验证，但完整F1/TYMP→CAF→TAM→CD8链条尚未成立。",
    )

    supported = find_paragraph_starting(document, "第一，cNMF F1与TYMP在患者内总体同向")
    supported.add_run(
        " 第六，增量筛查进一步支持SPP1髓系免疫抑制状态，并在GSE264586观察到小样本POSTN⁺ CAF—APOE⁺ TAM局部桥接。"
    )

    unsupported = find_paragraph_starting(document, "F1不对应稳定的RCC遗传易感性富集。")
    unsupported.add_run(
        " 新增复核还表明，GSE178481的TYMP—SPP1关联依赖暂定肿瘤样口径；F1/TYMP→POSTN⁺ CAF及TAM→CD8功能改变均未得到支持，因此不能把局部桥接外推为完整因果链。"
    )

    conclusion = find_paragraph_starting(document, "综合结论  第一部分26项分析已经全部完成。")
    conclusion.add_run(
        " Phase 52–53增量复核未改变原有26项主分析结论，但进一步限定了恶性细胞口径与机制外推：TYMP—SPP1仅获暂定肿瘤样层支持，SPP1免疫抑制TAM状态及POSTN⁺ CAF—APOE⁺ TAM局部连接可作为后续验证方向，完整F1/TYMP→CAF→TAM→CD8链条目前不成立。"
    )

    document.core_properties.modified = datetime.now(timezone.utc)
    document.save(OUTPUT)

    text = "\n".join(paragraph.text for paragraph in document.paragraphs)
    required = [
        "C06  新增患者三层复核与参考桥接筛查",
        "rho=0.929",
        "精确P=0.0167",
        "q=0.09380",
        "完整F1/TYMP→CAF→TAM→CD8链条目前不成立",
    ]
    missing = [item for item in required if item not in text]
    if missing:
        raise RuntimeError(f"Required report text missing after save: {missing}")

    status = {
        "status": "REPORT_UPDATED_PENDING_RENDER_QC",
        "source": str(SOURCE.relative_to(ROOT)),
        "source_sha256": sha256(SOURCE),
        "output": str(OUTPUT.relative_to(ROOT)),
        "output_bytes": OUTPUT.stat().st_size,
        "output_sha256": sha256(OUTPUT),
        "paragraphs": len(document.paragraphs),
        "tables": len(document.tables),
        "required_text_checks": {item: item in text for item in required},
        "edits": [
            "updated version note and statistical rule",
            "added executive-summary incremental evidence paragraph",
            "added C06 three-tier and reference-bridge section",
            "updated supported, unsupported, and final-conclusion paragraphs",
        ],
    }
    (QC / "report_update_status_v1.json").write_text(
        json.dumps(status, ensure_ascii=False, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
