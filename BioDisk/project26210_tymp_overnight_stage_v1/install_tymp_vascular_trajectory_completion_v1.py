#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import shutil
from datetime import datetime
from pathlib import Path


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
STAGE_V = Path("/mnt/BioDisk/wangmingyuan/project26210_tymp_vascular_niche_stage_v1")
STAGE_T = Path("/mnt/BioDisk/wangmingyuan/project26210_tymp_overnight_stage_v1")
DEST_V = ROOT / "results/part4_spatial/tymp_vascular_niche_v1"
DEST_T = ROOT / "results/part4_spatial/tymp_myeloid_state_ordering_v1"
DEST_S = ROOT / "results/part4_spatial/tymp_vascular_trajectory_synthesis_v1"
DOCS = {
    ROOT / "CODEX_HANDOFF.md": "c6f5160343cbf0b521e6d7be54b4bf246f9d3b9374e39afb1a6d9665a0abbab6",
    ROOT / "PROJECT_STATUS.md": "275f3bee090f52e197add891f00f405f85829460c4c58c2acdeab7636d2d91e7",
    ROOT / "RESULTS_MANIFEST.tsv": "8337a9610cad404398ed654ac614c6a61082f530ef819ec43f9f2226630f335e",
}
BACKUP = STAGE_T / "project_records_before_vascular_trajectory_install_v1"


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_file(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def manifest_row(result_id: str, sample: str, analysis: str, artifact_type: str, path: Path,
                 generator: str, inputs: str, parameters: str, version: str, qc_evidence: str,
                 notes: str) -> list[str]:
    return [
        result_id, "Part4", sample, analysis, artifact_type, str(path), generator, inputs,
        "batch_bench", parameters, version, "2026-08-04", str(path.stat().st_size), "passed",
        qc_evidence, "not_used", "", f"SHA-256={sha256(path)};{notes}",
    ]


def main() -> None:
    for path, expected in DOCS.items():
        observed = sha256(path)
        if observed != expected:
            raise RuntimeError(f"Project record changed since preflight: {path}; {observed} != {expected}")
    for destination in [DEST_V, DEST_T, DEST_S, BACKUP]:
        if destination.exists():
            raise RuntimeError(f"Refusing to overwrite: {destination}")
    if json.loads((STAGE_T / "joint_qc_v2/independent_verification_summary_v2.json").read_text())["status"] != "passed":
        raise RuntimeError("Numeric verification is not passed")
    if json.loads((STAGE_T / "joint_qc_v2/joint_figure_verification_summary_v1.json").read_text())["status"] != "passed":
        raise RuntimeError("Figure verification is not passed")

    BACKUP.mkdir(parents=True)
    for path in DOCS:
        shutil.copy2(path, BACKUP / path.name)

    shutil.copytree(STAGE_V / "output_v1", DEST_V)
    (DEST_V / "provenance").mkdir()
    for source in [
        STAGE_V / "run_tymp_vascular_niche_analysis_v1.py",
        STAGE_V / "run_tymp_vascular_niche_background_v1.sh",
        STAGE_V / "background_status_v1.tsv",
        STAGE_V / "logs/analysis_v1.log",
        STAGE_V / "logs/memory_gate_launcher_v1.log",
    ]:
        copy_file(source, DEST_V / "provenance" / source.name)

    shutil.copytree(STAGE_T / "trajectory_v1", DEST_T)
    (DEST_T / "provenance").mkdir()
    for source in [
        STAGE_T / "run_monocyte_spp1_trajectory_v1.py",
        STAGE_T / "run_overnight_v1.sh",
        STAGE_T / "overnight_status_v1.tsv",
        STAGE_T / "logs/trajectory_v1.log",
    ]:
        copy_file(source, DEST_T / "provenance" / source.name)

    DEST_S.mkdir(parents=True)
    shutil.copytree(STAGE_T / "joint_figures_v4", DEST_S / "figures")
    shutil.copytree(STAGE_T / "joint_qc_v2", DEST_S / "qc")
    shutil.copytree(STAGE_T / "preflight_verification_v1_failure", DEST_S / "preflight/v1_failures")
    for version in ["joint_figures_v1", "joint_figures_v2", "joint_figures_v3"]:
        shutil.copytree(STAGE_T / version, DEST_S / "preflight" / version)
    for source in [
        STAGE_T / "verify_tymp_vascular_trajectory_v1.py",
        STAGE_T / "render_tymp_vascular_trajectory_figures_v1.py",
        STAGE_T / "verify_joint_figures_v1.py",
    ]:
        copy_file(source, DEST_S / "scripts" / source.name)

    installed_files = sorted(path for base in [DEST_V, DEST_T, DEST_S] for path in base.rglob("*") if path.is_file())
    with (DEST_S / "installed_file_manifest_v1.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["path", "size_bytes", "sha256"])
        for path in installed_files:
            writer.writerow([str(path), path.stat().st_size, sha256(path)])

    timestamp = now()
    status_path = ROOT / "PROJECT_STATUS.md"
    status = status_path.read_text(encoding="utf-8")
    old_heading = "## 0. 当前权威状态（2026-08-03 17:13 +0800）"
    if old_heading not in status:
        raise RuntimeError("PROJECT_STATUS authoritative heading changed unexpectedly")
    status = status.replace(old_heading, "## 0. 当前权威状态（2026-08-04 10:16 +0800）", 1)
    insert = """

- **Part 4 TYMP血管生态位与髓系横断面状态排序已完成并通过独立验收。** 血管分析覆盖GLMF1/GLMF2共7,079 spots及TJ-RCC 12例34,609 spots，9条预设关系、4×4/5×5分块和每项9,999次置换；正式目录为`results/part4_spatial/tymp_vascular_niche_v1`。髓系状态排序覆盖参考对象中53,794个Monocyte/Other Macrophage/APOE+C1QC+ TAM/SPP1+ TAM细胞，63个达到样本内检验门槛的样本；正式目录为`results/part4_spatial/tymp_myeloid_state_ordering_v1`。
- **患者级稳定支持限于两条空间关系，不支持整条机制链。** `SPP1–vascular niche→T-cell balance`在4×4和5×5均为12/12患者预期负向，中位rho=`-0.456/-0.383`、全局BH FDR均=`0.00220`；reference-dependent `TYMP–MTA core→SPP1–vascular niche`两尺度均12/12正向，中位rho=`0.485/0.545`、FDR均=`0.00220`。SPP1 TAM→vascular compartment、cytotoxicity、hypoxia/angiogenesis、invasion/EMT及SPP1/monocyte balance均未获统一全局FDR支持。无治疗响应、生存或纵向标签，临床耐药不可估计。
- **DPT只支持横断面状态排序，不证明单核细胞向SPP1+ TAM谱系转化。** 以单核细胞根和scVI邻居图计算的DPT中，C1QC状态为57/63样本正向（方向FDR=`5.74e-11`），不含SPP1基因本身的SPP1-state为41/63正向（FDR=`0.0316`），而SPP1-state minus C1QC-state仅11/63正向（FDR=`3.90e-7`）。该结果更符合异质的髓系状态连续谱，不能替代RNA velocity、谱系示踪或干预实验。
- **独立数值验收35/35、正式图文件检查4/4、人工图检4/4通过。** 正式综合目录为`results/part4_spatial/tymp_vascular_trajectory_synthesis_v1`。验证器首轮2项预期错误、绘图v1中止、v2语义歧义和v3版式问题均保留在`preflight`；只有v4图集为正式版本。本轮未修改任何客户DOCX，现有报告仍为candidate、未外发。
"""
    marker = "\n- **用户于2026-08-03授权：若预计3天内可完成"
    if marker not in status:
        raise RuntimeError("PROJECT_STATUS insertion marker missing")
    status = status.replace(marker, insert + marker, 1)
    status_path.write_text(status, encoding="utf-8")

    handoff_path = ROOT / "CODEX_HANDOFF.md"
    handoff = handoff_path.read_text(encoding="utf-8")
    handoff += f"""

## 2026-08-04 Part 4 TYMP血管生态位与髓系状态排序终态

- 血管生态位正式目录：`{DEST_V}`；覆盖GLMF两例7,079 spots和TJ-RCC 12例34,609 spots。两种网格中稳定的患者级结果为：SPP1–vascular niche与T-cell balance 12/12负向（中位rho=-0.456/-0.383，FDR=0.00220），reference-dependent TYMP–MTA core与SPP1–vascular niche 12/12正向（中位rho=0.485/0.545，FDR=0.00220）。其余预设机制关系未形成统一FDR支持，临床耐药因无结局标签不可估计。
- 髓系横断面排序正式目录：`{DEST_T}`；53,794细胞，63个合格样本。C1QC-state 57/63随DPT上升，SPP1-state（排除SPP1本身）仅41/63弱上升，SPP1-minus-C1QC仅11/63正向。不得表述为单核细胞向SPP1+ TAM的谱系转化、RNA velocity或因果过程。
- 综合正式图和QC：`{DEST_S}`。独立数值验收35/35、图文件检查4/4、人工图检4/4通过；只有`figures`中的v4图集可正式使用。早期验证器和绘图失败均保留在`preflight`。
- 本轮环境为`batch_bench` Python 3.11.15；每个新步骤均通过50 GiB内存门禁。运行时`MemAvailable`约426–427 GiB；SwapFree约240 kB，记录为共享主机资源警告但未触发项目阻断规则。
- 本轮未修改客户DOCX。若要把这两项结果写入完整主报告或TYMP精选版，必须由用户明确授权后新建版本化候选并重新执行DOCX结构、渲染和人工版式验收。

## 本次交接更新时间

{timestamp}（Asia/Shanghai）
"""
    handoff_path.write_text(handoff, encoding="utf-8")

    manifest_path = ROOT / "RESULTS_MANIFEST.tsv"
    rows = [
        manifest_row("p4_tymp_vascular_patient_summary_v1_20260804", "12 TJ-RCC tumor cases",
                     "TYMP-SPP1 vascular-niche patient direction synthesis", "statistics_table",
                     DEST_V / "tables/tjrcc_cohort_direction_summary_v1.tsv",
                     str(DEST_V / "provenance/run_tymp_vascular_niche_analysis_v1.py"),
                     "14 spatial H5AD;TYMP state tables;BayesSpace computational boundaries",
                     "4x4 and5x5;9999 permutations;patient sign test;global BH", "v1",
                     str(DEST_S / "qc/independent_verification_v2.tsv"),
                     "12/12 negative T-cell-balance and 12/12 positive reference-dependent TYMP-niche routes;association only"),
        manifest_row("p4_tymp_vascular_spot_scores_v1_20260804", "GLMF1;GLMF2;12 TJ-RCC tumor cases",
                     "TYMP-SPP1 vascular-niche spot scores", "spot_table",
                     DEST_V / "tables/spot_mechanism_scores_v1.tsv.gz",
                     str(DEST_V / "provenance/run_tymp_vascular_niche_analysis_v1.py"),
                     "cell2location q05;expression signatures;computational boundary",
                     "41688 spots;residualized signatures;route definitions", "v1",
                     str(DEST_S / "qc/independent_verification_v2.tsv"),
                     "not independent biological replicates"),
        manifest_row("p4_tymp_myeloid_cell_scores_v1_20260804", "71 reference samples",
                     "myeloid cross-sectional DPT state ordering", "cell_table",
                     DEST_T / "tables/myeloid_trajectory_cell_scores_v1.tsv.gz",
                     str(DEST_T / "provenance/run_monocyte_spp1_trajectory_v1.py"),
                     "ref_signatures.h5ad;X_scVI;counts",
                     "53794 cells;30-NN;diffmap;DPT;monocyte root", "v1",
                     str(DEST_S / "qc/independent_verification_v2.tsv"),
                     "cross-sectional ordering;not velocity or lineage tracing"),
        manifest_row("p4_tymp_myeloid_direction_summary_v1_20260804", "63 eligible reference samples",
                     "sample-level DPT score direction synthesis", "statistics_table",
                     DEST_T / "tables/cohort_pseudotime_direction_summary_v1.tsv",
                     str(DEST_T / "provenance/run_monocyte_spp1_trajectory_v1.py"),
                     "myeloid_trajectory_cell_scores_v1.tsv.gz",
                     "within-sample Spearman;sample sign test;BH", "v1",
                     str(DEST_S / "qc/independent_verification_v2.tsv"),
                     "C1QC 57/63 positive;SPP1-state 41/63;SPP1-minus-C1QC 11/63"),
    ]
    for stem, sample, analysis in [
        ("vascular_niche_patient_effects_v1", "12 TJ-RCC tumor cases", "vascular-niche patient effects"),
        ("vascular_niche_cohort_summary_v1", "12 TJ-RCC tumor cases", "vascular-niche cohort summary"),
        ("myeloid_cross_sectional_state_overview_v2", "53794 reference myeloid cells", "myeloid cross-sectional state overview"),
        ("myeloid_state_patient_directions_v1", "63 eligible reference samples", "myeloid state direction synthesis"),
    ]:
        rows.append(manifest_row(f"p4_{stem}_20260804", sample, analysis, "main_figure",
                                 DEST_S / f"figures/{stem}.png",
                                 str(DEST_S / "scripts/render_tymp_vascular_trajectory_figures_v1.py"),
                                 "formal vascular or trajectory tables", "formal v4 PNG and PDF", "v4",
                                 str(DEST_S / "qc/joint_figure_verification_v1.tsv"),
                                 "manual review PASS;inference boundary embedded"))
    rows.extend([
        manifest_row("p4_tymp_vascular_trajectory_verification_v2_20260804", "joint Part4 extension",
                     "independent numeric and fingerprint verification", "qc_table",
                     DEST_S / "qc/independent_verification_v2.tsv",
                     str(DEST_S / "scripts/verify_tymp_vascular_trajectory_v1.py"),
                     "vascular and trajectory formal tables;14 current input hashes",
                     "recompute Spearman;sign tests;BH;counts;hashes", "v2",
                     "35/35", "corrected verifier;v1 expectation failure retained"),
        manifest_row("p4_tymp_vascular_trajectory_manual_v1_20260804", "4 formal figures",
                     "manual visual review", "manual_qc",
                     DEST_S / "qc/joint_figure_manual_review_v1.tsv",
                     "Codex visual inspection", "four formal v4 PNG",
                     "layout;labels;color;clipping;inference boundary", "v1",
                     "4/4", "only v4 figure set passed for formal use"),
    ])
    with manifest_path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)

    metadata = {
        "status": "completed",
        "completed_at": timestamp,
        "destinations": [str(DEST_V), str(DEST_T), str(DEST_S)],
        "installed_files": len(installed_files),
        "records_before": {str(path): expected for path, expected in DOCS.items()},
        "records_after": {str(path): sha256(path) for path in DOCS},
        "client_docx_modified": False,
    }
    (DEST_S / "installation_metadata_v1.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
