#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import os
from datetime import datetime
from pathlib import Path


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1"
TRAJ = ROOT / "09_metabolism_trajectory/trajectory_v1"
MARKER = "## 2026-08-05 Part 1 TF、代谢终态与拟时序启动"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def append_once(path: Path, block: str) -> None:
    original = path.read_text(encoding="utf-8")
    if MARKER in original:
        return
    temporary = path.with_name(path.name + ".trajectory_record_tmp")
    temporary.write_text(original.rstrip() + "\n\n" + block.strip() + "\n", encoding="utf-8")
    os.replace(temporary, path)


def append_manifest(rows: list[list[str]]) -> None:
    path = PROJECT / "RESULTS_MANIFEST.tsv"
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    if len(header) != 18:
        raise RuntimeError(f"Unexpected manifest schema: {len(header)} columns")
    known = {line.split("\t", 1)[0] for line in lines[1:] if line}
    additions = [row for row in rows if row[0] not in known]
    if not additions:
        return
    if any(len(row) != len(header) for row in additions):
        raise RuntimeError("Manifest row does not match header")
    temporary = path.with_name(path.name + ".trajectory_record_tmp")
    temporary.write_text("\n".join(lines + ["\t".join(row) for row in additions]) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pid", required=True, type=int)
    parser.add_argument("--started-at", required=True)
    args = parser.parse_args()
    now = datetime.now().astimezone().isoformat(timespec="seconds")
    scenic = ROOT / "07_pyscenic_v1/intermediate/auc_grnboost2_all_cells_v1.csv"
    metaviper = ROOT / "08_metaviper_v1/rds/metaviper_activity_entrez_by_cell_v1.rds"
    metabolism = ROOT / "09_metabolism_trajectory/metabolism_v1/source_data/scFEA_four_axis_cell_scores_v1.tsv.gz"
    env_manifest = Path("/mnt/BioDisk/wangmingyuan/project26210_envs/part1_trajectory_v1/file_manifest_v1.tsv")
    background = TRAJ / "background_status_v1.tsv"
    for path in (scenic, metaviper, metabolism, env_manifest, background):
        if not path.exists() or path.stat().st_size == 0:
            raise FileNotFoundError(path)

    common = f"""{MARKER}

- pySCENIC 双路线核心计算已于 2026-08-05 13:59 完成，`PYSCENIC_CORE_COMPLETE` 与 `OVERNIGHT_COMPLETE` 均存在；GENIE3/GRNBoost2 的邻接、regulon 和 11,950 细胞 AUCell 结果已落盘。当前只确认核心运行完成，尚未完成独立终验，状态保持 `partial`。
- metaVIPER 已于 2026-08-05 13:59 完成：5,566 个 regulator、11,921 个有可估计活性的细胞、20/20 个达到门槛的患者，日志终态为 `PART1_METAVIPER_V1_PASS`。独立终验尚待补，状态保持 `partial`。
- 代谢分析已完成 11,950 个细胞的四通路 AUCell、scFEA M168（168 modules，100 epochs）和患者配对统计；独立验证 21/21 通过。MTA/Kennedy 仅为转录/通量代理，没有匹配代谢组，不得解释为 5′-MTA 实测。
- 拟时序正式输入固定为 `eligible_for_subclone_selection_v1 == True` 的 11,950 个保守恶性清晰上皮细胞（26 名患者、5 个研究）。根方向锚点仅用患者内标准化的 CNV burden、hypoxia、cycle 与 base program，明确排除 TYMP、cNMF F1 和代谢分数，避免循环论证。
- 拟时序采用 Monocle3、Slingshot、CytoTRACE2、scTour 及 pseudotime-informed CellRank；复用冻结的 `X_scVI`/`X_umap`。由于没有 spliced/unspliced 层，结果只表示横断面状态排序和伪时间约束的 fate probability，不是 RNA velocity、谱系追踪或真实时间。
- 独立环境 `/mnt/BioDisk/wangmingyuan/project26210_envs/part1_trajectory_v1` 已通过 6/6 环境核验及 Monocle3/dplyr 两种加载顺序的 `graph_test` 兼容性测试；历史故障根因是 `libudunits2.so.0` 运行时路径缺失，不是 dplyr 本身冲突。未修改共享 Conda/R 环境。
- 拟时序后台于 `{args.started_at}` 启动，PID `{args.pid}`；正式输出目录为 `{TRAJ}`，主日志为 `{TRAJ / 'logs/background_master_v1.log'}`。启动时状态为 running，最终 `RUN_COMPLETE` 只有在自动核验全部通过后才会生成；正式图仍需人工查看后补终验记录。

## 本次更新时间

{now}（Asia/Shanghai）"""

    append_once(PROJECT / "PROJECT_STATUS.md", common)
    append_once(ROOT / "PART1_REBUILD_STATUS.md", common)
    append_once(PROJECT / "CODEX_HANDOFF.md", common)

    rows = [
        ["p1_pyscenic_core_v1_20260805", "Part1", "11950 malignant epithelial cells", "pySCENIC dual-GRN core run", "activity_matrix", str(scenic), str(ROOT / "07_pyscenic_v1/run_part1_scenic_metaviper_overnight_v1.sh"), "verified pySCENIC loom;human cisTarget databases", "part1_scenic_v1", "GENIE3 and GRNBoost2;AUCell all cells", "v1", "2026-08-05", str(scenic.stat().st_size), "partial", str(ROOT / "07_pyscenic_v1/PYSCENIC_CORE_COMPLETE"), "candidate", "", f"SHA-256={sha256(scenic)};core computation complete;independent terminal verification pending."],
        ["p1_metaviper_activity_v1_20260805", "Part1", "11921 scored cells;5566 regulators", "metaVIPER regulator activity", "activity_matrix", str(metaviper), str(ROOT / "08_metaviper_v1/run_part1_metaviper_v1.R"), "patient expression matrices;kidney ARACNe interactomes", "R_4.5.2;viper1.44;aracne.networks1.36", "20 eligible patients;metaVIPER", "v1", "2026-08-05", str(metaviper.stat().st_size), "partial", str(ROOT / "08_metaviper_v1/RUN_COMPLETE"), "candidate", "", f"SHA-256={sha256(metaviper)};20/20 eligible patients;independent terminal verification pending."],
        ["p1_metabolism_four_axis_v1_20260805", "Part1", "11950 malignant epithelial cells", "four-pathway AUCell and scFEA M168", "cell_table", str(metabolism), str(ROOT / "09_metabolism_trajectory/run_part1_metabolism_background_v1.sh"), "raw counts;KEGG pathway definitions;scFEA M168", "part1_metabolism_v1", "168 modules;100 epochs;patient-paired statistics", "v1", "2026-08-05", str(metabolism.stat().st_size), "passed", "21/21 independent checks", "candidate", "", f"SHA-256={sha256(metabolism)};MTA/Kennedy is a proxy without matched metabolomics."],
        ["p1_trajectory_environment_v1_20260805", "Part1", "Monocle3;Slingshot;CytoTRACE2;scTour;CellRank", "isolated trajectory environment", "environment_manifest", str(env_manifest), "/mnt/BioDisk/wangmingyuan/project26210_trajectory_env_stage_v1/install_part1_trajectory_env_v1.sh", "private R/Python trajectory environment", "R4.5.3;Monocle3 1.4.27;Slingshot 2.18", "two package load orders;graph_test compatibility", "v1", "2026-08-05", str(env_manifest.stat().st_size), "passed", "6/6 plus two load-order compatibility tests", "not_used", "", f"SHA-256={sha256(env_manifest)};shared environments unchanged;runtime libudunits2 path fixed in wrapper."],
        ["p1_trajectory_launch_v1_20260805", "Part1", "11950 cells;26 patients;5 studies", "multi-method cross-sectional trajectory", "background_run", str(background), str(ROOT / "09_metabolism_trajectory/run_part1_trajectory_master_v1.sh"), "CNV-annotated H5AD;cNMF;CNV clones;metabolism scores", "part1_trajectory_v1", "Monocle3;Slingshot;CytoTRACE2;scTour;CellRank;seed26210", "v1", "2026-08-05", str(background.stat().st_size), "partial", f"PID={args.pid};running at record update", "candidate", "", f"SHA-256={sha256(background)};cross-sectional state ordering only;automatic and manual terminal verification pending."],
    ]
    append_manifest(rows)
    print(f"PART1_RECORDS_UPDATED trajectory_pid={args.pid}")


if __name__ == "__main__":
    main()
