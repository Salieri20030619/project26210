#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import os
from datetime import datetime
from pathlib import Path


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1"
OUT = ROOT / "09_metabolism_trajectory/trajectory_v1"
MARKER = "## 2026-08-05 Part 1 拟时序 scTour 与 CellRank 恢复"
RESULT_ID = "p1_trajectory_sctour_cellrank_resume_v1_20260805"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def append_once(path: Path, block: str) -> None:
    text = path.read_text(encoding="utf-8")
    if MARKER in text:
        return
    tmp = path.with_name(path.name + ".trajectory_resume_tmp")
    tmp.write_text(text.rstrip() + "\n\n" + block.strip() + "\n", encoding="utf-8")
    os.replace(tmp, path)


def update_manifest(pid: int) -> None:
    path = PROJECT / "RESULTS_MANIFEST.tsv"
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines[0].split("\t")) != 18:
        raise RuntimeError("Unexpected RESULTS_MANIFEST schema")
    old_matches = 0
    for i in range(1, len(lines)):
        row = lines[i].split("\t")
        if row[0] == "p1_trajectory_launch_v1_20260805":
            old_matches += 1
            row[16] = RESULT_ID
            if "scTour initialization stopped" not in row[17]:
                row[17] += ";initial run completed R methods and CytoTRACE2, then scTour initialization stopped on missing n_genes_by_counts."
            lines[i] = "\t".join(row)
    if old_matches != 1:
        raise RuntimeError(f"Unexpected original launch row count: {old_matches}")
    if not any(line.startswith(RESULT_ID + "\t") for line in lines[1:]):
        status = OUT / "background_status_v1.tsv"
        row = [RESULT_ID, "Part1", "11950 cells", "resume scTour and pseudotime-informed CellRank", "background_run", str(status), str(ROOT / "09_metabolism_trajectory/run_part1_trajectory_master_v1.sh"), "completed Monocle3;Slingshot;CytoTRACE2 checkpoints", "part1_trajectory_v1", "scTour 200 epochs;CellRank 0.8 pseudotime plus 0.2 connectivity;seed26210", "v1", "2026-08-05", str(status.stat().st_size), "partial", f"tmux pane PID={pid};running at resume record update", "candidate", "", f"SHA-256={sha256(status)};reuses 11950-cell CytoTRACE2 checkpoint;automatic and manual terminal verification pending."]
        lines.append("\t".join(row))
    tmp = path.with_name(path.name + ".trajectory_resume_tmp")
    tmp.write_text("\n".join(lines) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pid", required=True, type=int)
    parser.add_argument("--started-at", required=True)
    parser.add_argument("--socket", required=True)
    args = parser.parse_args()
    now = datetime.now().astimezone().isoformat(timespec="seconds")
    block = f"""{MARKER}

- 已完成结果保持不变：Monocle3、Slingshot 和 5/5 研究 CytoTRACE2；CytoTRACE2 checkpoint 含 11,950/11,950 个有限评分，本次恢复不重算。
- 上次停止点为新建 trajectory AnnData 未生成 scTour 明确要求的 `n_genes_by_counts`。上游 QC 的 `n_genes_raw` 一直存在；本次按 scTour 实际使用的 2,500-HVG raw-count 矩阵计算 `n_genes_by_counts` 和 `total_counts`，不改变全局 QC 或筛选细胞。
- 新增 `SCTOUR_COMPLETE`、逐细胞 scTour 分数和 latent checkpoint；只有 scTour checkpoint 完整后才运行 pseudotime-informed CellRank，避免后续故障导致重复训练。
- 剩余两种方法于 `{args.started_at}` 在 tmux socket `{args.socket}` 恢复，pane PID `{args.pid}`。完成自动核验前状态仍为 `partial`；最终解释仍限于横断面状态排序和伪时间约束的 fate probability。

## 本次更新时间

{now}（Asia/Shanghai）"""
    for path in (PROJECT / "PROJECT_STATUS.md", ROOT / "PART1_REBUILD_STATUS.md", PROJECT / "CODEX_HANDOFF.md"):
        append_once(path, block)
    update_manifest(args.pid)
    print(f"PART1_TRAJECTORY_RESUME_RECORDS_UPDATED pid={args.pid}")


if __name__ == "__main__":
    main()
