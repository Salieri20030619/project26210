#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import os
import re
from pathlib import Path


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1"
OUT = ROOT / "09_metabolism_trajectory/trajectory_v1"
MARKER = "## 2026-08-05 Part 1 TF、代谢终态与拟时序启动"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def correct_document(path: Path, pid: int, started_at: str, socket: str) -> None:
    text = path.read_text(encoding="utf-8")
    if text.count(MARKER) != 1:
        raise RuntimeError(f"Unexpected trajectory marker count in {path}")
    replacement = (
        f"- 拟时序后台于 `{started_at}` 由独立 tmux socket `{socket}` 启动，pane PID `{pid}`；"
        f"正式输出目录为 `{OUT}`，主日志为 `{OUT / 'logs/background_master_v1.log'}`。"
        "启动时状态为 running，最终 `RUN_COMPLETE` 只有在自动核验全部通过后才会生成；正式图仍需人工查看后补终验记录。\n"
        "- 启动审计：首次尝试因 R 包装器不接受 `--version` 在读入矩阵前安全退出；第二次已完成输入准备，但 R `fread` 把 scVI latent 的空索引表头误识别为数据行，"
        "造成 11,951/11,950 行不一致并安全退出。现已显式指定 `header=TRUE`，两次失败记录均保留在 `trajectory_v1/provenance`，正式计算改由 tmux 持久托管。"
    )
    pattern = re.compile(r"- 拟时序后台于 `[^`]+` 启动，PID `\d+`；正式输出目录为 `[^`]+`，主日志为 `[^`]+`。启动时状态为 running，最终 `RUN_COMPLETE` 只有在自动核验全部通过后才会生成；正式图仍需人工查看后补终验记录。")
    text, count = pattern.subn(replacement, text)
    if count != 1:
        raise RuntimeError(f"Could not replace launch paragraph in {path}: {count}")
    temporary = path.with_name(path.name + ".trajectory_correct_tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def correct_manifest(pid: int) -> None:
    path = PROJECT / "RESULTS_MANIFEST.tsv"
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    matches = 0
    status = OUT / "background_status_v1.tsv"
    for i in range(1, len(lines)):
        row = lines[i].split("\t")
        if row[0] == "p1_trajectory_launch_v1_20260805":
            matches += 1
            row[12] = str(status.stat().st_size)
            row[14] = f"tmux pane PID={pid};running at corrected record update"
            row[17] = f"SHA-256={sha256(status)};tmux-persisted cross-sectional trajectory;automatic and manual terminal verification pending."
            lines[i] = "\t".join(row)
    if matches != 1 or len(header) != 18:
        raise RuntimeError(f"Unexpected trajectory manifest matches/schema: {matches}/{len(header)}")
    temporary = path.with_name(path.name + ".trajectory_correct_tmp")
    temporary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pid", required=True, type=int)
    parser.add_argument("--started-at", required=True)
    parser.add_argument("--socket", required=True)
    args = parser.parse_args()
    for path in (PROJECT / "PROJECT_STATUS.md", ROOT / "PART1_REBUILD_STATUS.md", PROJECT / "CODEX_HANDOFF.md"):
        correct_document(path, args.pid, args.started_at, args.socket)
    correct_manifest(args.pid)
    print(f"PART1_TRAJECTORY_LAUNCH_RECORDS_CORRECTED pid={args.pid}")


if __name__ == "__main__":
    main()
