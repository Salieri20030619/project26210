#!/usr/bin/env python3
"""Final audit wrapper: complete inventory, all Office content, targeted plain-text content."""
import os
import re
import runpy
from pathlib import Path

EXCLUDED_DIRS = {
    ".git", "__pycache__", ".cache", "env", "envs", ".venv", "venv",
    "site-packages", "node_modules", "_deps", "miniconda3", "conda-bld",
}
EVIDENCE_PATH = re.compile(r"(?i)mta|metabol|tymp|spp1|mtap|polyamine|44(?:x|fold|倍)|质谱|代谢组|非靶|湿实验|wet.?lab")
BASELINE_NAMES = {"project_status.md", "results_manifest.tsv", "codex_handoff.md"}
ORIGINAL_WALK = os.walk
ORIGINAL_READ_BYTES = Path.read_bytes


def evidence_walk(top, *args, **kwargs):
    for root, dirs, files in ORIGINAL_WALK(top, *args, **kwargs):
        dirs[:] = [d for d in dirs if d.lower() not in EXCLUDED_DIRS]
        yield root, dirs, files


def targeted_read_bytes(self):
    path_text = str(self).replace("\\", "/")
    if self.name.lower() in BASELINE_NAMES or EVIDENCE_PATH.search(path_text):
        return ORIGINAL_READ_BYTES(self)
    return b""


os.walk = evidence_walk
Path.read_bytes = targeted_read_bytes
runpy.run_path(str(Path(__file__).with_name("01_inventory_and_evidence_scan_v1.py")), run_name="__main__")
