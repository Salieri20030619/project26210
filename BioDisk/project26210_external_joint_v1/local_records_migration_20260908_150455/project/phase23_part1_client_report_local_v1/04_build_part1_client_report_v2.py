from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "03_build_part1_client_report_v1.py"
SPEC = importlib.util.spec_from_file_location("part1_client_report_builder", SOURCE)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Cannot load builder: {SOURCE}")
builder = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(builder)

builder.ASSETS = ROOT / "assets_v2"
builder.OUTPUT = ROOT / "PART1_EXTERNAL_SINGLE_CELL_CLIENT_REPORT_v2.docx"
builder.STATUS = ROOT / "build_status_v2.json"
builder.build_document()
