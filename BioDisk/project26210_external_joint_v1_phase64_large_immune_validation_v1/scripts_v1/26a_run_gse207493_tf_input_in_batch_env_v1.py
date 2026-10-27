#!/usr/bin/env python3
"""Environment bridge: prepare sparse input in batch_bench; loom is built separately."""

from __future__ import annotations

import runpy
import sys
import types
from pathlib import Path


stub = types.ModuleType("loompy")
stub.create = lambda *args, **kwargs: None
sys.modules["loompy"] = stub
script = Path(__file__).with_name("26_prepare_gse207493_tf_inputs_v1.py")
runpy.run_path(str(script), run_name="__main__")
