#!/usr/bin/env python3
"""NumPy compatibility wrapper for the frozen loom builder."""

from __future__ import annotations

import runpy
from pathlib import Path

import numpy as np


original_upper = np.char.upper
np.char.upper = lambda values: original_upper(np.asarray(values, dtype=str))
script = Path(__file__).with_name("28_build_gse207493_scenic_loom_v1.py")
runpy.run_path(str(script), run_name="__main__")
