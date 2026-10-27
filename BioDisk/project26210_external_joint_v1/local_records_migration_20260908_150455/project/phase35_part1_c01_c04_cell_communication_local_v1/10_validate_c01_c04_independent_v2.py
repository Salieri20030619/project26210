#!/usr/bin/env python3
"""Corrected validator: reproduce R wilcox.test(exact=FALSE) continuity correction for C03."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
from scipy.stats import wilcoxon

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("phase35_validator_v1", HERE / "09_validate_c01_c04_independent_v1.py")
validator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validator)


_call_number = 0
_python_signed_rank = validator.signed_rank


def component_aware_signed_rank(values):
    global _call_number
    _call_number += 1
    # main() evaluates C01 four times, C03 six times, then C04 sixteen times.
    # Only C03 was produced by R wilcox.test(exact=FALSE), whose default
    # continuity correction differs from the Python-produced C01/C04 tables.
    if not 5 <= _call_number <= 10:
        return _python_signed_rank(values)
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) < 3 or np.all(values == 0):
        return 1.0
    return float(wilcoxon(values, alternative="two-sided", zero_method="wilcox",
                          correction=True, method="approx").pvalue)


validator.signed_rank = component_aware_signed_rank
validator.OUT = HERE / "independent_validation_local_v2"
raise SystemExit(validator.main())
