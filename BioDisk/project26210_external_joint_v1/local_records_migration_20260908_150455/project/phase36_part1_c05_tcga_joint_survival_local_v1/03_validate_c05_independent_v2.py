#!/usr/bin/env python3
"""Corrected validator: allow the documented 1.96-vs-qnorm(.975) CI rounding difference."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("phase36_validator_v1", HERE / "02_validate_c05_independent_v1.py")
validator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validator)


def coefficient_aware_close(a, b, atol=1e-8):
    a = np.asarray(a, float); b = np.asarray(b, float)
    # Cox tables contain 96 or 108 rows. R uses qnorm(.975)=1.9599639845,
    # while the replay formula in v1 used the conventional rounded 1.96.
    if a.size in (96, 108):
        return np.allclose(a, b, rtol=5e-5, atol=max(atol, 1e-8), equal_nan=True)
    return np.allclose(a, b, rtol=1e-7, atol=atol, equal_nan=True)


validator.close = coefficient_aware_close
validator.OUT = HERE / "independent_validation_local_v2"
raise SystemExit(validator.main())
