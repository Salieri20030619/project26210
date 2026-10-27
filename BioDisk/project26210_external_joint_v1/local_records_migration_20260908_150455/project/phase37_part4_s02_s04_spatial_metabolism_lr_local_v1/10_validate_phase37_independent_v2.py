#!/usr/bin/env python3
"""Tie-aware exact signed-rank validator, independent of SciPy version defaults."""
import importlib.util,itertools
from pathlib import Path
import numpy as np
from scipy.stats import rankdata

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("v1",HERE/"09_validate_phase37_independent_v1.py"); v1=importlib.util.module_from_spec(spec); spec.loader.exec_module(v1)
def exact_signed_rank(v):
 v=np.asarray(v,float); v=v[np.isfinite(v)]; v=v[v!=0]
 if len(v)<3:return 1.
 ranks=rankdata(np.abs(v),method="average"); observed=abs(np.sum(np.sign(v)*ranks))
 distribution=np.fromiter((abs(np.sum(np.asarray(s)*ranks)) for s in itertools.product((-1.,1.),repeat=len(v))),dtype=float)
 return float(np.mean(distribution>=observed-1e-12))
v1.sr=exact_signed_rank; v1.OUT=HERE/"independent_validation_local_v2"; raise SystemExit(v1.main())
