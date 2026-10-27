from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("immune_report_v2", SCRIPT_DIR / "03_build_report_v2.py")
v2 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(v2)
base = v2.base
original_add_figure = base.add_figure


def add_figure_v3(doc, path, width_in, caption, alt_title, alt_description):
    name = Path(path).name
    if name == "panel_D_cross_cohort_replication.png":
        width_in = 5.85
    elif name == "panel_E_myeloid_ordering.png":
        width_in = 4.90
    return original_add_figure(doc, path, width_in, caption, alt_title, alt_description)


base.add_figure = add_figure_v3


if __name__ == "__main__":
    base.main()
