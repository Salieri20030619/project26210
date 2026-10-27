from __future__ import annotations

import pickle
from pathlib import Path


ROOT = Path("/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/07_pyscenic_v1/intermediate")
FILES = {
    "GENIE3": ROOT / "regulons_genie3_v1.dat",
    "GRNBoost2": ROOT / "regulons_grnboost2_v1.dat",
}
TARGET_TFS = {"ILF2", "HOXD3", "HNF4A"}

print("method\ttf\tregulon\ttarget\tweight")
for method, path in FILES.items():
    with path.open("rb") as handle:
        regulons = pickle.load(handle)
    for regulon in regulons:
        tf = getattr(regulon, "transcription_factor", None)
        if tf is None:
            tf = str(getattr(regulon, "name", "")).split("(")[0]
        if tf not in TARGET_TFS:
            continue
        gene2weight = getattr(regulon, "gene2weight", {})
        for target, weight in sorted(gene2weight.items(), key=lambda item: abs(float(item[1])), reverse=True):
            print(f"{method}\t{tf}\t{getattr(regulon, 'name', tf)}\t{target}\t{float(weight):.10g}")
