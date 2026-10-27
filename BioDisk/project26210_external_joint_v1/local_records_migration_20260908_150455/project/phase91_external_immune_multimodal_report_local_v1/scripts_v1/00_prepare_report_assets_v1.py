from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

from PIL import Image


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: 00_prepare_report_assets_v1.py <main_figure_png> <output_dir>")
    source = Path(sys.argv[1])
    output = Path(sys.argv[2])
    if output.exists():
        raise SystemExit(f"Refusing existing output: {output}")
    output.mkdir(parents=True)
    image = Image.open(source).convert("RGB")
    w, h = image.size
    if (w, h) != (4058, 5622):
        raise SystemExit(f"Unexpected main figure size: {(w, h)}")

    crops = {
        "main_figure_v12_full.png": (0, 0, w, h),
        "panel_A_analysis_architecture.png": (80, 0, w - 80, 900),
        "panel_B_three_cohort_atlas.png": (60, 820, 2080, 2720),
        "panel_C_cellular_source_programs.png": (1940, 820, w - 50, 2720),
        "panel_D_cross_cohort_replication.png": (50, 1900, w - 50, 3820),
        "panel_E_myeloid_ordering.png": (50, 3600, 2100, h - 30),
        "panel_F_spatial_communication.png": (1960, 3600, w - 30, h - 30),
        "panel_F_spatial_maps.png": (1960, 3600, w - 30, 4750),
        "panel_F_cellchat_chord_specificity.png": (1960, 4480, w - 30, h - 30),
    }
    records = []
    for name, box in crops.items():
        target = output / name
        crop = image.crop(box)
        crop.save(target, optimize=True)
        records.append({
            "asset": name,
            "crop_box": list(box),
            "width": crop.width,
            "height": crop.height,
            "bytes": target.stat().st_size,
            "sha256": sha256(target),
        })
    (output / "ASSET_MANIFEST_v1.json").write_text(
        json.dumps({"source": str(source), "source_sha256": sha256(source), "assets": records}, indent=2),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
