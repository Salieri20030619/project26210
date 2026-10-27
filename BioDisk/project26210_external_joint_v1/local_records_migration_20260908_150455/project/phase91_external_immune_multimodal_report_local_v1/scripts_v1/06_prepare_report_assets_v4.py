from __future__ import annotations

import hashlib
import json
import shutil
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
    if len(sys.argv) != 4:
        raise SystemExit("Usage: 06_prepare_report_assets_v4.py <main_figure_png> <assets_v3> <output_dir>")
    source = Path(sys.argv[1])
    assets_v3 = Path(sys.argv[2])
    output = Path(sys.argv[3])
    if output.exists():
        raise SystemExit(f"Refusing existing output: {output}")
    shutil.copytree(assets_v3, output)
    image = Image.open(source).convert("RGB")
    target = output / "panel_B_three_cohort_atlas.png"
    crop_box = (60, 720, 2020, 1960)
    crop = image.crop(crop_box)
    crop.save(target, optimize=True)
    status = {
        "status": "PASS",
        "source": str(source),
        "source_sha256": sha256(source),
        "adjustment": "Panel B lower boundary expanded to retain the complete legend without entering Panel D",
        "panel_B_crop_box": list(crop_box),
        "panel_B_width": crop.width,
        "panel_B_height": crop.height,
        "panel_B_sha256": sha256(target),
    }
    (output / "ASSET_ADJUSTMENT_v4.json").write_text(json.dumps(status, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
