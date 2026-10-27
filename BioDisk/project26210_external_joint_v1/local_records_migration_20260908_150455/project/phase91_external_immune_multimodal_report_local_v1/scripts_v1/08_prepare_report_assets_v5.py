from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

from PIL import Image


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit(
            "Usage: 08_prepare_report_assets_v5.py <main_figure_png> <assets_v4> <output_dir>"
        )
    source = Path(sys.argv[1])
    assets_v4 = Path(sys.argv[2])
    output = Path(sys.argv[3])
    if output.exists():
        raise SystemExit(f"Refusing existing output: {output}")
    shutil.copytree(assets_v4, output)

    image = Image.open(source).convert("RGB")
    target = output / "panel_B_three_cohort_atlas.png"
    crop_box = (60, 660, 2110, 1960)
    crop = image.crop(crop_box)
    crop.save(target, optimize=True)
    status = {
        "status": "PASS",
        "source": str(source),
        "source_sha256": sha256(source),
        "adjustment": (
            "Panel B top and right boundaries expanded to retain its full title and annotation; "
            "the lower boundary still retains the complete legend."
        ),
        "panel_B_crop_box": list(crop_box),
        "panel_B_width": crop.width,
        "panel_B_height": crop.height,
        "panel_B_sha256": sha256(target),
    }
    (output / "ASSET_ADJUSTMENT_v5.json").write_text(
        json.dumps(status, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
