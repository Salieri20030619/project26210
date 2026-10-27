from __future__ import annotations

import sys
from pathlib import Path

from PIL import Image


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: 09_make_preview_v1.py <input_image> <output_jpg>")
    source = Path(sys.argv[1])
    output = Path(sys.argv[2])
    output.parent.mkdir(parents=True, exist_ok=True)
    with Image.open(source) as image:
        rgb = image.convert("RGB")
        rgb.thumbnail((1360, 1760), Image.Resampling.LANCZOS)
        rgb.save(output, quality=84, optimize=True, progressive=True)


if __name__ == "__main__":
    main()
