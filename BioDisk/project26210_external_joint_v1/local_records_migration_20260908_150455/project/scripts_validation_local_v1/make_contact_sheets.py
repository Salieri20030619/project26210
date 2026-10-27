from __future__ import annotations

import sys
from pathlib import Path

from PIL import Image, ImageDraw


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("usage: make_contact_sheets.py PAGE_DIR OUTPUT_DIR")
    page_dir = Path(sys.argv[1]).resolve()
    output_dir = Path(sys.argv[2]).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    pages = sorted(page_dir.glob("page-*.png"))
    if not pages:
        raise SystemExit(f"no page PNG files in {page_dir}")
    for sheet_index, start in enumerate(range(0, len(pages), 2), start=1):
        batch = pages[start : start + 2]
        images = [Image.open(path).convert("RGB") for path in batch]
        label_height = 30
        gap = 20
        width = sum(image.width for image in images) + gap * (len(images) - 1)
        height = max(image.height for image in images) + label_height
        sheet = Image.new("RGB", (width, height), "#d9d9d9")
        draw = ImageDraw.Draw(sheet)
        x = 0
        for path, image in zip(batch, images):
            draw.text((x + 8, 7), path.stem, fill="black")
            sheet.paste(image, (x, label_height))
            x += image.width + gap
        target = output_dir / f"sheet-{sheet_index:03d}.png"
        sheet.save(target, optimize=True)
        for image in images:
            image.close()
    print(f"pages={len(pages)} sheets={(len(pages)+1)//2}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
