from __future__ import annotations

import sys
from pathlib import Path

from PIL import Image, ImageDraw


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: 07_make_contact_sheets_v1.py <render_dir> <output_dir>")
    source = Path(sys.argv[1])
    output = Path(sys.argv[2])
    output.mkdir(parents=True, exist_ok=True)
    pages = sorted(source.glob("page-*.png"))
    if not pages:
        raise RuntimeError(f"No rendered pages found in {source}")

    thumb_width = 760
    label_height = 34
    columns = 2
    pages_per_sheet = 4
    with Image.open(pages[0]) as first:
        ratio = thumb_width / first.width
        thumb_height = int(round(first.height * ratio))

    for sheet_index, start in enumerate(range(0, len(pages), pages_per_sheet), start=1):
        group = pages[start : start + pages_per_sheet]
        rows = (len(group) + columns - 1) // columns
        sheet = Image.new(
            "RGB", (columns * thumb_width, rows * (thumb_height + label_height)), "white"
        )
        draw = ImageDraw.Draw(sheet)
        for index, page in enumerate(group):
            with Image.open(page) as image:
                thumb = image.convert("RGB").resize(
                    (thumb_width, thumb_height), Image.Resampling.LANCZOS
                )
            x = (index % columns) * thumb_width
            y = (index // columns) * (thumb_height + label_height)
            sheet.paste(thumb, (x, y + label_height))
            draw.text((x + 10, y + 8), page.stem, fill="#202225")
        sheet.save(output / f"render_contact_sheet_{sheet_index:02d}.png")


if __name__ == "__main__":
    main()
