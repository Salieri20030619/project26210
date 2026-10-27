#!/usr/bin/env python3
from pathlib import Path

from PIL import Image, ImageDraw


ROOT = Path("/mnt/BioDisk/wangmingyuan/project26210_tymp_curated_report_stage_v2")
PAGE_DIR = ROOT / "page_png_v2"
OUT_DIR = ROOT / "contact_sheets_v2"


def build_sheet(paths: list[Path], output: Path) -> None:
    thumb_width = 510
    margin = 28
    label_height = 34
    thumbs = []
    for path in paths:
        with Image.open(path) as src:
            rgb = src.convert("RGB")
            height = round(rgb.height * thumb_width / rgb.width)
            thumbs.append((path, rgb.resize((thumb_width, height), Image.Resampling.LANCZOS)))

    cell_height = max(image.height for _, image in thumbs) + label_height
    sheet = Image.new("RGB", (margin * 3 + thumb_width * 2, margin * 3 + cell_height * 2), "white")
    draw = ImageDraw.Draw(sheet)
    for index, (path, image) in enumerate(thumbs):
        row, col = divmod(index, 2)
        x = margin + col * (thumb_width + margin)
        y = margin + row * (cell_height + margin)
        draw.text((x, y), path.stem, fill="black")
        sheet.paste(image, (x, y + label_height))
        draw.rectangle(
            (x, y + label_height, x + image.width - 1, y + label_height + image.height - 1),
            outline="#777777",
            width=2,
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    sheet.save(output, quality=92)


def main() -> None:
    pages = sorted(PAGE_DIR.glob("page-*.png"))
    if not pages:
        raise SystemExit("No rendered pages found")
    for start in range(0, len(pages), 4):
        batch = pages[start : start + 4]
        build_sheet(batch, OUT_DIR / f"pages_{start + 1:02d}_{start + len(batch):02d}.jpg")


if __name__ == "__main__":
    main()
