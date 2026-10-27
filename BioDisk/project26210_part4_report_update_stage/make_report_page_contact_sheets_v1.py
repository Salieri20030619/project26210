#!/usr/bin/env python3
from pathlib import Path
from PIL import Image, ImageDraw

SOURCE = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_report_full_update_stage/page_png_v2")
OUT = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_report_full_update_stage/contact_sheets_v2")
OUT.mkdir(parents=True, exist_ok=True)

pages = sorted(SOURCE.glob("page-*.png"), key=lambda p: int(p.stem.split("-")[-1]))
cover = next(SOURCE.glob("cover-*.png"))
groups = [[cover], pages[0:5], pages[5:10], pages[10:15], pages[15:20]]
for index, group in enumerate(groups, 1):
    thumbs = []
    for path in group:
        image = Image.open(path).convert("RGB")
        image.thumbnail((510, 660), Image.Resampling.LANCZOS)
        tile = Image.new("RGB", (530, 700), "white")
        tile.paste(image, ((530 - image.width) // 2, 30))
        ImageDraw.Draw(tile).text((12, 8), path.stem, fill="black")
        thumbs.append(tile)
    sheet = Image.new("RGB", (530 * len(thumbs), 700), "#d9d9d9")
    for column, tile in enumerate(thumbs):
        sheet.paste(tile, (530 * column, 0))
    sheet.save(OUT / f"report_pages_contact_{index:02d}_v1.png", dpi=(150, 150))
print(f"created={len(groups)}")
