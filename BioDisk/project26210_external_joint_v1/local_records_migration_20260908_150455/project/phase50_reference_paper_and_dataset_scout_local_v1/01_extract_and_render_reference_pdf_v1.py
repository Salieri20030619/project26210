#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import os
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEPS = HERE / "_deps"
sys.path.insert(0, str(DEPS))

import pandas as pd
import pypdfium2 as pdfium
from PIL import Image, ImageOps, ImageDraw
from pypdf import PdfReader


if len(sys.argv) != 2:
    raise SystemExit("Usage: python 01_extract_and_render_reference_pdf_v1.py <source.pdf>")
SOURCE = Path(sys.argv[1]).resolve()
FINAL = HERE / "paper_review_v1"
PART = HERE / "paper_review_v1.part"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


if FINAL.exists() or PART.exists():
    raise FileExistsError("Refusing existing paper-review output")
if not SOURCE.exists():
    raise FileNotFoundError(SOURCE)
PART.mkdir()
(PART / "pages").mkdir()
(PART / "contact_sheets").mkdir()

reader = PdfReader(SOURCE)
pdf = pdfium.PdfDocument(str(SOURCE))
if len(reader.pages) != len(pdf):
    raise RuntimeError("Page-count disagreement between PDF readers")

page_rows = []
text_sections = []
rendered = []
for index, page in enumerate(reader.pages, start=1):
    text = page.extract_text() or ""
    text_sections.append(f"\n===== PDF PAGE {index} =====\n{text.rstrip()}\n")
    nonempty = [line.strip() for line in text.splitlines() if line.strip()]
    page_rows.append(
        {
            "pdf_page": index,
            "characters": len(text),
            "first_nonempty_line": nonempty[0] if nonempty else "",
            "second_nonempty_line": nonempty[1] if len(nonempty) > 1 else "",
        }
    )
    bitmap = pdf[index - 1].render(scale=1.35)
    image = bitmap.to_pil().convert("RGB")
    target = PART / "pages" / f"page_{index:02d}.png"
    image.save(target, format="PNG", optimize=True)
    rendered.append(target)

(PART / "full_text_v1.txt").write_text("".join(text_sections), encoding="utf-8")
pd.DataFrame(page_rows).to_csv(PART / "page_index_v1.tsv", sep="\t", index=False, encoding="utf-8", lineterminator="\n")

thumb_width = 520
for sheet_index, start in enumerate(range(0, len(rendered), 4), start=1):
    items = []
    for path in rendered[start : start + 4]:
        image = Image.open(path).convert("RGB")
        height = round(image.height * thumb_width / image.width)
        image = image.resize((thumb_width, height), Image.Resampling.LANCZOS)
        canvas = Image.new("RGB", (thumb_width + 20, height + 50), "white")
        canvas.paste(image, (10, 40))
        draw = ImageDraw.Draw(canvas)
        draw.text((12, 12), path.stem.replace("_", " ").upper(), fill="black")
        items.append(canvas)
    sheet_height = max(item.height for item in items)
    sheet = Image.new("RGB", (sum(item.width for item in items), sheet_height), "#d9d9d9")
    x = 0
    for item in items:
        sheet.paste(ImageOps.pad(item, (item.width, sheet_height), color="white"), (x, 0))
        x += item.width
    sheet.save(PART / "contact_sheets" / f"pages_{start + 1:02d}_{min(start + 4, len(rendered)):02d}.png", format="PNG", optimize=True)

marker = (
    "status=PASS\n"
    f"source_sha256={sha256(SOURCE)}\n"
    f"source_bytes={SOURCE.stat().st_size}\n"
    f"pages={len(reader.pages)}\n"
    f"rendered_pages={len(rendered)}\n"
    f"text_characters={sum(row['characters'] for row in page_rows)}\n"
)
(PART / "PAPER_EXTRACTION_RENDER_PASS_v1").write_text(marker, encoding="utf-8")
os.replace(PART, FINAL)
print(marker, end="")
