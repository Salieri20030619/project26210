from pathlib import Path

from PIL import Image, ImageDraw


PHASE = Path(__file__).resolve().parent
SOURCE = PHASE / "rendered_local_v1"
OUTPUT = PHASE / "qc_local_v1" / "render_contact_sheet_v1.png"
pages = sorted(SOURCE.glob("page-*.png"))
thumb_width = 300
label_height = 28
columns = 4
rows = (len(pages) + columns - 1) // columns

with Image.open(pages[0]) as first:
    ratio = thumb_width / first.width
    thumb_height = int(round(first.height * ratio))

sheet = Image.new("RGB", (columns * thumb_width, rows * (thumb_height + label_height)), "white")
draw = ImageDraw.Draw(sheet)
for index, page in enumerate(pages):
    with Image.open(page) as image:
        thumb = image.convert("RGB").resize((thumb_width, thumb_height), Image.Resampling.LANCZOS)
    x = (index % columns) * thumb_width
    y = (index // columns) * (thumb_height + label_height)
    sheet.paste(thumb, (x, y + label_height))
    draw.text((x + 8, y + 6), page.stem, fill="#202225")

sheet.save(OUTPUT)
