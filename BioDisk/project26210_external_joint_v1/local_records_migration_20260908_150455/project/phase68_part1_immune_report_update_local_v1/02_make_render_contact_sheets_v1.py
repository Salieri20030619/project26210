from pathlib import Path

from PIL import Image, ImageDraw


PHASE = Path(__file__).resolve().parent
SOURCE = PHASE / "rendered_local_v1"
OUTPUT = PHASE / "qc_local_v1"
PAGES_PER_SHEET = 4
THUMB_WIDTH = 760
LABEL_HEIGHT = 34
COLUMNS = 2


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    pages = sorted(SOURCE.glob("page-*.png"))
    if not pages:
        raise RuntimeError(f"No rendered pages found in {SOURCE}")
    with Image.open(pages[0]) as first:
        ratio = THUMB_WIDTH / first.width
        thumb_height = int(round(first.height * ratio))
    for sheet_index, start in enumerate(range(0, len(pages), PAGES_PER_SHEET), start=1):
        group = pages[start : start + PAGES_PER_SHEET]
        rows = (len(group) + COLUMNS - 1) // COLUMNS
        sheet = Image.new(
            "RGB",
            (COLUMNS * THUMB_WIDTH, rows * (thumb_height + LABEL_HEIGHT)),
            "white",
        )
        draw = ImageDraw.Draw(sheet)
        for index, page in enumerate(group):
            with Image.open(page) as image:
                thumb = image.convert("RGB").resize(
                    (THUMB_WIDTH, thumb_height), Image.Resampling.LANCZOS
                )
            x = (index % COLUMNS) * THUMB_WIDTH
            y = (index // COLUMNS) * (thumb_height + LABEL_HEIGHT)
            sheet.paste(thumb, (x, y + LABEL_HEIGHT))
            draw.text((x + 10, y + 8), page.stem, fill="#202225")
        sheet.save(OUTPUT / f"render_contact_sheet_{sheet_index:02d}_v1.png")


if __name__ == "__main__":
    main()
