from __future__ import annotations

import argparse
import base64
import io

from PIL import Image


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("image")
    parser.add_argument("--quality", type=int, default=88)
    args = parser.parse_args()
    with Image.open(args.image) as source:
        rgb = source.convert("RGB")
        buffer = io.BytesIO()
        rgb.save(buffer, format="JPEG", quality=args.quality, optimize=True)
    print(base64.b64encode(buffer.getvalue()).decode("ascii"))


if __name__ == "__main__":
    main()
