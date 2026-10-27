from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import win32com.client


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("docx")
    parser.add_argument("output_dir")
    parser.add_argument("--first-page", type=int)
    parser.add_argument("--last-page", type=int)
    parser.add_argument("--dpi", type=int, default=144)
    args = parser.parse_args()

    docx = Path(args.docx).resolve()
    out = Path(args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    pdf = out / f"{docx.stem}.pdf"

    word = None
    document = None
    try:
        word = win32com.client.DispatchEx("Word.Application")
        word.Visible = False
        word.DisplayAlerts = 0
        document = word.Documents.Open(str(docx), ReadOnly=True, AddToRecentFiles=False)
        document.Repaginate()
        pages = int(document.ComputeStatistics(2))  # wdStatisticPages
        document.ExportAsFixedFormat(str(pdf), 17, OpenAfterExport=False)
    finally:
        if document is not None:
            document.Close(False)
        if word is not None:
            word.Quit()

    gs = Path(r"D:\aaa\gs10.07.1\bin\gswin64c.exe")
    command = [
        str(gs),
        "-q",
        "-dSAFER",
        "-dBATCH",
        "-dNOPAUSE",
        "-sDEVICE=png16m",
        f"-r{args.dpi}",
        "-dTextAlphaBits=4",
        "-dGraphicsAlphaBits=4",
    ]
    if args.first_page:
        command.append(f"-dFirstPage={args.first_page}")
    if args.last_page:
        command.append(f"-dLastPage={args.last_page}")
    command.extend([f"-sOutputFile={out / 'page-%03d.png'}", str(pdf)])
    subprocess.run(command, check=True)

    pngs = sorted(out.glob("page-*.png"))
    status = {
        "status": "PASS" if pdf.exists() and pngs else "FAIL",
        "docx": str(docx),
        "pdf": str(pdf),
        "document_pages": pages,
        "rendered_pngs": len(pngs),
        "first_page": args.first_page,
        "last_page": args.last_page,
        "dpi": args.dpi,
    }
    (out / "render_status.json").write_text(
        json.dumps(status, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(status, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
