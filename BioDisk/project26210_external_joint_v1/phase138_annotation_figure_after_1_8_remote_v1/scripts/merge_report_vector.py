from pathlib import Path
import argparse
import json
import hashlib
from pypdf import PdfReader, PdfWriter, Transformation

parser = argparse.ArgumentParser()
parser.add_argument('--report', required=True)
parser.add_argument('--figure', required=True)
parser.add_argument('--output', required=True)
parser.add_argument('--page', type=int, default=4)
args = parser.parse_args()
report = PdfReader(args.report)
figure = PdfReader(args.figure)
placement = json.loads(Path(args.report).with_name('placement.json').read_text())['placements'][0]
target = placement['page'] - 1
assert len(figure.pages) == 1
assert '1.8' in report.pages[target-1].extract_text()
writer = PdfWriter()
for i, page in enumerate(report.pages):
    if i == target:
        fp = figure.pages[0]
        x, y = placement['x'], placement['y']
        scale = placement['width'] / float(fp.mediabox.width)
        assert y + float(fp.mediabox.height)*scale < float(page.mediabox.height)
        page.merge_transformed_page(fp, Transformation().scale(scale).translate(x, y), over=True)
    writer.add_page(page)
writer.add_metadata({'/Title': 'Project 26210 TYMP report with annotation Figure 1',
                     '/Author': 'Project 26210',
                     '/Subject': 'Annotation composite after section 1.8'})
output = Path(args.output)
output.parent.mkdir(parents=True, exist_ok=True)
with output.open('wb') as stream:
    writer.write(stream)
check = PdfReader(output)
assert len(check.pages) == len(report.pages)
assert len(check.pages[target].images) == len(figure.pages[0].images)
assert 'Single-cell identities' in check.pages[target].extract_text()
assert '2 分析结果' not in check.pages[target-1].extract_text()
assert '图1' in (check.pages[target].extract_text()+check.pages[target+1].extract_text())
print(json.dumps({'pages': len(check.pages), 'figure_page': args.page,
                  'vector_overlay': True, 'figure_scale': scale,
                  'figure_images_preserved': len(check.pages[target].images),
                  'sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
                  'bytes': output.stat().st_size}))
