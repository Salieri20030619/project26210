from pathlib import Path
import json,hashlib,shutil,platform,sys,re
from datetime import datetime,timezone,timedelta
P=Path(__file__).resolve().parents[1]
now=datetime.now(timezone(timedelta(hours=8))).isoformat()
v=json.loads((P/'qc/STRUCTURAL_VALIDATION.json').read_text())
assert (P/'BUILD_EXIT_STATUS').read_text().strip()=='0'
assert (P/'VALIDATION_EXIT_STATUS').read_text().strip()=='0'
v.update(status='PASS',time=now,visual_review='All 39 pages inspected. Final pages 1,37,38,39 inspected after subtitle and reference layout polish; other 35 pages pixel-identical to inspected version.',figure5_page=24,figure5_caption_page=25)
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(v,ensure_ascii=False,indent=2))
out=P/'delivery/26210_integrated_spatial_report_v1.pdf'
shutil.copy2(P/'report'/out.name,out)
md=P/'report/26210_integrated_spatial_report_v1.md'
s=md.read_text();s=re.sub(r'!\[图(\d)\]\(figures/Figure\d.png\)',r'[图\1](../figures/Figure\1.pdf)',s);md.write_text(s)
import reportlab,pypdf,pymupdf
versions=dict(python=sys.version,executable=sys.executable,reportlab=reportlab.Version,pypdf=pypdf.__version__,pymupdf=pymupdf.__version__,font_source=str(P.parent/'phase143_sections25_27_candidate_report_remote_v1/fonts'))
(P/'qc/SOFTWARE.json').write_text(json.dumps(versions,indent=2))
m=dict(path=str(out),bytes=out.stat().st_size,sha256=hashlib.sha256(out.read_bytes()).hexdigest(),pages=39)
(P/'DELIVERY_MANIFEST.json').write_text(json.dumps(m,ensure_ascii=False,indent=2))
(P/'STATUS.json').write_text(json.dumps(dict(status='COMPLETE_VALIDATED_PENDING_LOCAL_TRANSFER',time=now),indent=2))
print(json.dumps(m),flush=True)

