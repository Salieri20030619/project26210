from pathlib import Path
import io
from pypdf import PdfReader,PdfWriter
from reportlab.pdfgen import canvas
r=Path.cwd();d=r/'phase206_all_mainfigures_A4_20260918_v10';mm=72/25.4
for n in [1,2,3]:
 f=d/'figures'/f'Main_Figure_{n}_A4_v10.pdf';rd=PdfReader(f);buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=(210*mm,297*mm));c.setFillColorRGB(1,1,1);c.rect(9*mm,278.5*mm,5*mm,6*mm,fill=1,stroke=0);c.setStrokeColorRGB(0,0,0);c.setLineWidth(.52*mm);c.setLineCap(0);p=c.beginPath();p.moveTo(10.2*mm,280*mm);p.lineTo(11.55*mm,283.2*mm);p.lineTo(12.9*mm,280*mm);c.drawPath(p);c.line(10.7*mm,281.15*mm,12.4*mm,281.15*mm);c.save();buf.seek(0);wr=PdfWriter();pg=wr.add_page(rd.pages[0]);pg.merge_page(PdfReader(buf).pages[0]);b=io.BytesIO();wr.write(b);f.write_bytes(b.getvalue())
print('VECTOR_LETTERS_READY')
