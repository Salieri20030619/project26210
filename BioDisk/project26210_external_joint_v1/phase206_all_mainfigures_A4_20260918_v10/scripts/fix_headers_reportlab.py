from pathlib import Path
import io
from pypdf import PdfReader,PdfWriter
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from matplotlib import font_manager
r=Path.cwd();d=r/'phase206_all_mainfigures_A4_20260918_v10';mm=72/25.4;pdfmetrics.registerFont(TTFont('Header',font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans',weight='bold'))))
for n,title in [(1,'Cell atlas'),(2,'Metabolic pathways'),(3,'Macrophage markers')]:
 p=d/'figures'/f'Main_Figure_{n}_A4_v10.pdf';reader=PdfReader(p);buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=(210*mm,297*mm));c.setFillColorRGB(1,1,1);c.rect(9*mm,(297-18.5)*mm,93*mm,6.5*mm,fill=1,stroke=0);c.setFillColorRGB(0,0,0);c.setFont('Header',10);c.drawString(10*mm,(297-17)*mm,'A');c.setFont('Header',7.5);c.drawString(15*mm,(297-17)*mm,title);c.save();buf.seek(0)
 writer=PdfWriter();page=writer.add_blank_page(210*mm,297*mm);page.merge_page(reader.pages[0]);page.merge_page(PdfReader(buf).pages[0]);tmp=p.with_name(p.stem+'_header.pdf')
 with tmp.open('wb') as f:writer.write(f)
 p.write_bytes(tmp.read_bytes());tmp.unlink()
print('HEADERS_REPAIRED')
