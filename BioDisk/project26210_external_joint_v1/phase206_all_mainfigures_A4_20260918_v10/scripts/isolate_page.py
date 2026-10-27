from pathlib import Path
import fitz
from matplotlib import font_manager
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
from PIL import Image
r=Path.cwd();d=r/'phase206_all_mainfigures_A4_20260918_v10';mm=72/25.4;hf=font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans',weight='bold'))
for n in [1,2,3]:
 p=d/'figures'/f'Main_Figure_{n}_A4_v10.pdf';old=fitz.open(p);new=fitz.open();pg=new.new_page(width=210*mm,height=297*mm);pg.show_pdf_page(pg.rect,old,0)
 pg.draw_rect(fitz.Rect(9.5*mm,13.3*mm,14*mm,18*mm),color=None,fill=(1,1,1));shape=pg.new_shape()
 for poly in TextPath((0,0),'A',size=10,prop=FontProperties(fname=hf)).to_polygons():
  shape.draw_polyline([fitz.Point(10*mm+xx,17*mm-yy) for xx,yy in poly])
 shape.finish(color=None,fill=(0,0,0),even_odd=True);shape.commit()
 tmp=p.with_name(p.stem+'_fixed.pdf');new.save(tmp,garbage=4,deflate=True);old.close();p.write_bytes(tmp.read_bytes());tmp.unlink()
 pg.get_pixmap(matrix=fitz.Matrix(1.8,1.8)).save(str(d/'qc'/f'Fig{n}_preview.png'))
 pg.get_pixmap(matrix=fitz.Matrix(4,4),clip=fitz.Rect(9*mm,12*mm,58*mm,20*mm)).save(str(d/'qc'/f'header{n}.png'))
print('GRAPHICS_STATE_ISOLATED')
