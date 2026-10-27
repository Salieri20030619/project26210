from pathlib import Path
import fitz
r=Path.cwd();d=r/'phase206_all_mainfigures_A4_20260918_v10';mm=72/25.4
for n in [1,2,3]:
 doc=fitz.open(d/'figures'/f'Main_Figure_{n}_A4_v10.pdf');p=doc[0];p.get_pixmap(matrix=fitz.Matrix(4,4),clip=fitz.Rect(9*mm,12*mm,58*mm,20*mm)).save(str(d/'qc'/f'saved_header{n}.png'));p.get_pixmap(matrix=fitz.Matrix(1.8,1.8)).save(str(d/'qc'/f'Final_Fig{n}.png'))
print('SAVED_PDF_RENDER_COMPLETE')
