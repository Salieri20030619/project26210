from pathlib import Path
import json,hashlib,subprocess
from PIL import Image,ImageOps,ImageDraw
from pypdf import PdfReader
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase192_sections28_212_figures_20260917_v1'
qc=[];thumbs=[]
for f in sorted((d/'figures').glob('*.pdf')):
 p=PdfReader(f);assert len(p.pages)==1
 png=f.with_suffix('.png');im=Image.open(png)
 assert min(im.info.get('dpi',(0,0)))>599
 subprocess.run(['pdftoppm','-r','110','-singlefile','-jpeg',str(f),str(d/'qc'/('pdf_'+f.stem))],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
 preview=Image.open(d/'qc'/('pdf_'+f.stem+'.jpg')).convert('RGB')
 tile=Image.new('RGB',(700,710),'white');preview.thumbnail((690,680));tile.paste(preview,((700-preview.width)//2,5))
 ImageDraw.Draw(tile).text((15,690),f.stem,fill='#555555');thumbs.append(tile)
 qc.append({'file':f.name,'pages':1,'width_points':float(p.pages[0].mediabox.width),'height_points':float(p.pages[0].mediabox.height),'png_size':im.size,'dpi':im.info.get('dpi')})
sheet=Image.new('RGB',(1400,2840),'white')
for i,t in enumerate(thumbs):sheet.paste(t,((i%2)*700,(i//2)*710))
sheet.save(d/'qc/PDF_CONTACT_SHEET.jpg',quality=90)
(d/'qc/FILE_VALIDATION.json').write_text(json.dumps(qc,indent=2))
sources=[r/'part2_rebuild_v1/08_regulons/input_v1/CELL_METADATA.tsv.gz',r/'part2_rebuild_v1/08_regulons/input_v1/GENES_BY_CELLS_RAW.npz',r/'part2_rebuild_v1/03_macrophage_clustering/run_v4_exploratory/UMAP_COORDINATES.tsv.gz',r/'part2_rebuild_v1/07_trajectory/monocle3_v3/ALL_CELL_PSEUDOTIME.tsv.gz',r/'part2_rebuild_v1/07_trajectory/dynamics_v2/CELL_EXPRESSION_PROGRAM_SCORES.tsv.gz']
sources+=list((d/'scripts').glob('*'))+list((d/'tables').glob('*'))
(d/'qc/SOURCE_AND_CODE_SHA256.json').write_text(json.dumps({str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in sources},indent=2))
print('Validated eight PDFs and 600-dpi PNGs')
