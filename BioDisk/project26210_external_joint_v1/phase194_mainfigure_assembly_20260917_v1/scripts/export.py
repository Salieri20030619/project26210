from pathlib import Path
import subprocess,json,hashlib,shutil,datetime
from pypdf import PdfReader
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase194_mainfigure_assembly_20260917_v1'
# Export from final vector PDFs.
for f in sorted((d/'figures').glob('*.pdf')):
 assert len(PdfReader(f).pages)==1
 subprocess.run(['pdftoppm','-r','600','-singlefile','-png',str(f),str(f.with_suffix(''))],check=True)
 subprocess.run(['pdftocairo','-svg',str(f),str(f.with_suffix('.svg'))],check=True)
 print('EXPORTED',f.stem,flush=True)
