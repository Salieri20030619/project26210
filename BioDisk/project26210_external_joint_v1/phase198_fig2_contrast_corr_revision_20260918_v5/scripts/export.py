from pathlib import Path
import subprocess
d=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase198_fig2_contrast_corr_revision_20260918_v5')
for p in sorted((d/'figures').glob('*.pdf')):
 subprocess.run(['pdftoppm','-scale-to','1800','-singlefile','-jpeg',str(p),str(d/'qc'/p.stem)],check=True)
 subprocess.run(['pdftoppm','-r','600','-singlefile','-png',str(p),str(p.with_suffix(''))],check=True)
 subprocess.run(['pdftocairo','-svg',str(p),str(p.with_suffix('.svg'))],check=True)
 print('EXPORTED',p.name,flush=True)
