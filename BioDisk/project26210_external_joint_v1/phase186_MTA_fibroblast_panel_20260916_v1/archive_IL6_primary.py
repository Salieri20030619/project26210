from pathlib import Path
import subprocess,xml.etree.ElementTree as E
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase186_MTA_fibroblast_panel_20260916_v1')
z=subprocess.run(['curl','-fLsS','--connect-timeout','10','--max-time','35','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11271635/fullTextXML','-o',str(r/'sources/PMC11271635.xml')],capture_output=True,text=True)
print('XML',z.returncode,z.stderr)
if z.returncode==0:
 x=E.parse(r/'sources/PMC11271635.xml')
 for e in x.iter('p'):
  t=' '.join(e.itertext())
  if any(w in t for w in ['IL-6 neutralization','recombinant IL-6','monocyte differentiation','GSE']) and len(t)<10000:print(t[:2800])
