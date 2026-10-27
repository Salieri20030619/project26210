from pathlib import Path
import subprocess,concurrent.futures,struct,json,xml.etree.ElementTree as E
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase186_MTA_fibroblast_panel_20260916_v1')
jobs=[('zip_tail.bin','https://zenodo.org/api/records/8130216/files/PNAS_Data.zip/content',['-H','Range: bytes=30392672-30523743']),('PMC8976074.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8976074/fullTextXML',[])]
def get(j):
 n,u,extra=j
 z=subprocess.run(['curl','-fLsS','--connect-timeout','10','--max-time','45',*extra,u,'-o',str(r/'sources'/n)],capture_output=True,text=True)
 print(n,z.returncode,z.stderr,flush=True)
 if z.returncode==0 and n.endswith('xml'):
  x=E.parse(r/'sources'/n)
  for e in x.iter('p'):
   t=' '.join(e.itertext())
   if any(k in t.lower() for k in ['wnt5a','chemerin','tgfb','tgf-','co-cultur']):print(t[:4000],flush=True)
 if z.returncode==0 and n.endswith('bin'):
  b=(r/'sources'/n).read_bytes(); print('TAIL_BYTES',len(b))
  k=b.rfind(b'PK\x05\x06')
  if k>=0: print('EOCD',struct.unpack_from('<4s4H2LH',b,k))
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(get,jobs))
