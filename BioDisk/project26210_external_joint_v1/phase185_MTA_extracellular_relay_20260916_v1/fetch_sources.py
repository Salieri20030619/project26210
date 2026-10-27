from pathlib import Path
import urllib.request, concurrent.futures, xml.etree.ElementTree as E, json
root=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase185_MTA_extracellular_relay_20260916_v1')
(root/'sources').mkdir(parents=True,exist_ok=True)
def fetch(pmc):
 u=f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML'
 try:
  b=urllib.request.urlopen(u,timeout=45).read(); x=E.fromstring(b)
  (root/'sources'/f'{pmc}.xml').write_bytes(b)
  txt='\n'.join(' '.join(e.itertext()) for e in x.iter() if e.tag in ['p','title','article-title','table-wrap','sec'] and e.tag!='sec')
  (root/'sources'/f'{pmc}.txt').write_text(txt)
  print(pmc,'OK',len(b))
  for line in txt.splitlines():
   if any(w in line.lower() for w in ['conditioned','hypoxanth','spp1','osteopontin','data availability','data and code','gse','methylthioadenosine','polarization']): print(pmc,line[:4500])
 except Exception as e: print(pmc,'FAILED',str(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:list(ex.map(fetch,['PMC10756296','PMC11973769','PMC7072758']))
