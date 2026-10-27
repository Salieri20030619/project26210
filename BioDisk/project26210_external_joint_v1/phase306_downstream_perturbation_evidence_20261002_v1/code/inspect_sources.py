from pathlib import Path
import requests,json,hashlib,pandas as pd,concurrent.futures,xml.etree.ElementTree as ET,openpyxl
from bs4 import BeautifulSoup
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase306_downstream_perturbation_evidence_20261002_v1')
def session():
 s=requests.Session();s.trust_env=False;s.proxies={'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'};return s
def get_sample(acc):
 f=P/'metadata'/f'{acc}.xml'
 if not f.exists():
  r=session().get('https://www.ebi.ac.uk/ena/browser/api/xml/'+acc,timeout=45);r.raise_for_status();f.write_bytes(r.content)
 e=ET.fromstring(f.read_bytes()).find('SAMPLE');return {'sample_accession':acc,'alias':e.get('alias'),**{a.findtext('TAG'):a.findtext('VALUE') for a in e.findall('.//SAMPLE_ATTRIBUTE')}}
a=pd.read_csv(P/'metadata/ENA_PRJNA1223824.tsv',sep='\t')
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:out=list(ex.map(get_sample,a.sample_accession.unique()))
df=a.merge(pd.DataFrame(out),on='sample_accession');df.to_csv(P/'tables/RCC_ALL_RUN_METADATA.tsv',sep='\t',index=False)
print('SAMPLES',df[['run_accession','alias','age','collection_date']].to_string(index=False),flush=True)
for ident in [30255397]:
 r=session().get(f'https://api.figshare.com/v2/articles/{ident}',timeout=45);r.raise_for_status();j=r.json();(P/f'metadata/figshare_{ident}.json').write_text(json.dumps(j,indent=2))
 for f in j['files']:
  rr=session().get(f['download_url'],timeout=45);rr.raise_for_status();assert hashlib.md5(rr.content).hexdigest()==f['computed_md5'];(P/'downloads'/f['name']).write_bytes(rr.content);print('SUPPLEMENT',f['name'],len(rr.content),flush=True)
b=BeautifulSoup((P/'metadata/Trehan_article.html').read_text(),'html.parser')
for fig in b.find_all('figure'):
 cap=fig.get_text(' ',strip=True)
 if any(x in cap[:25] for x in ['Fig. 8','Fig. 9']):print('CAP',cap,flush=True)
w=openpyxl.load_workbook(P/'downloads/Trehan_SourceData.xlsx',read_only=True,data_only=True)
for sn in ['Figure 8','Figure 9','Supplementary Figures']:
 print('SHEET',sn,flush=True)
 for row in w[sn]:
  strings=[f'{c.coordinate}:{c.value}' for c in row if isinstance(c.value,str)]
  if strings:print(' | '.join(strings),flush=True)
