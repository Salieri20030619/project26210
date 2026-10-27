from pathlib import Path
import json,subprocess,concurrent.futures,xml.etree.ElementTree as ET,datetime
P=Path(__file__).resolve().parents[1]
urls=[
('ebi_info.json','https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-12767/info'),
('sdrf_fire.txt','https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/767/E-MTAB-12767/Files/E-MTAB-12767.sdrf.txt'),
('sdrf_arrayexpress.txt','https://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/experiment/MTAB/E-MTAB-12767/E-MTAB-12767.sdrf.txt'),
('sdrf_public.txt','https://www.ebi.ac.uk/biostudies/files/E-MTAB-12767/E-MTAB-12767.sdrf.txt')]
def fetch(x):
 name,url=x
 r=subprocess.run(['curl','-q','-fLsS','--proxy','http://127.0.0.1:17899','--noproxy','','--connect-timeout','12','--max-time','25',url],capture_output=True,timeout=30)
 if r.returncode==0:(P/'sources'/name).write_bytes(r.stdout)
 return dict(name=name,url=url,exit=r.returncode,bytes=len(r.stdout),error=r.stderr.decode()[-180:])
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:results=list(ex.map(fetch,urls))
(P/'SDRF_ACCESS.json').write_text(json.dumps(results,indent=2));print(json.dumps(results),flush=True)
for name,_ in urls:
 f=P/'sources'/name
 if f.exists():print(name,f.read_text()[:1000])
root=ET.parse(P/'sources/mRCC_paper.xml')
for x in root.iter('ref'):
 if x.get('id')=='R48':print('REFERENCE48',ET.tostring(x,encoding='unicode'))

