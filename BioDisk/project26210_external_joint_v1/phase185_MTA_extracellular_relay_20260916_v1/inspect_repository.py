from pathlib import Path
import urllib.request, json, zipfile, hashlib, xml.etree.ElementTree as E
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase185_MTA_extracellular_relay_20260916_v1')
u='https://zenodo.org/records/8130216/files/PNAS_Data.zip?download=1'
try:
 with urllib.request.urlopen(u,timeout=45) as a,open(r/'sources/PNAS_Data.zip','wb') as b:
  while True:
   c=a.read(1024*1024)
   if not c:break
   b.write(c)
 p=r/'sources/PNAS_Data.zip'; h=hashlib.md5(p.read_bytes()).hexdigest()
 assert h=='2a4fc1c9809ff58ea712ec4897b1998b',h
 z=zipfile.ZipFile(p); assert z.testzip() is None
 names=z.namelist(); (r/'ZENODO_FILE_INVENTORY.json').write_text(json.dumps(names,indent=2))
 print('ZIP_OK',p.stat().st_size,h)
 print('\n'.join(names)[:12000])
except Exception as e:print('ZIP_FAILED',str(e))
x=E.parse(r/'sources/PMC10756296.xml')
for e in x.iter('media'):print('MEDIA',e.attrib)
