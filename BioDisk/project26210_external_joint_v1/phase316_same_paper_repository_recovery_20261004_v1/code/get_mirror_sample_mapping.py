from pathlib import Path
import requests,json,struct,zlib,hashlib,subprocess
P=Path(__file__).resolve().parents[1];j=json.loads((P/'resources/FIGURE34_ZIP_INDEX.json').read_text());s=requests.Session();s.trust_env=False
with s.get(j['url'],headers={'Range':'bytes=0-109598'},timeout=(15,40),stream=True) as r:
 assert r.status_code==206
 data=r.raw.read(109600);assert len(data)==109599
(P/'resources/FIGURE34_code_range.bin').write_bytes(data)
rows=[]
for x in j['members']:
 if x['offset']>=109599 or not x['name'].endswith('.R'):continue
 off=x['offset'];h=struct.unpack_from('<4s5H3I2H',data,off);assert h[0]==b'PK\x03\x04';start=off+30+h[-2]+h[-1];raw=zlib.decompress(data[start:start+x['compressed']],-15)
 assert len(raw)==x['bytes'] and zlib.crc32(raw)==x['CRC']
 name=x['name'].replace('/','__');(P/'resources'/name).write_bytes(raw)
 rows.append({'member':x['name'],'CRC':x['CRC'],'sha256':hashlib.sha256(raw).hexdigest()})
(P/'qc/ZIP_CODE_MEMBER_VALIDATION.json').write_text(json.dumps(rows,indent=2))
for f in (P/'resources').glob('Figure34*R'):
 lines=f.read_text(errors='replace').splitlines()
 for i,l in enumerate(lines):
  if any(k in l for k in ['readRDS','P7','P6','CRC8','Human_CRC','ST_data','HRA000979']):print(f.name,i+1,l[:200])
