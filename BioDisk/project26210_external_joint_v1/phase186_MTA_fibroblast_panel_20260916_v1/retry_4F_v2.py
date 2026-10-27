import concurrent.futures
from pathlib import Path
import json,subprocess,struct,zlib,zipfile,io,xml.etree.ElementTree as E
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase186_MTA_fibroblast_panel_20260916_v1')
rows=[x for x in json.loads((r/'ZIP_DIRECTORY.json').read_text()) if x['name'].endswith('.xlsx') and any(x['name'].startswith(y) for y in ['Figure 4/4C','Figure 4/4E','Figure 4/4F','Figure 5/5D','Supporting Figures/S20/','Supporting Figures/S21/','Supporting Figures/S24/'])]
rows=[x for x in rows if x['name']=='Figure 4/4F/Figure 4F.xlsx']
d=r/'sources/qPCR';d.mkdir(exist_ok=True)
def get(x):
 st=x['offset'];end=st+30+len(x['name'].encode())+x['compressed']+1000
 p=d/(x['name'].split('/')[-1]+'.range')
 z=subprocess.run(['curl','-fLsS','--connect-timeout','10','--max-time','35','-H',f'Range: bytes={st}-{end}','https://zenodo.org/api/records/8130216/files/PNAS_Data.zip/content','-o',str(p)],capture_output=True,text=True)
 try:
  b=p.read_bytes();a=struct.unpack_from('<4s5H3L2H',b,0);assert a[0]==b'PK\x03\x04'
  off=30+a[-2]+a[-1];data=zlib.decompress(b[off:off+x['compressed']],-15)
  assert len(data)==x['uncompressed'] and zlib.crc32(data)==x['crc32']
  dest=d/x['name'].split('/')[-1];dest.write_bytes(data)
  zz=zipfile.ZipFile(io.BytesIO(data)); texts=[]
  for n in zz.namelist():
   if n.endswith('.xml'):
    xml=E.fromstring(zz.read(n))
    texts.extend(e.text for e in xml.iter() if e.tag.endswith('}t') and e.text)
  (d/(dest.stem+'_strings.json')).write_text(json.dumps(texts,ensure_ascii=False,indent=2))
  return {'file':x['name'],'crc_pass':True,'strings':sorted(set(texts))}
 except Exception as e:return {'file':x['name'],'error':str(e),'curl_error':z.stderr}
out=list(concurrent.futures.ThreadPoolExecutor(max_workers=7).map(get,rows));(r/'QPCR_ENDPOINT_RETRY.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
for x in out: print(json.dumps(x,ensure_ascii=False)[:2500])
