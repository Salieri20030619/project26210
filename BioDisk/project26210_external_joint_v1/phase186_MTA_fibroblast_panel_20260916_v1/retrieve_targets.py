from pathlib import Path
import subprocess,json,concurrent.futures,struct,zlib,zipfile,io,xml.etree.ElementTree as E
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase186_MTA_fibroblast_panel_20260916_v1')
panel=json.loads((r/'PRESET_PANEL.json').read_text())['genes']
def mapping(g):
 p=r/'sources'/f'hgnc_{g}.json'
 z=subprocess.run(['curl','-fLsS','--connect-timeout','10','--max-time','25','-H','Accept: application/json',f'https://rest.genenames.org/fetch/symbol/{g}','-o',str(p)],capture_output=True,text=True)
 try:
  a=json.loads(p.read_text())['response']['docs'];assert len(a)==1
  return {'symbol':g,'ensembl_id':a[0]['ensembl_gene_id'],'hgnc_id':a[0]['hgnc_id']}
 except Exception as e:return {'symbol':g,'error':str(e)}
out=list(concurrent.futures.ThreadPoolExecutor(max_workers=6).map(mapping,panel));(r/'GENE_MAPPING_HGNC.json').write_text(json.dumps(out,indent=2));print('MAPPING',json.dumps(out),flush=True)
rows=[x for x in json.loads((r/'ZIP_DIRECTORY.json').read_text()) if x['name'].endswith('.xlsx') and any(x['name'].startswith(y) for y in ['Figure 4/4C','Figure 4/4E','Figure 4/4F','Figure 5/5D','Supporting Figures/S20/','Supporting Figures/S21/','Supporting Figures/S24/'])]
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
out=list(concurrent.futures.ThreadPoolExecutor(max_workers=7).map(get,rows));(r/'QPCR_ENDPOINT_AUDIT.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
for x in out: print(json.dumps(x,ensure_ascii=False)[:2500])
