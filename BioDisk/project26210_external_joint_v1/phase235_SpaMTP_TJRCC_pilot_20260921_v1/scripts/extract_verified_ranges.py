from pathlib import Path
import struct,zlib,json,hashlib,datetime
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase235_SpaMTP_TJRCC_pilot_20260921_v1')
meta=json.loads((P/'metadata/ZIP_DIRECTORY.json').read_text());checks=[]
for sample,start in [('R114_T',108652358),('Y12_T',274095356)]:
 path=P/'downloads'/(sample+'_archive_range.bin');buf=path.read_bytes()
 for e in meta['entries']:
  if not e['name'].startswith(sample.replace('_','-')+'-'):continue
  pos=e['offset']-start;head=struct.unpack_from('<IHHHHHIIIHH',buf,pos)
  assert head[0]==0x04034b50 and head[3]==8 and not(head[2]&1)
  nl,el=head[-2:];name=buf[pos+30:pos+30+nl].decode()
  assert name==e['name'] and '/' not in name and '\\' not in name
  data=zlib.decompress(buf[pos+30+nl+el:pos+30+nl+el+e['compressed']],-15)
  crc=zlib.crc32(data)&0xffffffff
  assert len(data)==e['bytes'] and crc==e['crc'],name
  dest=P/'extracted'/name
  assert not dest.exists()
  dest.write_bytes(data)
  checks.append(dict(file=name,bytes=len(data),crc32=crc,crc_match=True,sha256=hashlib.sha256(data).hexdigest(),archive_offset=e['offset'],transport='local7897proxy range -> authorized sgs10 SCP'))
  print('VERIFIED',name,len(data),flush=True)
(P/'audit/EXTRACTED_FILE_INTEGRITY.json').write_text(json.dumps(checks,indent=2))
(P/'PILOT_DOWNLOAD_VERIFIED.json').write_text(json.dumps(dict(time=datetime.datetime.now().astimezone().isoformat(),files=len(checks),samples=2,all_entry_crc_pass=True,whole_archive_downloaded=False,whole_archive_MD5_not_claimed=True,proxy_port=7897),indent=2))
import importlib.util
for m in ['pyimzml','lxml','numpy','pandas','scipy','skimage']:
 print('MODULE',m,bool(importlib.util.find_spec(m)))
for f in sorted((P/'extracted').glob('*.imzML')):
 print('XML_HEADER',f.name,f.read_text()[:5000])
