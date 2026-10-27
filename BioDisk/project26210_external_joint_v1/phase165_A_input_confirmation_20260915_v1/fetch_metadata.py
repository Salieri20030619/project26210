from pathlib import Path
import requests,gzip,json,hashlib
O=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase165_A_input_confirmation_20260915_v1');(O/'metadata').mkdir(parents=True,exist_ok=True)
u='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE175nnn/GSE175540/soft/GSE175540_family.soft.gz'
r=requests.get(u,timeout=45);r.raise_for_status();(O/'metadata/GSE175540_family.soft.gz').write_bytes(r.content)
s=gzip.decompress(r.content).decode('utf-8','replace')
(O/'metadata/GSE175540_family.soft').write_text(s,encoding='utf-8')
samples=[]
for block in s.split('^SAMPLE = ')[1:]:
 lines=block.splitlines();gsm=lines[0].strip()
 if not gsm.startswith('GSM59240'):continue
 vals={}
 for line in lines[1:]:
  if line.startswith('!Sample_') and ' = ' in line:
   key,val=line.split(' = ',1);vals.setdefault(key,[]).append(val)
 samples.append(dict(gsm=gsm,title=vals.get('!Sample_title'),source=vals.get('!Sample_source_name_ch1'),characteristics=vals.get('!Sample_characteristics_ch1'),description=vals.get('!Sample_description')))
(O/'metadata/SAMPLE_METADATA_EXTRACT.json').write_text(json.dumps(samples,ensure_ascii=False,indent=2))
print(json.dumps(samples,ensure_ascii=False,indent=2));print('SHA256',hashlib.sha256(r.content).hexdigest())

