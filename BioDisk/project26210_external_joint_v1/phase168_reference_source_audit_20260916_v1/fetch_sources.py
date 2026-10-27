from pathlib import Path
import json,requests,zipfile,hashlib,gzip,concurrent.futures
O=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase168_reference_source_audit_20260916_v1')
z=json.loads((O/'zenodo_record_response.json').read_text())['files'][0]
urls={'zenodo_code.zip':z['links']['self'],'GSE210041_family.soft.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE210nnn/GSE210041/soft/GSE210041_family.soft.gz','published_supplement.pdf':'https://www.science.org/doi/suppl/10.1126/sciadv.aed0733/suppl_file/sciadv.aed0733_sm.pdf','author_repo_tree.json':'https://api.github.com/repos/qgggt61-commits/001/git/trees/v1?recursive=1'}
def fetch(item):
 name,url=item
 try:
  r=requests.get(url,timeout=(15,45));r.raise_for_status();(O/name).write_bytes(r.content)
  return dict(file=name,url=url,status=r.status_code,bytes=len(r.content),sha256=hashlib.sha256(r.content).hexdigest())
 except Exception as e:return dict(file=name,url=url,error=str(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:logs=list(pool.map(fetch,urls.items()))
(O/'FETCH_LOG.json').write_text(json.dumps(logs,indent=2));print(json.dumps(logs))
if (O/'zenodo_code.zip').exists():
 assert hashlib.md5((O/'zenodo_code.zip').read_bytes()).hexdigest()==z['checksum'].split(':')[1]
 with zipfile.ZipFile(O/'zenodo_code.zip') as f:
  names=f.namelist();(O/'ZIP_CONTENTS.json').write_text(json.dumps(names,indent=2))
  D=O/'zenodo_code';D.mkdir(exist_ok=True)
  for x in f.infolist():
   dest=(D/x.filename).resolve()
   assert dest.is_relative_to(D.resolve())
   f.extract(x,D)
 print('ZIP_FILES',json.dumps(names))
if (O/'GSE210041_family.soft.gz').exists():
 txt=gzip.decompress((O/'GSE210041_family.soft.gz').read_bytes()).decode()
 (O/'GSE210041_family.soft').write_text(txt)
 for block in txt.split('^SAMPLE = ')[1:]:
  print('\nSAMPLE',block.splitlines()[0])
  print('\n'.join(l for l in block.splitlines() if l.startswith(('!Sample_title','!Sample_characteristics','!Sample_source','!Sample_relation'))))
