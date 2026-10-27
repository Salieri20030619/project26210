from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import requests,json,gzip,hashlib,shutil,traceback,os
P=Path(__file__).resolve().parents[1];D=P/'reference';D.mkdir(exist_ok=True)
urls={
'GRCh38.genome.fa.gz':'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_21/GRCh38.genome.fa.gz',
'gencode.v38.basic.annotation.gtf.gz':'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.basic.annotation.gtf.gz'
}
def run(k,u):
 dest=D/k;partial=D/(k+'.partial');sha=hashlib.sha256()
 with requests.get(u,stream=True,timeout=(15,120),proxies={'https':'http://127.0.0.1:17897'}) as r:
  r.raise_for_status();size=int(r.headers.get('Content-Length',0));n=0
  with partial.open('wb') as f:
   for b in r.iter_content(1024*1024):
    f.write(b);sha.update(b);n+=len(b)
  assert not size or n==size,(k,n,size)
 partial.rename(dest)
 out=D/k.removesuffix('.gz')
 with gzip.open(dest,'rb') as a,out.open('wb') as b:shutil.copyfileobj(a,b,length=1024*1024)
 print('REFERENCE_READY',k,n,out.stat().st_size,flush=True)
 return dict(file=str(dest),url=u,bytes=n,sha256=sha.hexdigest(),gzip_crc_verified=True,expanded=str(out))
try:
 with ThreadPoolExecutor(max_workers=2) as pool:records=list(pool.map(lambda x:run(*x),urls.items()))
 (P/'metadata/REFERENCE_DOWNLOADS.json').write_text(json.dumps(records,indent=2));(P/'REFERENCE_EXIT_STATUS').write_text('0\n')
except BaseException:
 (P/'REFERENCE_EXIT_STATUS').write_text('1\n');traceback.print_exc();raise

