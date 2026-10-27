from pathlib import Path
import requests,tarfile,struct,json,zlib,binascii,hashlib,traceback,time,os
P=Path(__file__).resolve().parents[1]
url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207493/suppl/GSE207493_scATAC-seq.tar.gz'
out=P/'downloads/RCC100';out.mkdir(parents=True,exist_ok=True)
manifest=[];start=time.time();compressed_total=0
def exact(f,n):
 b=f.read(n)
 if len(b)!=n:raise EOFError((n,len(b)))
 return b
try:
 with requests.get(url,stream=True,timeout=(15,120),proxies={'https':'http://127.0.0.1:17897','http':'http://127.0.0.1:17897'}) as response:
  response.raise_for_status()
  with tarfile.open(fileobj=response.raw,mode='r|gz') as tar:
   member=next(iter(tar));assert member.name=='scATAC-seq.zip'
   f=tar.extractfile(member)
   while True:
    head=exact(f,30);h=struct.unpack('<IHHHHHIIIHH',head)
    assert h[0]==0x04034b50,h
    flags,method,crc,csize,usize,namelen,extralen=h[2],h[3],h[6],h[7],h[8],h[9],h[10]
    name=exact(f,namelen).decode('utf-8');extra=exact(f,extralen)
    if not name.startswith('RCC100/'):print('NEXT_SAMPLE_STOP',name,flush=True);break
    assert not flags&8,'data descriptor unsupported'
    ep=0
    while ep+4<=len(extra):
     tag,n=struct.unpack('<HH',extra[ep:ep+4]);v=extra[ep+4:ep+4+n];ep+=4+n
     if tag==1:
      q=iter(struct.unpack('<'+'Q'*(len(v)//8),v))
      if usize==0xffffffff:usize=next(q)
      if csize==0xffffffff:csize=next(q)
    if name.endswith('/'):assert csize==0;continue
    assert compressed_total+csize<5_000_000_000,(name,csize)
    dest=(P/'downloads'/name).resolve()
    assert dest.is_relative_to(out.resolve()) and dest.name not in ['.','..']
    dest.parent.mkdir(parents=True,exist_ok=True)
    tmp=dest.with_name(dest.name+'.partial')
    dec=zlib.decompressobj(-15) if method==8 else None
    assert method in [0,8]
    print('EXTRACT_START',name,csize,usize,flush=True)
    left=csize;written=0;digest=hashlib.sha256();crc_obs=0;last=0
    with tmp.open('wb') as w:
     while left:
      b=exact(f,min(left,1024*1024));left-=len(b);compressed_total+=len(b)
      v=dec.decompress(b) if dec else b;w.write(v);written+=len(v);digest.update(v);crc_obs=binascii.crc32(v,crc_obs)
      if time.time()-last>60:
       (P/'DOWNLOAD_PROGRESS.json').write_text(json.dumps(dict(member=name,compressed_read=csize-left,compressed_total=compressed_total,member_size=csize,elapsed_seconds=time.time()-start,status='RUNNING'),indent=2));last=time.time()
     if dec:
      v=dec.flush();w.write(v);written+=len(v);digest.update(v);crc_obs=binascii.crc32(v,crc_obs);assert dec.eof
    assert written==usize,(name,written,usize)
    assert crc_obs&0xffffffff==crc,(name,crc_obs,crc)
    tmp.rename(dest)
    manifest.append(dict(name=name,path=str(dest),bytes=written,sha256=digest.hexdigest(),zip_crc32=hex(crc),crc_verified=True))
    (P/'metadata/RCC100_DOWNLOAD_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
    print('EXTRACT_DONE',name,written,flush=True)
 (P/'DOWNLOAD_STATUS.json').write_text(json.dumps(dict(status='RCC100_COMPLETE_CRC_VERIFIED',files=len(manifest),compressed_bytes=compressed_total,elapsed_seconds=time.time()-start),indent=2))
 (P/'DOWNLOAD_EXIT_STATUS').write_text('0\n')
except BaseException as e:
 (P/'DOWNLOAD_STATUS.json').write_text(json.dumps(dict(status='FAILED_PARTIAL_PRESERVED',error=repr(e),files=len(manifest)),indent=2))
 (P/'DOWNLOAD_EXIT_STATUS').write_text('1\n');traceback.print_exc();raise

