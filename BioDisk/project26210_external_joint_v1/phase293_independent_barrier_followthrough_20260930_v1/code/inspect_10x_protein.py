from pathlib import Path
import subprocess,re,json,io,zipfile,hashlib
P=Path(__file__).resolve().parents[1]
urls=json.loads((P/'resources/BIOIMAGE_URLS.json').read_text())
u=[x for x in urls if x.endswith('_outs.zip')][0]
head=subprocess.run(['curl','-L','-sS','-I','--connect-timeout','15','--max-time','55',u],capture_output=True,timeout=65)
(P/'resources/10x_archive_HEAD_protein.txt').write_bytes(head.stdout+head.stderr)
assert head.returncode==0,head.stderr
h=head.stdout.decode(errors='replace');sizes=re.findall(r'(?im)^content-length:\s*(\d+)',h);assert sizes,h
size=int(sizes[-1]);assert size>0
class Remote(io.RawIOBase):
 def __init__(self):self.pos=0;self.calls=0
 def seekable(self):return True
 def readable(self):return True
 def tell(self):return self.pos
 def seek(self,pos,whence=0):
  self.pos=pos if whence==0 else self.pos+pos if whence==1 else size+pos
  return self.pos
 def read(self,n=-1):
  if n<0:n=size-self.pos
  if n==0:return b''
  assert 0<=self.pos<size and n<=10485760,(self.pos,n,size)
  end=min(size-1,self.pos+n-1)
  q=subprocess.run(['curl','-L','--fail','-sS','--connect-timeout','15','--max-time','90','--max-filesize','10485760','-r',f'{self.pos}-{end}',u],capture_output=True,timeout=100)
  assert q.returncode==0,q.stderr
  b=q.stdout;assert len(b)==end-self.pos+1,(len(b),self.pos,end)
  self.calls+=1
  (P/'resources'/f'10x_zip_protein_range_{self.calls}.bin').write_bytes(b)
  print('RANGE',self.pos,end,len(b),flush=True);self.pos=end+1;return b
r=Remote()
with zipfile.ZipFile(r) as z:
 entries=[dict(name=x.filename,bytes=x.file_size,compressed=x.compress_size,offset=x.header_offset) for x in z.infolist()]
 (P/'resources/10x_ARCHIVE_MEMBERS_protein.json').write_text(json.dumps(entries,indent=2))
 chosen=[x for x in z.infolist() if x.filename.endswith('protein_panel.json')]
 for x in chosen:
  assert x.file_size<5000000
  b=z.read(x);name=Path(x.filename).name
  (P/'inputs'/('10x_'+name)).write_bytes(b)
  print('EXTRACTED',name,len(b),hashlib.sha256(b).hexdigest(),flush=True)
print(json.dumps(dict(status='PANEL_ACQUIRED',archive_bytes=size,range_calls=r.calls,selected=len(chosen))))
