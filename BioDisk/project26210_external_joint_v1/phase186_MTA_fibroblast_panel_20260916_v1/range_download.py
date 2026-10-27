from pathlib import Path
import json,subprocess,concurrent.futures,gzip
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase186_MTA_fibroblast_panel_20260916_v1')
m=json.loads((r/'sources/zenodo_metadata.json').read_text())
print('ZENODO_FILES',[(f['key'],f['size'],f['links']) for f in m['files']],flush=True)
url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE249nnn/GSE249082/suppl/GSE249082_AC-59_RawCounts.txt.gz'
size=5577063; chunk=262144
d=r/'sources/count_chunks';d.mkdir(exist_ok=True)
def get(i):
 start=i*chunk;end=min(size-1,start+chunk-1); p=d/f'{i:03}.bin';h=d/f'{i:03}.headers'
 z=subprocess.run(['curl','-fLsS','--connect-timeout','10','--max-time','65','-H',f'Range: bytes={start}-{end}','-D',str(h),url,'-o',str(p)],capture_output=True,text=True)
 ok=z.returncode==0 and p.exists() and p.stat().st_size==end-start+1 and f'content-range: bytes {start}-{end}/{size}' in h.read_text().lower()
 return {'i':i,'start':start,'end':end,'ok':ok,'error':z.stderr}
out=list(concurrent.futures.ThreadPoolExecutor(max_workers=10).map(get,range((size+chunk-1)//chunk)))
(r/'RANGE_RETRIEVAL.json').write_text(json.dumps(out,indent=2))
print('CHUNKS',sum(x['ok'] for x in out),len(out),flush=True)
if all(x['ok'] for x in out):
 p=r/'sources/counts_complete.txt.gz';p.write_bytes(b''.join((d/f'{x["i"]:03}.bin').read_bytes() for x in out))
 data=gzip.decompress(p.read_bytes()).decode()
 print('HEADER',data.splitlines()[0]);print('FIRST',data.splitlines()[1][:800]);print('ROWS',len(data.splitlines()))
