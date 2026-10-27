from pathlib import Path
import requests,json,hashlib,os,datetime,concurrent.futures
P=Path(__file__).resolve().parents[1]
def fetch(item):
 name,url=item;s=requests.Session();s.trust_env=False;s.proxies={'https':'http://127.0.0.1:7897'}
 f=P/'downloads'/name
 for attempt in range(3):
  try:
   r=s.get(url,stream=True,timeout=(25,80));r.raise_for_status();size=int(r.headers.get('Content-Length',0));assert 0<size<2_000_000_000
   if not f.exists():
    t=f.with_suffix(f.suffix+'.partial')
    with t.open('wb') as out:
     for b in r.iter_content(4*1024*1024):out.write(b)
    assert t.stat().st_size==size;os.replace(t,f)
   r.close();assert f.stat().st_size==size
   rec={'name':name,'url':url,'bytes':size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'route':'proxy','http_etag':r.headers.get('ETag')}
   print('FETCHED',name,size,flush=True);return rec
  except Exception as e:
   print('FETCH_RETRY',name,attempt,type(e).__name__,flush=True)
   if attempt==2:raise
base='https://cf.10xgenomics.com/samples/xenium/4.0.0/Xenium_V1_Human_Clear_Cell_Renal_Cell_Carcinoma_FFPE_Protein/Xenium_V1_Human_Clear_Cell_Renal_Cell_Carcinoma_FFPE_Protein_'
items=[('tenx_'+x,base+x) for x in ['cells.csv.gz','cell_feature_matrix.h5','cell_groups.csv']]
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool: records=list(pool.map(fetch,items))
(P/'qc/TENX_DOWNLOAD_VALIDATION.json').write_text(json.dumps({'status':'PASS','time':datetime.datetime.now().astimezone().isoformat(),'files':records},indent=2))
import pandas as pd,h5py
for name in ['tenx_cells.csv.gz','tenx_cell_groups.csv']:
 df=pd.read_csv(P/'downloads'/name);print('SCHEMA',name,df.shape,list(df.columns),flush=True);print(df.head(4).to_string(index=False),flush=True)
 for c in df:
  if df[c].nunique()<100:print('VALUES',c,df[c].value_counts().to_dict(),flush=True)
with h5py.File(P/'downloads/tenx_cell_feature_matrix.h5') as f:
 print('H5_FEATURES',list(f['matrix/features']),flush=True)
 for c in ['feature_type','name']:
  print('FEATURE',c,[x.decode() for x in f['matrix/features'][c][:]],flush=True)
