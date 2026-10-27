from pathlib import Path
import requests,json,hashlib,datetime,concurrent.futures
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
P=R/'phase306_downstream_perturbation_evidence_20261002_v1'
for d in ['code','metadata','downloads','resources','tables','qc','logs','delivery']:(P/d).mkdir(parents=True,exist_ok=True)
urls={
'RCC_article.html':'https://pmc.ncbi.nlm.nih.gov/articles/PMC12485372/',
'OMIX010399.html':'https://ngdc.cncb.ac.cn/omix/release/OMIX010399',
'OMIX010400.html':'https://ngdc.cncb.ac.cn/omix/release/OMIX010400',
'Trehan_article.html':'https://www.nature.com/articles/s41467-025-59529-0',
'Trehan_SourceData.xlsx':'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-025-59529-0/MediaObjects/41467_2025_59529_MOESM4_ESM.xlsx',
'ENA_PRJNA1223824.tsv':'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1223824&result=read_run&fields=run_accession,sample_accession,experiment_accession,sample_title,experiment_title,fastq_bytes,fastq_ftp&format=tsv',
}
def session():
 s=requests.Session();s.trust_env=False;s.proxies={'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'};return s
def fetch(k,u):
 try:
  r=session().get(u,timeout=55);dest=P/('downloads' if k.endswith('.xlsx') else 'metadata')/k;dest.write_bytes(r.content)
  out={'file':str(dest.relative_to(P)),'url':u,'status':r.status_code,'bytes':len(r.content),'sha256':hashlib.sha256(r.content).hexdigest()};print(json.dumps(out),flush=True);return out
 except Exception as e:return {'file':k,'url':u,'error':str(e)}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:records=list(ex.map(lambda x:fetch(*x),urls.items()))
try:
 r=session().post('https://api.figshare.com/v2/articles/search',json={'resource_doi':'10.1158/2326-6066.CIR-24-1146','page_size':100},timeout=45);r.raise_for_status();a=r.json();(P/'metadata/RCC_figshare_search.json').write_text(json.dumps(a,indent=2));print('FIGSHARE',[(x['id'],x['title']) for x in a],flush=True)
except Exception as e:records.append({'figshare_error':str(e)})
(P/'metadata/FETCH_MANIFEST.json').write_text(json.dumps(records,indent=2))
print('FETCH_COMPLETE',flush=True)
