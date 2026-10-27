from pathlib import Path
import requests,json,hashlib,concurrent.futures
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase307_spatial_perturbation_and_RCC_protein_20261002_v1')
for d in ['code','downloads','metadata','resources','tables','qc','logs','delivery']:(P/d).mkdir(parents=True,exist_ok=True)
def session():
 s=requests.Session();s.trust_env=False;s.proxies={'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'};return s
urls={
'RCC_SourceData.xlsx':'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-026-76858-w/MediaObjects/41467_2026_76858_MOESM10_ESM.xlsx',
'RCC_article.html':'https://www.nature.com/articles/s41467-026-76858-w',
'SPAC_portal.html':'https://spac.pku-genomics.org/',
'SPAC_code_tree.json':'https://api.github.com/repos/zenglab-pku/SPAC-seq/git/trees/main?recursive=1',
'SPAC_article.html':'https://www.cell.com/cell/fulltext/S0092-8674(26)00516-7',
'RCC_NCBI_SRR32360655.xml':'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=SRR32360655&rettype=full&retmode=xml',
'RCC_SRA_run.csv':'https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=SRR32360655',
}
def fetch(item):
 k,u=item
 try:
  r=session().get(u,timeout=50);f=P/('downloads' if k.endswith('.xlsx') else 'metadata')/k;f.write_bytes(r.content);o={'file':str(f.relative_to(P)),'url':u,'HTTP':r.status_code,'bytes':len(r.content),'sha256':hashlib.sha256(r.content).hexdigest()};print(json.dumps(o),flush=True);return o
 except Exception as e:return {'file':k,'url':u,'error':str(e)}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:records=list(ex.map(fetch,urls.items()))
(P/'metadata/FETCH_MANIFEST.json').write_text(json.dumps(records,indent=2));print('FETCH_DONE',flush=True)
