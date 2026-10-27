from pathlib import Path
import sys,json,datetime,hashlib,subprocess,traceback,zipfile,urllib.request,xml.etree.ElementTree as ET,re,concurrent.futures
import pandas as pd
P=Path(__file__).resolve().parents[1];ROOT=P.parent
def now():return datetime.datetime.now().astimezone().isoformat()
def dump(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2))
def fetch(item):
 name,url=item
 try:
  req=urllib.request.Request(url,headers={'User-Agent':'project26210-source-audit/1.0'})
  with urllib.request.urlopen(req,timeout=20) as r:body=r.read(12*1024*1024+1);status=r.status
  if len(body)>12*1024*1024:raise ValueError('metadata size limit')
  (P/'sources'/name).write_bytes(body)
  return dict(name=name,url=url,status='FETCHED',http=status,bytes=len(body),sha256=hashlib.sha256(body).hexdigest(),at=now())
 except Exception as e:return dict(name=name,url=url,status='UNAVAILABLE',error=str(e),at=now())
def main():
 r=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
 (P/'logs/memory_gate.log').write_text(r.stdout+r.stderr);assert r.returncode==0
 dump(P/'STATUS.json',dict(status='SOURCE_AUDIT',at=now()))
 urls=[
 ('figshare24599295.json','https://api.figshare.com/v2/articles/24599295'),
 ('GSE294006.soft.txt','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294006&targ=self&form=text&view=full'),
 ('GSE294006_files.html','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE294nnn/GSE294006/suppl/'),
 ('S-BIAD1661.json','https://www.ebi.ac.uk/biostudies/api/v1/studies/S-BIAD1661'),
 ('RCC_IMC_README.md','https://raw.githubusercontent.com/BodenmillerGroup/RCC_publication/main/README.md'),
 ('RCC_IMC_TREE.json','https://api.github.com/repos/BodenmillerGroup/RCC_publication/git/trees/main?recursive=1'),
 ('RCC_SPP1_CD44_fulltext.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11687132/fullTextXML'),
 ('RCC_FAP_mIF_pubmed.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:39751643&format=json'),
 ('SpaMTP_README.md','https://raw.githubusercontent.com/GenomicsMachineLearning/SpaMTP/main/README.md')]
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:res=list(pool.map(fetch,urls))
 dump(P/'sources/SOURCE_ACCESS_LOG.json',res)
 print(json.dumps(res,ensure_ascii=False),flush=True)
 archive=ROOT/'phase241_spatial_MSI_download_20260924_v1/downloads/spatial_metabolomics.zip'
 inventory=[];mass=[];snips=[]
 with zipfile.ZipFile(archive) as z:
  for f in z.infolist():
   inventory.append(dict(member=f.filename,bytes=f.file_size,compressed_bytes=f.compress_size,CRC=f.CRC))
   suffix=Path(f.filename).suffix.lower()
   if suffix=='.imzml':
    raw=z.read(f.filename);root=ET.fromstring(raw);cv=[c.attrib for c in root.iter() if c.tag.endswith('cvParam')]
    levels=sorted(set(c.get('value','') for c in cv if c.get('accession')=='MS:1000511'))
    tandem=[c for c in cv if any(k in c.get('name','').lower() for k in ['collision energy','selected ion m/z','activation','precursor'])]
    mass.append(dict(member=f.filename,ms_levels=';'.join(levels),precursor_or_fragmentation_cv=len(tandem),file_role='deposited processed MSI mass-intensity data'))
   elif suffix in ['.csv','.tsv','.txt'] and f.file_size<5*1024*1024:
    txt=z.read(f.filename).decode('utf-8',errors='replace')
    matches=[line[:1000] for line in txt.splitlines() if re.search(r'spermine|spermidine|MS/MS|MS2|standard|annotation|optical|registration',line,re.I)]
    snips.append(dict(member=f.filename,matches=matches[:30],first_lines=txt.splitlines()[:4]))
 pd.DataFrame(inventory).to_csv(P/'tables/MSI_ARCHIVE_ALL_MEMBERS.tsv',sep='\t',index=False)
 pd.DataFrame(mass).to_csv(P/'tables/MSI_FRAGMENTATION_EVIDENCE.tsv',sep='\t',index=False)
 dump(P/'tables/MSI_TEXT_MEMBER_EVIDENCE.json',snips)
 # Inspect author annotation and source-tree assets, preserving identifiers.
 treepath=ROOT/'phase235_SpaMTP_TJRCC_pilot_20260921_v1/metadata/TJRCC_AUTHOR_CODE_TREE.json'
 tree=json.loads(treepath.read_text())
 hits=[x for x in tree.get('tree',[]) if re.search(r'annot|align|regist|hist|optic|image|png|tif|supp|xlsx|csv|rds|rda',x.get('path',''),re.I)]
 dump(P/'tables/AUTHOR_SOURCE_RELEVANT_FILES.json',hits)
 source=ROOT/'phase235_SpaMTP_TJRCC_pilot_20260921_v1/metadata/figshare_article_24599295_cached.json'
 article=json.loads(source.read_text())
 dump(P/'tables/FIGSHARE_EXISTING_METADATA.json',{'title':article.get('title'),'description':article.get('description'),'references':article.get('references'),'files':article.get('files')})
 report=['# 空间独立验证资格与源材料核查','',f'检查时间：{now()}','',
 '## 可完成与尚缺的环节',
 '|环节|本轮行动|资格定位|','|---|---|---|',
 '|SPP1巨噬及接收端身份|Phase246以分离的谱系标记核对既有映射，保留CD8B/T谱系共同检出资格|RNA一致性，不是独立蛋白或细胞分割验证|',
 '|组织结构|固定400um邻域及表达富集区最短图距离；不同来源患者分别统计|表达区域结构，不改称CD8细胞排斥或CAF接触|',
 '|复合指标|原T-cell balance、去SPP1分母及直接T细胞丰度完整比较|用于检验共享成分依赖|',
 '|MSI配准|核查官方完整归档和作者材料；不按分子方向选择朝向|独立光学/变换依据需实际取得才能解锁|',
 '|质量峰身份|检查归档MS层级、碎片/前体信息和注释配套|准确质量候选不能升级结构确证|',
 '|独立空间队列|核查新公开GEO与空间蛋白来源的数据和所测标记资格|同一已用数据不重复算独立验证|',
 '|C4空间身份|保留Phase219/220留研究混合恢复失败证据|不放宽假阳性/研究门槛，不重训cell2location|','',
 '## C4当前数值门槛',
 'Phase220主深度2000 UMI完整面板零C4误检率16.67%（门槛10%）、最差研究巨噬内误差25百分点（门槛20），未通过。提高至10000 UMI后误检率8.97%，但研究误差仍25百分点。已有技术拟合不能作为C4空间丰度；原SPP1-state仍按原身份使用。','',
 '本轮源材料、访问失败及全部文件清单存于sources/和tables/；未联系作者、未申请受控数据，未将网页或论文展示图当作可复算独立原始数据。']
 (P/'REPORT_CN.md').write_text('\n'.join(report))
 dump(P/'FINAL_VALIDATION.json',dict(status='PASS_SOURCE_INVENTORY',at=now(),source_requests=len(urls),source_success=sum(x['status']=='FETCHED' for x in res),archive_members=len(inventory),imzml_inspected=len(mass),no_image_modification=True,no_identity_upgrade=True))
 dump(P/'STATUS.json',dict(status='SOURCE_AUDIT_COMPLETE_REVIEW_PENDING',at=now()));(P/'EXIT_STATUS').write_text('0\n')
 print('MEMBERS',len(inventory),'IMZML',len(mass))
try:main()
except BaseException:
 (P/'logs/ERROR.txt').write_text(traceback.format_exc());dump(P/'STATUS.json',dict(status='FAILED',error=traceback.format_exc()));(P/'EXIT_STATUS').write_text('1\n');raise

