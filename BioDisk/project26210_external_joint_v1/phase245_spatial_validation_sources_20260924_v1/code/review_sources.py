from pathlib import Path
import json,datetime,re,xml.etree.ElementTree as ET,zipfile,io,hashlib,subprocess
import pandas as pd
from scipy.io import whosmat
P=Path(__file__).resolve().parents[1];ROOT=P.parent
def now():return datetime.datetime.now().astimezone().isoformat()
def dump(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2))
r=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True);(P/'logs/review_memory.log').write_text(r.stdout+r.stderr);assert r.returncode==0
members=pd.read_csv(P/'tables/MSI_ARCHIVE_ALL_MEMBERS.tsv',sep='\t')
mass=pd.read_csv(P/'tables/MSI_FRAGMENTATION_EVIDENCE.tsv',sep='\t',dtype={'ms_levels':str})
assert len(members)==75 and len(mass)==28
cache=[]
with zipfile.ZipFile(ROOT/'phase241_spatial_MSI_download_20260924_v1/downloads/spatial_metabolomics.zip') as z:
 for name in z.namelist():
  if not name.endswith('.ibh'):continue
  b=z.read(name)
  try:fields=str(whosmat(io.BytesIO(b)))
  except Exception as e:fields='format inspection: '+str(e)
  cache.append(dict(member=name,bytes=len(b),header=b[:80].decode('ascii',errors='replace'),fields=fields,interpretation='MSiReader parsing cache; not an independent optical image',source='https://msireader.com/docs/MSiReader_UserManual/MSiReader_User_Manual-21.html'))
pd.DataFrame(cache).to_csv(P/'tables/IBH_CACHE_REVIEW.tsv',sep='\t',index=False)
bio=json.loads((P/'sources/S-BIAD1661.json').read_text())
strings=[]
def walk(x):
 if isinstance(x,dict):
  for k,v in x.items():walk(v)
 elif isinstance(x,list):
  for v in x:walk(v)
 elif isinstance(x,str):strings.append(x)
walk(bio)
relevant=[s for s in strings if re.search(r'CD68|CD3|CD8|SPP1|osteopontin|multiplex|CD11b|antibod',s,re.I)]
dump(P/'tables/BIOIMAGE_PANEL_EVIDENCE.json',relevant)
xml=ET.parse(P/'sources/RCC_SPP1_CD44_fulltext.xml')
paras=[' '.join(x.itertext()) for x in xml.iter() if x.tag in ['p','sec'] and any(k in ' '.join(x.itertext()).lower() for k in ['data availability','gse210041','available upon request','available from'])]
dump(P/'tables/RCC_SPP1_PAPER_ACCESS_EVIDENCE.json',paras[:20])
alltext=' '.join(xml.getroot().itertext());accessions=sorted(set(re.findall(r'GSE\d+',alltext)))
rows=[
 dict(source='GSE294006',material='3 public Visium FFPE treated tumors; 27 processed files',eligibility='INDEPENDENT_DESCRIPTIVE_CASES',reason='new source versus prior GLMF/GSE175540/TJRCC; processed data downloaded and Phase247 analyses; three cases and mixed spots cannot validate cell-specific contacts',url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294006'),
 dict(source='RCC_publication_IMC',material='498-patient study code and README',eligibility='RAW_DATA_NOT_PUBLIC_PER_AUTHOR',reason='author README says raw data released upon publication; no accessible raw segmented-cell table obtained',url='https://github.com/BodenmillerGroup/RCC_publication'),
 dict(source='S-BIAD1661',material='tumoroid mIF and DESI-MSI',eligibility='WRONG_MARKER_PANEL_FOR_SPP1_CD8',reason='available panel lists CD56,CD20,CD68,CD3,PCK,CD11b; no SPP1/CD8 listed in retrieved panel; not direct validation of SPP1-CD8 relationship',url='https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD1661'),
 dict(source='PMC11687132',material='published SPP1-CD44 RCC study',eligibility='LITERATURE_SUPPORT_NOT_NEW_RAW_VALIDATION',reason='retrieved accession list '+','.join(accessions)+'; this accession is already in the current spatial cohort; published selected images do not provide our patient-level raw spatial identity validation',url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11687132/fullTextXML'),
 dict(source='FAP_mIF_39751643',material='FAP/lymphoid/myeloid spatial study',eligibility='NO_DIRECT_SPP1_PANEL_VERIFIED',reason='available publication description supports FAP/CD8/CD68, not verified SPP1 co-measurement and reusable patient-level raw data; no false independent validation claim',url='https://pubmed.ncbi.nlm.nih.gov/39751643/')]
pd.DataFrame(rows).to_csv(P/'tables/EXTERNAL_VALIDATION_QUALIFICATION.tsv',sep='\t',index=False)
conclusions={'at':now(),'archive_members':75,'imzml':28,'ibd':28,'ibh_cache':len(cache),'other_file_types':sorted(set(Path(n).suffix for n in members.member)-{'.imzML','.ibd','.ibh'}),'ms_level_values':sorted(mass.ms_levels.unique().tolist()),'precursor_fragmentation_CV_count':int(mass.precursor_or_fragmentation_cv.sum()),'independent_optical_or_author_transform_found_in_archive':False,'molecular_identity':'mass-compatible candidate only; no structural confirmation','optical_registration':'still geometry-only, no independent anatomical confirmation','C4':'original mixture-recovery gate remains not passed','source_access_failure':'current figshare API403 retained; previously downloaded official archive and cached original metadata used; no bypass attempted'}
dump(P/'SOURCE_REVIEW_CONCLUSIONS.json',conclusions)
report=(P/'REPORT_CN.md').read_text()
report+='\n\n## 源核查结果与剩余材料\n完整MSI归档只有28个imzML、28个ibd和19个ibh缓存。28份imzML的ms_level字段均为0，未提供可用前体/碎片CV记录；不能将该字段解释为经确认MS1，更不能据此完成MS/MS结构确证。未找到独立MSI光学图、作者变换或标准品验证文件。当前正离子精胺/亚精胺继续使用质量相容候选名称。\n\n.ibh是MSiReader为加速imzML加载保存的解析缓存，本轮检查文件头与变量清单，不将其当成独立光学图。官方说明：https://msireader.com/docs/MSiReader_UserManual/MSiReader_User_Manual-21.html 。\n\n新增GSE294006三病例已进入Phase247；498患者IMC源目前原始数据未公开；S-BIAD1661有CD68/CD3等染色但检索到的面板缺少SPP1和CD8，不能直接回答本项目所需关系；旧SPP1-CD44论文所列GSE175540已在当前空间队列中，不能作为新增独立队列。\n\n尚需的独立材料：同片SPP1与巨噬、CD8谱系标记及细胞分割/坐标；血管/成纤维标记与病理审核边界；MSI采集切片光学图或作者配准变换；与候选峰对应的MS/MS或标准品证据。这些材料未取得，不能通过程序打分或放宽模型门槛补成已完成。客户既有湿实验不重新分析或索取原始数据。\n'
(P/'REPORT_CN.md').write_text(report)
dump(P/'FINAL_VALIDATION.json',dict(status='PASS_QUALIFICATION_REVIEW',at=now(),source_requests=9,archive_members=75,imzml_metadata=28,ibh_caches=len(cache),external_sources_reviewed=5,independent_staining_or_pathology_obtained=False,independent_registration_obtained=False,structural_peak_confirmation=False,new_figures=0))
dump(P/'STATUS.json',dict(status='COMPLETE_WITH_DOCUMENTED_EVIDENCE_GAPS',at=now()));(P/'RUN_FINISHED').write_text(now())
print(json.dumps(conclusions));print('ACCESSIONS',accessions);print('BIO_PANEL',[s[:600] for s in relevant])

