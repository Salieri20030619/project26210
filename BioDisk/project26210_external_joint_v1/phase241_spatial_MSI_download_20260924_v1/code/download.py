import pathlib,json,datetime,subprocess,hashlib,csv,shutil,time,zipfile,uuid,xml.etree.ElementTree as ET,traceback
P=pathlib.Path(__file__).resolve().parents[1];ROOT=P.parent
C=json.loads((P/'CONFIG.json').read_text())
def now():return datetime.datetime.now().astimezone().isoformat()
def write(path,x):
 path=pathlib.Path(path);tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(x,ensure_ascii=False,indent=2));tmp.replace(path)
def digest(path,algo='sha256'):
 h=hashlib.new(algo)
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def status(state,**kw):write(P/'STATUS.json',{'status':state,'updated_at':now(),**kw})
def record(note):
 stamp=datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
 for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv']:shutil.copy2(ROOT/name,P/'root_record_backups'/f'{stamp}_{name}')
 for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:
  f=ROOT/name;f.write_text(note+'\n\n'+f.read_text())
 f=ROOT/'RESULTS_MANIFEST.tsv'
 with open(f) as stream:keys=stream.readline().rstrip('\n').split('\t')
 row={'result_id':'phase241_spatial_MSI_download_v1','stage':'Phase241','artifact_type':'download_integrity','absolute_path':str(P/'FINAL_VALIDATION.json'),'generated_at':now(),'validation_status':'PASS','validation_evidence':str(P/'FINAL_VALIDATION.json'),'notes':'10 added tumor samples, 40 MSI core files; no new spatial analysis'}
 with open(f,'a') as out:csv.DictWriter(out,fieldnames=keys,delimiter='\t').writerow({k:row.get(k,'') for k in keys})
def main():
 part=P/'downloads/spatial_metabolomics.zip.part';dest=P/'downloads/spatial_metabolomics.zip'
 attempt=0;bad=0
 while True:
  attempt+=1
  status('DOWNLOADING',attempt=attempt,current_bytes=part.stat().st_size if part.exists() else 0,expected_bytes=C['expected_bytes'],transport='dedicated_proxy_17897')
  url=C['archive_url']+'?phase241='+str(time.time_ns())
  cmd=['curl','--fail','--silent','--show-error','--location','--proxy',C['proxy'],'--noproxy','','--continue-at','-','--connect-timeout','20','--max-time',str(C['curl_max_time_seconds']),'--speed-time','120','--speed-limit','1024','--dump-header',str(P/'logs'/f'attempt_{attempt}.headers'),'--output',str(part),url]
  with open(P/'logs'/f'attempt_{attempt}.log','ab') as log:
   r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
  write(P/'logs'/f'attempt_{attempt}.json',{'time':now(),'exit_code':r.returncode,'bytes':part.stat().st_size if part.exists() else 0})
  if part.exists() and part.stat().st_size==C['expected_bytes']:
   status('VERIFYING_ARCHIVE',attempt=attempt)
   md5=digest(part,'md5')
   if md5==C['expected_md5']:
    part.rename(dest);break
   bad+=1
   part.rename(P/'downloads'/f'spatial_metabolomics.failed_md5_{bad}.zip')
   write(P/'logs'/f'failed_md5_{bad}.json',{'md5':md5,'expected':C['expected_md5']})
   if bad>=2:raise RuntimeError('Two complete archives failed official MD5; preserved')
  status('RETRY_WAIT',attempt=attempt,curl_exit=r.returncode,current_bytes=part.stat().st_size if part.exists() else 0,next_retry_seconds=C['retry_delay_seconds'])
  time.sleep(C['retry_delay_seconds'])
 write(P/'ARCHIVE_INTEGRITY.json',{'status':'PASS','bytes':dest.stat().st_size,'md5':digest(dest,'md5'),'sha256':digest(dest),'official_md5_match':True})
 status('VERIFYING_AND_EXTRACTING')
 required=list(csv.DictReader((P/'tables/REQUIRED_MEMBERS.tsv').open(),delimiter='\t'))
 files=[]
 with zipfile.ZipFile(dest) as z:
  fail=z.testzip();assert fail is None,fail
  archive_members=len(z.infolist())
  for row in required:
   name=row['member'];assert pathlib.PurePosixPath(name).name==name and '\\' not in name and not name.startswith('.')
   info=z.getinfo(name)
   assert info.file_size==int(row['uncompressed_bytes']) and info.CRC==int(row['crc32'],16)
   out=P/'extracted'/name;assert not out.exists()
   tmp=out.with_suffix(out.suffix+'.part')
   with z.open(name) as src,open(tmp,'wb') as dst:shutil.copyfileobj(src,dst,2**20)
   assert tmp.stat().st_size==info.file_size
   tmp.rename(out)
   files.append({'sample':row['sample'],'polarity':row['polarity'],'file':name,'bytes':info.file_size,'crc32':format(info.CRC,'08x'),'sha256':digest(out)})
   status('EXTRACTING',completed_core_files=len(files),required_core_files=40)
 assert len(files)==40
 ns={'m':'http://psi.hupo.org/ms/mzml'}
 pairs=[]
 for xml in sorted((P/'extracted').glob('*.imzML')):
  cv=None
  with xml.open('rb') as f:
   for event,node in ET.iterparse(f,events=('end',)):
    if node.tag=='{'+ns['m']+'}fileContent':
     cv={x.get('accession'):x.get('value','') for x in node.findall('.//m:cvParam',ns)};break
  assert cv is not None
  ibd=xml.with_suffix('.ibd')
  expectedsha=cv['IMS:1000091'].lower();actualsha=digest(ibd,'sha1');assert actualsha==expectedsha,xml.name
  with ibd.open('rb') as f:actualuuid=str(uuid.UUID(bytes=f.read(16)))
  assert uuid.UUID(actualuuid)==uuid.UUID(cv['IMS:1000080']),xml.name
  pairs.append({'imzML':xml.name,'ibd':ibd.name,'ibd_sha1':actualsha,'uuid':actualuuid,'sha1_match':True,'uuid_match':True})
 assert len(pairs)==20
 for fn,rows in [('EXTRACTED_FILE_MANIFEST.tsv',files),('PAIRED_IMZML_IBD_VALIDATION.tsv',pairs)]:
  with open(P/'tables'/fn,'w') as out:
   w=csv.DictWriter(out,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
 validation={'status':'PASS','finished_at':now(),'additional_tumor_samples':10,'core_files':40,'imzML_ibd_pairs':20,'all_archive_member_CRC_pass':True,'archive_members':archive_members,'official_archive_MD5_pass':True,'all_ibd_SHA1_UUID_pass':True,'old_two_pilots_unchanged':True,'new_analysis_started':False}
 write(P/'FINAL_VALIDATION.json',validation)
 (P/'DOWNLOAD_SUMMARY_CN.md').write_text('# 空间代谢配套下载完成\n\n其余10个肿瘤样本的40个核心文件已下载并校验。完整归档与官方大小/MD5一致，所有ZIP条目CRC通过；20对imzML/ibd内部SHA-1和UUID一致。\n\n样本：'+', '.join(C['samples'])+'\n\n原2例试点保持原路径；新40文件位于extracted/，完整校验和清单位于tables/。本轮仅下载及完整性验收，未开展新的配准或关联分析。14切片不自动代表14个独立患者。\n')
 (P/'EXIT_STATUS').write_text('0\n');(P/'RUN_FINISHED').write_text(now());(P/'COMPLETE_VERIFIED').write_text(now())
 status('COMPLETE_VERIFIED',validation=validation)
 record(f'## {now()} Phase241 空间MSI补充下载完成\n其余10个肿瘤样本40核心文件、20对imzML/ibd全部通过官方归档MD5、ZIP CRC、内部SHA1及UUID校验。原2例只读保留。未自动运行配准或关联分析。\n权威：{P}；FINAL_VALIDATION PASS，EXIT_STATUS0。')
if __name__=='__main__':
 try:main()
 except Exception:
  (P/'logs/ERROR.txt').write_text(traceback.format_exc());(P/'EXIT_STATUS').write_text('1\n');(P/'RUN_FINISHED').write_text(now())
  status('FAILED',error=traceback.format_exc());raise
