from pathlib import Path
import os,sys,json,datetime,subprocess,hashlib,traceback,shutil,csv
P=Path(__file__).resolve().parents[1];ROOT=P.parent
def now():return datetime.datetime.now().astimezone().isoformat()
def write(f,x):
 f=Path(f);tmp=f.with_suffix(f.suffix+'.tmp');tmp.write_text(json.dumps(x,ensure_ascii=False,indent=2));tmp.replace(f)
def status(s,**kw):write(P/'STATUS.json',{'status':s,'updated_at':now(),**kw})
def main():
 steps=['preflight','parse_spectra','prepare_masks','register','export_candidates','analyze_regions','patient_summary']
 for i,name in enumerate(steps,1):
  status('RUNNING',stage=name,stage_index=i,total_stages=len(steps))
  with open(P/'logs'/f'{name}.log','wb') as log:
   r=subprocess.run([sys.executable,'-u',str(P/'code'/f'{name}.py')],stdout=log,stderr=subprocess.STDOUT)
  (P/'logs'/f'{name}.exit_code').write_text(str(r.returncode)+'\n')
  if r.returncode:raise RuntimeError(name+' failed; see stage log')
 status('FINAL_VALIDATION')
 qa_names=['PREFLIGHT_VALIDATION.json','RAW_SPECTRA_VALIDATION.json','REGISTRATION_VALIDATION.json','ANALYSIS_QA.json','PATIENT_SUMMARY_QA.json']
 q={name:json.loads((P/name).read_text()) for name in qa_names}
 assert all(v['status']=='PASS' for v in q.values())
 assert q['RAW_SPECTRA_VALIDATION.json']['spectra_files']==24
 for x in json.loads((P/'audit/RAW_INPUT_MANIFEST.json').read_text()):
  f=Path(x['path']);assert f.stat().st_size==x['bytes'] and f.stat().st_mtime_ns==x['mtime_ns'],str(f)
 for x in json.loads((P/'audit/OLD_RDS_BEFORE.json').read_text()):
  f=Path(x['path']);assert f.stat().st_size==x['bytes'] and f.stat().st_mtime_ns==x['mtime_ns'],str(f)
 assert not list(P.rglob('*.rds'))
 outfiles=[f for folder in ['tables','transforms','audit'] for f in (P/folder).glob('*') if f.is_file()]
 def sha(f):
  h=hashlib.sha256()
  with f.open('rb') as s:
   for b in iter(lambda:s.read(2**20),b''):h.update(b)
  return h.hexdigest()
 write(P/'OUTPUT_SHA256.json',{str(f.relative_to(P)):sha(f) for f in outfiles})
 final={'status':'PASS','finished_at':now(),'patients':12,'spectra':24,'stages_passed':7,'registration_candidates':q['REGISTRATION_VALIDATION.json']['retained_candidates'],'regional_effect_records':q['ANALYSIS_QA.json']['all_effect_records'],'patient_endpoints_planned':12,'old_RDS_unchanged':True,'raw_inputs_unchanged':True,'new_figures':0,'anatomically_validated':False,'output_manifest_files':len(outfiles)}
 write(P/'FINAL_VALIDATION.json',final)
 stamp=datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
 for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv']:shutil.copy2(ROOT/name,P/'root_record_backups'/f'{stamp}_{name}')
 note=f'## {now()} Phase242 十二例空间代谢分析完成\n12登记肿瘤患者、24正负谱，谱图QC、形态候选配准、区域关联及患者层汇总完成。配准未获独立解剖验证，质量峰为候选；保留所有阴性、朝向变号及不可估计结果。旧RDS/原两例和原研究结果保留，无新图。\n权威报告：{P}/REPORT_CN.md；完整表、候选变换和原两例重复比较均已保存；FINAL_VALIDATION PASS。'
 for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:
  f=ROOT/name;f.write_text(note+'\n\n'+f.read_text())
 f=ROOT/'RESULTS_MANIFEST.tsv'
 with f.open() as stream:keys=stream.readline().rstrip('\n').split('\t')
 row={'result_id':'phase242_12tumor_spatial_analysis_v1','stage':'Phase242','artifact_type':'analysis_report','absolute_path':str(P/'REPORT_CN.md'),'generated_at':now(),'validation_status':'PASS','validation_evidence':str(P/'FINAL_VALIDATION.json'),'notes':'12 tumor patients; candidate-only geometric registration; negative results retained; no new figures'}
 with f.open('a') as out:csv.DictWriter(out,fieldnames=keys,delimiter='\t').writerow({k:row.get(k,'') for k in keys})
 (P/'EXIT_STATUS').write_text('0\n');(P/'RUN_FINISHED').write_text(now());(P/'COMPLETE_VERIFIED').write_text(now())
 status('COMPLETE_VERIFIED',validation=final)
if __name__=='__main__':
 try:main()
 except Exception:
  (P/'logs/ERROR.txt').write_text(traceback.format_exc());(P/'EXIT_STATUS').write_text('1\n');(P/'RUN_FINISHED').write_text(now())
  status('FAILED',error=traceback.format_exc());raise
