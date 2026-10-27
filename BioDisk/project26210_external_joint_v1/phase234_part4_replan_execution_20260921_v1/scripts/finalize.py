from common import *
import shutil,csv,zipfile
def main():
 val=json.loads((P/'FINAL_VALIDATION.json').read_text());assert val['status']=='PASS'
 assert json.loads((P/'SUPPLEMENT_VALIDATION.json').read_text())['status']=='PASS'
 for name in ['A0','finalize_freeze','broad','narrow','validation']:
  assert (P/(name+'_EXIT_STATUS')).read_text().strip()=='0',name
 # Clarify zero-eligibility separately from small-n descriptive evidence.
 t=P/'tables/ALL_ENDPOINT_EVIDENCE_REGISTER.tsv';d=pd.read_csv(t,sep='\t',float_precision='round_trip')
 d.loc[d.n==0,'evidence_level']='NOT_ESTIMABLE';d.to_csv(t,sep='\t',index=False)
 note='''# 实施与验收记录
- 原始spot直径元数据与标准Visium物理间距不一致；在查看新关联前通过array_row/array_col及图像仿射核验完成标定，保留初始直径诊断。
- 宽层目录增加绝对血管q05目标，避免将源绝对丰度—目标构成误写成双绝对量关联；加入时间早于A1计算。
- 配对效应变化的FDR最初按来源与层级合池，交付前修正为来源×原端点家族×层级，原始P和效应未改；初稿留档。
- 独立验收初次发现TSV默认浮点读取可改变Wilcoxon绝对秩并列；改用round_trip读取后与存储前统计完全一致。未改变原始统计。
- 另一验收差异源于校验器未应用至少3个不同取值门槛；补齐原固定资格规则后通过，未修改效应或放宽资格。
- 正式验收PASS：独立精确符号检验1001项，数值效应9448项、最大误差<1e-15，92项原始计数程序重建，38图及333旧RDS保持。
- 完整受体复合体列被SPP1/POSTN共用时，其两配体资格已核对一致（38片均无不同门槛结果）；未造成覆盖冲突。
- 无图；所有非估计项与历史结果保留。可选MISTy未新增；旧COMMOT结果继承，当前LR仅为冻结方法的表达兼容性。
'''
 (P/'audit/IMPLEMENTATION_AND_VALIDATION_NOTES_CN.md').write_text(note)
 # Old error file remains auditable but is explicitly superseded by final validation and exit status.
 if (P/'ERROR_VALIDATION.txt').exists():
  shutil.copy2(P/'ERROR_VALIDATION.txt',P/'audit/LAST_VALIDATION_FAILURE_SUPERSEDED.txt')
 (P/'ERROR_RESOLUTION.json').write_text(json.dumps(dict(all_recorded_validation_failures_resolved=True,current_validation='FINAL_VALIDATION.json',exit_status=0)))
 required=['SUMMARY_CN.md','PART4_RESULTS_CN.md','MANUSCRIPT_SECTION4_DRAFT_CN.md','PROGRAM_DEFINITIONS.json','IMPLEMENTATION_FREEZE_BEFORE_ASSOCIATION.json','A0_FINAL_FREEZE.json','FINAL_VALIDATION.json','SUPPLEMENT_VALIDATION.json']
 for name in required:assert (P/name).stat().st_size>0
 status=dict(status='COMPLETE_VERIFIED',time=stamp(),scope='Phase233 executable A0-A5; optional MISTy not run; cell identity/histology-dependent endpoints explicit not estimable; A6 drawings await content review',sections=38,registered_patient_units=37,spots=114640,section_model_records=21128,estimable_section_models=19645,patient_summary_records=1668,old_RDS_preserved=333,new_figures=0,download_job_untouched=True)
 save('STATUS.json',status);(P/'EXIT_STATUS').write_text('0\n');(P/'COMPLETE_VERIFIED').write_text(status['time']+'\n')
 # All statistical outputs, relevant inherited tables, scripts and validation in one portable review archive.
 dest=P/'delivery';dest.mkdir(exist_ok=True)
 for name in ['SUMMARY_CN.md','PART4_RESULTS_CN.md']:shutil.copy2(P/name,dest/name)
 files=[P/n for n in required]+[P/'STATUS.json',P/'EXIT_STATUS',P/'COMPLETE_VERIFIED',P/'SOFTWARE.json']
 for sub in ['tables','legacy_tables','scripts']:
  files.extend(f for f in (P/sub).rglob('*') if f.is_file() and '__pycache__' not in str(f))
 files.extend(P/'audit'/n for n in ['PHYSICAL_CALIBRATION.tsv','INDEPENDENT_PATIENT_SIGN_CHECKS.tsv','NUMERICAL_EFFECT_CHECKS.tsv','RAW_SCORE_CHECKS.tsv','GRAPH_INTEGRITY.tsv','IMPLEMENTATION_AND_VALIDATION_NOTES_CN.md','SOURCE_MANIFEST.json'])
 hashes=[dict(path=str(f.relative_to(P)),bytes=f.stat().st_size,sha256=sha(f)) for f in files]
 save('OUTPUT_SHA256.json',hashes);files.append(P/'OUTPUT_SHA256.json')
 archive=dest/'Phase234_complete_results_and_statistics.zip'
 with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for f in files:z.write(f,f.relative_to(P).as_posix())
 with zipfile.ZipFile(archive) as z:
  assert z.testzip() is None
  for h in hashes:assert hashlib.sha256(z.read(h['path'])).hexdigest()==h['sha256']
 delivery=[dict(name=f.name,bytes=f.stat().st_size,sha256=sha(f)) for f in dest.iterdir() if f.is_file() and f.name!='DELIVERY_SHA256.json']
 (dest/'DELIVERY_SHA256.json').write_text(json.dumps(delivery,ensure_ascii=False,indent=2))
 # Save exact previous root record versions, and prepend new handoff without rewriting historical results.
 record=f"""## {status['time']} Phase234：按第四部分重订方案完成可执行分析与文字统计
用户恢复sgs10的26210项目并授权按Phase233开展分析。已完成A0输入/定义/物理尺度冻结、A1全组织及4x4/5x5、A2固定800/400/200μm及第一圈、A3资格/调整、A4合格表达边界和11条LR兼容性、A5来源整合与正文/图件内容方案。38片114640spot，37登记患者单位，GLMF2病例、GSE17554023患者、TJ12患者分开；70宽层端点。旧RDS/C4/原正文不改，无图。
比例尺诊断发现原spot直径元数据隐含约85μm邻点间距，已在新关联前用标准100μm Visium阵列行列坐标与图像仿射核验，38片通过，图不跨组织分量。
主结果：MTA转录背景—SPP1绝对q05，全组织GSE23/23正rho0.2256/P2.384e-7/FDR7.153e-7，TJ12/12正rho0.1634/P0.0004883/FDR0.005859；400μm仍同向。巨噬内比例未同步支持。绝对血管共同分布存在，比例/条件模型减弱或反向；CAF/ECM未得一致SPP1邻近增强。V4功能背景与V3方向分开；表达边界400μm肿瘤带SPP1比例偏低/基质偏高，合格患者<10，仅探索性描述。CD44无统一支持，整合素部分名义线索，完整复合体不足按不可估计。
所有正负、P/Wilcoxon/FDR/患者n、分母/尺度变化与旧证据保留。正式病理、CAF/CD8身份及细胞类型LR不合格单列；MISTy为可选补充未新增，旧COMMOT继承不改称新兼容性。A6须内容审核后绘图。
独立核验1001患者精确符号检验、9448效应（最大误差<1e-15）、92原始程序、38图、333旧RDS PASS，退出0。权威{P.name}/PART4_RESULTS_CN.md；SUMMARY_CN.md和delivery统计zip可交付。下载3002854/监听3002855本轮未操作。
当前断点：本轮可执行统计和文字已完成，待内容审阅；不增加阈值/半径追显著，不自动绘图。所有本轮计算已结束。
"""
 for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
  original=R/name;backup=P/'root_record_backups'/name
  if not backup.exists():shutil.copy2(original,backup)
  s=original.read_text()
  if 'Phase234：按第四部分重订方案完成' not in s:original.write_text(record+'\n\n'+s)
 manifest=R/'RESULTS_MANIFEST.tsv'
 with manifest.open('a',newline='') as f:
  w=csv.writer(f,delimiter='\t')
  for key in ['PART4_RESULTS_CN.md','SUMMARY_CN.md','tables/ALL_ENDPOINT_EVIDENCE_REGISTER.tsv','FINAL_VALIDATION.json','delivery/Phase234_complete_results_and_statistics.zip']:
   w.writerow(['phase234_'+Path(key).stem,'phase234_part4_replan_execution','report' if key.endswith('.md') else 'artifact',str(P/key),status['time'],'PASS',str(P/'FINAL_VALIDATION.json'),'Phase233 finite executable plan; preserve old results; no figures'])
 print(json.dumps(dict(status=status,delivery=delivery),ensure_ascii=False),flush=True)
if __name__=='__main__':main()
