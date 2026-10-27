from pathlib import Path
import json,hashlib,datetime,shutil,csv,io,fcntl,os
p=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase256_spatial_region_function_20260926_v1');root=p.parent
qc=json.loads((p/'qc/FINAL_VALIDATION.json').read_text());assert qc['status']=='PASS'
assert (p/'retry_EXIT_STATUS').read_text().strip()=='0'
now=datetime.datetime.now().astimezone().isoformat()
entry=f"""## {now} Phase256 SPP1相关区域功能分析完成（参考Fig.2G，3G未启动）
用户授权先做类似参考论文2G的分析，随后考虑3G。本轮38片/37登记患者：GLMF2、GSE17554023、TJRCC12分别报告；全部114个区域比较通过资格。主定义旧SPP1-state q05前20%对其余spot；组成/深度调整和巨噬内比例为敏感性。排除TYMP/SPP1/GPNMB/LGALS3/FABP5/CD9/LPL，保留完整基因版本；GO_BP/Hallmark/KEGG完整库，不只选免疫抑制条目。
最有力的功能图候选为溶酶体和胆固醇代谢。调整后溶酶体GSE17554022/23、TJRCC12/12正，患者中位Z效应+0.0201/+0.0253；胆固醇22/23、11/12正，+0.0303/+0.0325，患者层与GSEA均FDR支持。吞噬体/氧化磷酸化原区域阳性保留，但调整后平均功能差异减弱；不能只挑GSEA曲线。参考Fig2G十类免疫负调控条目在两来源主GSEA均未通过完整GO家族FDR；原始区域均值广泛升高与调整后反向结果同时保留。
区域伪bulk辅助支持两外部来源溶酶体/胆固醇方向；GLMF溶酶体2/2正、胆固醇1正1负、吞噬体2/2负，不能隐藏。巨噬内比例存在来源差异，不写成普遍状态诱导、免疫屏障或新C4空间定位；本轮无实测局部MTA。
GSEA P是给定患者排序中的基因集随机化；患者P为独立患者符号检验，两层分别分来源/功能库/主或敏感池校正，GLMF不推断。初轮780条数值NA已以nPermSimple100000定向复核，仍{qc['GSEA_not_estimable']}条保留NA，不作阴性；排序/集合/ES和假设家族未改。自建NPZ基因名格式失败和修复留档，数值未变。
权威目录：{p}。REPORT_CN.md、SUMMARY_CN.md及delivery/报告/关键表/完整统计复现ZIP已生成；{qc['functional_rows']}条功能记录、{qc['gsea_rows']}条GSEA记录，54组程序/ES独立核验，输入身份前后不变，FINAL_VALIDATION PASS，提取/分析/重试/验收退出0。Python batch_bench、R R_4.5.2，无新依赖。
新图0，3G尚未启动，旧RDS及既有Phase253报告/Phase255图件保持不变；当前待审核的是本轮文字统计及后续功能图选题。无运行的本轮计算任务。进度仅记录服务器。
"""
lock=(root/'phase256_record_update.lock').open('a+')
fcntl.flock(lock,fcntl.LOCK_EX)
backup=p/'qc/record_backup';backup.mkdir(exist_ok=True)
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:
 f=root/name;old=f.read_text()
 if 'Phase256 SPP1相关区域功能分析完成' in old:raise RuntimeError('Entry already present')
 shutil.copyfile(f,backup/name)
 tmp=f.with_name(name+'.phase256.tmp');tmp.write_text(entry+'\n'+old,encoding='utf-8');os.replace(tmp,f)
f=root/'RESULTS_MANIFEST.tsv';shutil.copyfile(f,backup/f.name)
with f.open(newline='') as h:header=next(csv.reader(h,delimiter='\t'))
assert header==['result_id','stage','artifact_type','absolute_path','generated_at','validation_status','validation_evidence','notes'],header
with f.open('a',newline='') as h:
 w=csv.writer(h,delimiter='\t',lineterminator='\n')
 for i,m in enumerate(json.loads((p/'DELIVERY_MANIFEST.json').read_text()),1):
  w.writerow([f'phase256_region_function_{i}','phase256_spatial_region_function_20260926_v1','functional_analysis_package',m['path'],now,'PASS_WITH_EXPLICIT_NUMERICAL_NA',str(p/'qc/FINAL_VALIDATION.json'),f"SPP1-region function; no new figures; 3G not started; SHA256={m['sha256']}"])
fcntl.flock(lock,fcntl.LOCK_UN)
(p/'qc/RECORD_UPDATE_VALIDATION.json').write_text(json.dumps({'status':'PASS','time':now,'updated':['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv'],'entries':4,'historical_records_backup':str(backup)},ensure_ascii=False,indent=2))
print('Four authoritative records updated; prior content preserved.')
