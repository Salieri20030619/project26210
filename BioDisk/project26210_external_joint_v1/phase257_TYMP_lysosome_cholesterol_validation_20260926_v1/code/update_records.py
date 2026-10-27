from pathlib import Path
import json,hashlib,datetime,shutil,csv,fcntl,os
p=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase257_TYMP_lysosome_cholesterol_validation_20260926_v1');root=p.parent
qc=json.loads((p/'qc/FINAL_VALIDATION.json').read_text());assert qc['status']=='PASS'
assert (p/'analysis_EXIT_STATUS').read_text().strip()=='0' and (p/'validation_EXIT_STATUS').read_text().strip()=='0'
now=datetime.datetime.now().astimezone().isoformat()
entry=f"""## {now} Phase257 TYMP—溶酶体/胆固醇—SPP1定向验证及独立小节完成
用户提出可加独立小节，先验证机制联系。本轮38片37患者空间检验+16位同文库细胞区室伪bulk；固定Phase256的Lysosome/Cholesterol metabolism，去TYMP/SPP1及五状态基因，并去LIPA/NPC1/NPC2/SORT1共享基因作敏感性。无新图、无3G、无重训。
关键新结果：肿瘤TYMP与肿瘤内胆固醇程序研究调整rho0.729125，P0.00905/FDR0.01810；去共享基因rho0.704266，P0.0152/敏感池FDR0.03445。肿瘤内溶酶体rho0.552044/P0.06345未获支持。巨噬SPP1—胆固醇rho0.647233/P0.02815/FDR0.05630，只称名义显著；巨噬SPP1—溶酶体P0.0723未获支持。
空间TYMP—溶酶体GSE175540/TJRCC中位rho0.02975/0.03026、18/23及11/12正；TYMP—胆固醇rho0.01804/0.02014、17/23及12/12正，主家族FDR均支持但效应弱，去共享后方向保留。4x4调整未获患者层支持。两功能与旧SPP1-state丰度/巨噬内比例联系更一致，不等同新C4或单基因RNA。
同文库肿瘤TYMP与巨噬溶酶体/胆固醇/SPP1均未获支持（rho-0.02705/0.24497/-0.10810）。16患者中10来自GSE207493、4来自GSE156632、另两研究各1；仅14研究内可交换，删除主要研究后6例模型自由度不足，不称独立多研究验证。不把分段相关写成肿瘤TYMP驱动巨噬功能的因果中介。
独立草稿SECTION_DRAFT_CN.md含方法1.16、结果2.20、讨论3.11，建议标题为“TYMP与胆固醇处理程序相关，SPP1相关区域呈现相应功能特征”。联系TYMP缺陷/酶活恢复、核苷-溶酶体和OPN脂质调节文献，明确不同疾病模型；原MTA主线及既有报告不变，尚未把草稿写入原文。
权威目录：{p}。REPORT_CN.md、独立小节、主统计及完整统计复现ZIP已生成于delivery。1976空间相关(20项资格不足NA)、156队列汇总、35配对关联；76程序复算、1904患者合并、35矩阵偏相关及P/BH/置换核验通过，121输入记录前后不变。首次float32精度检查失败已以float64计算修复，不放宽容差，日志保留。
FINAL_VALIDATION PASS，分析和验收退出0，Python batch_bench/BLAS1，无新依赖。四权威记录只存服务器，本地仅交付。无运行中的本轮计算；后续可审阅独立小节，图件/3G尚未启动。
"""
lock=(root/'phase257_record_update.lock').open('a+');fcntl.flock(lock,fcntl.LOCK_EX)
backup=p/'qc/record_backup';backup.mkdir(exist_ok=True)
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:
 f=root/name;old=f.read_text()
 if 'Phase257 TYMP—溶酶体/胆固醇—SPP1定向验证及独立小节完成' in old:raise RuntimeError('already recorded')
 shutil.copyfile(f,backup/name);temp=f.with_name(name+'.phase257.tmp');temp.write_text(entry+'\n'+old,encoding='utf-8');os.replace(temp,f)
f=root/'RESULTS_MANIFEST.tsv';shutil.copyfile(f,backup/f.name)
with f.open('a',newline='') as h:
 w=csv.writer(h,delimiter='\t',lineterminator='\n')
 for i,m in enumerate(json.loads((p/'DELIVERY_MANIFEST.json').read_text()),1):
  w.writerow([f'phase257_TYMP_function_{i}',p.name,'targeted_validation_and_section',m['path'],now,'PASS_EXPLORATORY_ASSOCIATION',str(p/'qc/FINAL_VALIDATION.json'),f"Separate section draft; cross-compartment negative retained; no figures; SHA256={m['sha256']}"])
fcntl.flock(lock,fcntl.LOCK_UN)
(p/'qc/RECORD_UPDATE_VALIDATION.json').write_text(json.dumps({'status':'PASS','time':now,'updated':['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv'],'entries':4},indent=2))
print('Phase257 four authoritative records updated.')
