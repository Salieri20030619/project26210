from pathlib import Path
import json,datetime,os,shutil,csv,io
P=Path(__file__).resolve().parents[1];R=P.parent
assert json.loads((P/'qc/FINAL_VALIDATION.json').read_text())['status']=='PASS'
for name in ['sequential','validation','delivery']:assert (P/f'logs/{name}.exit').read_text().strip()=='0'
now=datetime.datetime.now().astimezone().isoformat()
entry=f'''## Phase303终态（{now[:10]}）：新单细胞筛选与三段主线首轮分析完成
用户“先找到合适的单细胞，然后依次按这个主线推进”。已获取GSE210038/GSE224630/GSE222703原始计数并完成QC/身份核查，共15肿瘤文库、81535细胞。权威目录：{P.name}。
首选GSE210038五例未治疗患者作方向复核：四类细胞可同患者配对，但严格成纤维候选>=20仅1例；成纤维/壁细胞共标记池单列，不据此断言纯CAF或谱系转分化。T1/T3既往治疗另列；其GLMF空间病例不能当独立跨组学验证。
GSE224630五例仅35成熟巨噬且CD8不足，限肿瘤—基质补充；GSE222703三例成熟巨噬61/116/44、严格成纤维候选19/3/1，限肿瘤—CD8补充。GSE202813骨转移/骨髓情境未并主集；GSE159115已用；GSE171306供者独立性未清；HRA007082为更适配的19例ccRCC/22例RCC受控候选，未取得原始对象、未发送访问申请。
依次完成1 TYMP—SPP1、2基质程序、3免疫功能，冻结方案后计算。GSE210038患者内TYMP—SPP1：n4/4正，中位rho0.069585；IFN/缺氧调整4/4正、0.061150。T4仅9个SPP1检出细胞，不合患者内>=10门槛但合患者汇总>=5门槛。既往Phase280 C4/WholeMacro阳性不改、不与新增4例混合重检验。
GSE210038跨患者n5：肿瘤样TYMP—巨噬SPP1 rho-0.2；SPP1—共标记基质胶原成熟/黏着斑/收缩-0.3/+0.6/+0.5；胶原成熟—CD8细胞毒/功能障碍+0.2/+0.1。SPP1—CD8两项+0.3/+0.6。没有形成一致屏障链，不能只取正方向。严格成纤维链配对不足属于不可评估。各新研究有效n<10，全部P/FDR空，只报告效应/方向，未用细胞P补患者n。
交付delivery/26210_Phase303_单细胞筛选与主线首轮结果.txt及html、13张完整表/冻结规则/代码/验收；120个研究×端点组合含全部不可评估项。独立原文库汇总/全部rho/患者资格重算1273项通过，qc/FINAL_VALIDATION.json；正式程序exit0，初次作者gzip/样本映射及shell转义失败保留。原C4/旧对象/空间定义/正式文章及Phase302客户Word未改，无新增生物学图件。
本轮首轮筛选和可评估三段计算完成；尚未找到足够纯CAF配对的完整新验证队列，不把该缺口写为已完成验证。下一增量是新的合格配对输入，不重复当前5例换阈值。无计算后台，本轮自建SSH代理隧道已关闭；进度权威只在服务器，结果包交付本地Downloads/26210/Phase303及workspace。

'''
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 f=R/name;before=f.read_bytes();bak=P/'root_record_backups'/('COMPLETE_'+name)
 with bak.open('xb') as o:o.write(before)
 tmp=R/(name+'.phase303.complete.tmp')
 with tmp.open('xb') as o:o.write(entry.encode()+before)
 assert f.read_bytes()==before;shutil.copymode(f,tmp);os.replace(tmp,f)
f=R/'RESULTS_MANIFEST.tsv';before=f.read_bytes();bak=P/'root_record_backups/COMPLETE_RESULTS_MANIFEST.tsv'
with bak.open('xb') as o:o.write(before)
headers=before.decode('utf-8-sig').splitlines()[0].split('\t');s=io.StringIO();wr=csv.DictWriter(s,fieldnames=headers,delimiter='\t',lineterminator='\n')
for k,rel,typ in [('report','delivery/26210_Phase303_单细胞筛选与主线首轮结果.txt','report'),('edges','tables/ALL_SEQUENTIAL_EDGES.tsv','table'),('within','tables/STEP1_WITHIN_MACRO_SUMMARY.tsv','table'),('lineage','tables/LINEAGE_COUNTS_BY_PATIENT.tsv','table'),('validation','qc/FINAL_VALIDATION.json','validation')]:
 wr.writerow(dict(result_id='phase303_'+k,stage='phase303',artifact_type=typ,absolute_path=str(P/rel),generated_at=now,validation_status='PASS',validation_evidence=str(P/'qc/FINAL_VALIDATION.json'),notes='New single-cell qualification and sequential descriptive analyses; pure CAF chain remains underpowered; no patient inference'))
tmp=R/'RESULTS_MANIFEST.tsv.phase303.complete.tmp';tmp.write_bytes(before+(b'' if before.endswith(b'\n') else b'\n')+s.getvalue().encode());assert f.read_bytes()==before;shutil.copymode(f,tmp);os.replace(tmp,f)
(P/'STATUS.json').write_text(json.dumps(dict(status='FIRST_ROUND_COMPLETE_WITH_FULL_CHAIN_INPUT_GAP',at=now,completed_steps=['source qualification','TYMP-SPP1','macrophage-stroma','immune function','independent validation','delivery'],unfinished_scientific_question='No sufficiently replicated strict fibroblast complete chain; no causal/physical barrier proof',background_jobs=[]),ensure_ascii=False,indent=2))
print('SERVER_RECORDS_UPDATED')
