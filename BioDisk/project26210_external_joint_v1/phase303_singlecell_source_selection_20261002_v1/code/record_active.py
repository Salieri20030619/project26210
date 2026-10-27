from pathlib import Path
import json,datetime,os,shutil
P=Path(__file__).resolve().parents[1];R=P.parent
now=datetime.datetime.now().astimezone().isoformat()
entry=f'''## Phase303（{now[:10]}）：用户授权先找合适单细胞再依次推进主线
用户明确“先找到合适的单细胞，然后依次按这个主线推进”。执行顺序：新数据资格与谱系身份→TYMP—SPP1→巨噬—基质过程→免疫功能。允许本轮必要公开数据获取及版本化分析；旧C4与原对象保持，先文字统计不出图。
当前目录：{P.name}。官方元数据已核查4候选，原始计数已完整下载GSE210038/GSE224630/GSE222703（三归档共1.20GB）；GSE222703作者对象也已取得。GSE202813骨转移/骨髓语境与已用GSE178481混合公开资源需去重，当前不作为原发肾脏主输入。
GSE210038有7肾肿瘤患者，其中T1既往Nivo+Ipi、T3既往Sunitinib，主集合仅5例Naive；N1是三样本池不作独立患者。GSE224630为5未治疗患者，RCC4f冻存技术样本及RCC5血栓不加n。GSE222703为3患者、分选富集，不将比例当原组织比例；治疗在GEO样本未声明，保持Unknown。GSE210038与GLMF空间病例有来源联系，不能称独立多组学验证。
输入资格规则位于resources/QUALIFICATION_PROTOCOL.json，样本选择在tables/SAMPLE_PLAN.tsv。QC和身份分群进行中，TYMP/SPP1不用于HVG身份分群，尚未检验关联。先确认真实成纤维、壁细胞、巨噬及T/CD8，按身份和样本资格冻结，不按显著性筛样。
运行：code/qc_sources.py 与 code/export_author.R；日志logs/qc.log、logs/author.log，完成必须检查各exit及qc标记。作者RDS有两层gzip，首轮读取失败已诊断并独立解压，原始下载及日志保留。专用SSH代理隧道由本轮启动，下载通过服务器127.0.0.1:7897。
本轮执行尚未完成，后续必须接续，不仅交付计划。

'''
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 f=R/name;before=f.read_bytes();bak=P/'root_record_backups'/('START_'+name)
 with bak.open('xb') as o:o.write(before)
 tmp=R/(name+'.phase303.start.tmp')
 with tmp.open('xb') as o:o.write(entry.encode()+before)
 assert f.read_bytes()==before;shutil.copymode(f,tmp);os.replace(tmp,f)
(P/'STATUS.json').write_text(json.dumps(dict(status='RUNNING_QUALIFICATION',at=now,steps=['data qualification','TYMP-SPP1','macrophage-stroma','immune function'],current='data qualification'),indent=2))
print('ACTIVE_STATE_RECORDED')
