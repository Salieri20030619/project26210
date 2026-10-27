from pathlib import Path
import json,shutil,hashlib
P=Path(__file__).resolve().parents[1]
f=P/'report/report_materialized.json';b=json.loads(f.read_text());log=[]
shutil.copy2(f,P/'qc/materialized_before_final_fixes.json')
shutil.copy2(P/'report/26210_manuscript_immune_integrated_v2.pdf',P/'qc/first_reviewed_full.pdf')
(P/'qc/FIRST_RENDER_HASHES.json').write_text(json.dumps({x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in (P/'qc/render').glob('page-*.png')},indent=2))
for x in b:
 if 'text' in x:
  t=x['text'];new=t.replace('[24,25,19]','[19,24,25]').replace('[25,19]','[19,25]')
  if new.startswith('外部单细胞数据采用Obradovic'):new+=' 数据来源见文献[28,29]。'
  new=new.replace('临床数据沿用全文的CheckMate和JAVELIN Renal 101队列。','临床数据包括CheckMate和JAVELIN Renal 101队列，具体分析见1.15节。')
  new=new.replace('空间数据沿用全文的切片、细胞状态读数和组织坐标。','空间免疫分析使用1.16节所述切片、细胞状态读数和组织坐标。')
  if x.get('_base_index')==251:new=new.replace('直接T细胞绝对丰度没有相应下降','直接T细胞绝对丰度在上述SPP1邻域比较中没有相应下降')
  if new.startswith('正文区分本项目肾癌数据、公开肿瘤模型干预及文献引用。'):
   new+=' 其他癌种的SPP1与CD44功能源数据、宿主Spp1缺失单细胞资料分别保留为补充证据[26,27]。'
  if new!=t:log.append({'field':'text','before':t,'after':new});x['text']=new
 if x.get('_immune_table_id')==2:
  old=x['source'];x['source']='原224项免疫发现家族，11名患者；完整程序结果见配套统计资料。';log.append({'field':'source','before':old,'after':x['source']})
f.write_text(json.dumps(b,ensure_ascii=False,indent=2))
rf=P/'code/render_report.py';s=rf.read_text()
# Renderer-only line break keeps the accession intact without changing table data.
s=s.replace("def para(t,sty='p'):\n","def para(t,sty='p'):\n if t=='公开组织bulk PRJNA1223824':t='公开组织bulk\\nPRJNA1223824'\n")
rf.write_text(s)
(P/'qc/FINAL_EDITORIAL_FIXES.json').write_text(json.dumps(log,ensure_ascii=False,indent=2))
print('FINAL_FIXES',len(log))
