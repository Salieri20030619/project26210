from pathlib import Path
import json,re,shutil,hashlib,subprocess,os
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');old=r/'phase332_abstract_barrier_focus_20261006_v5';d=r/'phase335_Fig2_manuscript_integration_20261020_v6'
for s in ['report','figures','code','qc','delivery']:(d/s).mkdir(parents=True,exist_ok=True)
b=json.loads((old/'report/report_materialized.json').read_text(encoding='utf-8'));before=json.loads(json.dumps(b));mapping={'C':'A','D':'C','E':'D','F':'E','G':'F'}
for x in b:
 for key in ['text','caption']:
  if key in x:x[key]=re.sub(r'图\s*2([C-G])',lambda m:'图2'+mapping[m[1]],x[key])
b[112]['text']=b[112]['text'].replace('，主线图谱摘要见图2A','')
b[113]['text']=b[113]['text'].replace('染色体层面的CNV代表性示例见图2B。','上皮细胞恶性身份及分析集定位见图2A，染色体层面的CNV代表性示例见图2B。')
# Verify the edge in full source, independently of the displayed subnetwork.
import csv
source=r/'phase127_section2_4_pyscenic_redraw_remote_v1/data/CANDIDATE_REGULON_TARGETS.tsv'
rows=list(csv.DictReader(source.open(encoding='utf-8'),delimiter='\t'));edge=[x for x in rows if x['method']=='GRNBoost2' and x['tf']=='ILF2' and x['target']=='TYMP'];assert len(edge)==1
b[133]['text']='调控子展示子网络包含71条唯一TF→靶基因边，其中18条由两种算法共同支持（图2E）。完整调控子结果另包含GRNBoost2识别的ILF2→TYMP预测连接，该边未纳入图2E的展示子网络。'
b[141]['text']='A：左侧为38,670个肿瘤组织上皮细胞的原CNV身份：恶性共同支持10,844个、仅inferCNV支持恶性1,106个、非恶性共同支持498个、方法不一致22,216个、未定4,006个。前两类构成11,950个恶性细胞分析集。右侧三图在该恶性集的同一原UMAP中，依次显示七个cNMF主导程序、F1使用量和归一化TYMP RNA；颜色分别表示程序身份、使用量和log1p表达。'
for j,letter in [(144,'C'),(145,'D'),(146,'E'),(147,'F')]:b[j]['text']=letter+b[j]['text'][1:]
b[146]['text']=b[146]['text'].replace('灰点为靶基因','灰蓝点为靶基因').replace('网络为预测调控连接','水平标签及左右分列仅为排版。网络为预测调控连接')
del b[143]
oldtables=[x for x in before if x['type']=='table'];newtables=[x for x in b if x['type']=='table']
assert len(newtables)==len(oldtables)==41
assert all(x['headers']==y['headers'] and x['rows']==y['rows'] for x,y in zip(oldtables,newtables))
for f in (old/'figures').iterdir():shutil.copy2(f,d/'figures'/f.name)
shutil.copy2(r/'phase334_Fig2_reassembled_20261020_v1/delivery/Fig2_reassembled_A4_v1.pdf',d/'figures/Figure2.pdf')
(d/'report/report_materialized.json').write_text(json.dumps(b,ensure_ascii=False,indent=2),encoding='utf-8')
s=(old/'code/render_report.py').read_text(encoding='utf-8').replace(old.name,d.name).replace('26210_manuscript_immune_figure_integrated_v5','26210_manuscript_Fig2_updated_v6')
s=s.replace("2:set('ABCDEFG')","2:set('ABCDEF')")
s=s.replace("  elif n==1:","  elif n==2:\n   anchors=[('A',10,19),('B',10,106),('C',76,103),('D',142,103),('E',10,202),('F',120,202)]\n  elif n==1:")
(d/'code/render_report.py').write_text(s,encoding='utf-8')
qc={'tables':41,'all_table_headers_rows_unchanged':True,'figure2_source_sha256':hashlib.sha256((d/'figures/Figure2.pdf').read_bytes()).hexdigest(),'other_figures_unchanged':all((d/'figures'/f.name).read_bytes()==f.read_bytes() for f in (old/'figures').iterdir() if f.name!='Figure2.pdf'),'panel_mapping':mapping,'cross_reference_correction':{'issue':'ILF2->TYMP is in full regulon source, not in displayed graph','source':str(source),'source_edge':edge[0],'resolved':'separate full-network claim from displayed-subnetwork reference'},'changes':[{'old_index':j,'before':a.get('text',a.get('caption')),'after':next((x.get('text',x.get('caption')) for x in b if x.get('_phase327_index')==a.get('_phase327_index')),None)} for j,a in enumerate(before) if a not in b]}
(d/'qc/CONTENT_VALIDATION.json').write_text(json.dumps(qc,ensure_ascii=False,indent=2),encoding='utf-8')
env=os.environ.copy();env['PYTHONPATH']=str(r/'phase189_integrated_TYMP_MTA_macrophage_report_20260917_v1/vendor')
subprocess.run(['python3',str(d/'code/render_report.py')],env=env,check=True)
shutil.copy2(d/'report/26210_manuscript_Fig2_updated_v6.pdf',d/'delivery/26210_manuscript_Fig2_updated_v6.pdf')
print('INTEGRATED: tables unchanged, Fig2 references and anchors updated')
