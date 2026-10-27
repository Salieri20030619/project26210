from pathlib import Path
import requests,json,hashlib,datetime,xml.etree.ElementTree as ET,tarfile,io
P=Path(__file__).resolve().parents[1]
S=requests.Session();S.trust_env=False;S.proxies={'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'}
plan='''用户授权：先找到合适单细胞，再依次按TYMP—SPP1、巨噬—基质过程、免疫功能主线推进。
先独立核查公开原始计数、患者/文库/部位/治疗信息、四类细胞身份与有效配对。候选优先ccRCC全组织、具成纤维与免疫测量；独立性先于关联，既有C4冻结，旧对象只读。
候选初始：GSE210038（7肿瘤、2邻近正常，myCAF背景）；GSE202813（骨转移/骨髓，须与既有GSE178481去重）；GSE222703（3患者、分选富集，组成比例不可直接比较）；GSE224630（5例未治疗肿瘤，资格待官方核查）。
依据资格决定数据，不根据TYMP/SPP1或关联P选择队列。取得后先QC/谱系身份和同患者资格，再冻结分析方案，之后按顺序计算。各研究、原发/转移、治疗时间独立记录；生物学样本为患者，n<10按项目约定方向和效应描述。推断以目标效应原P<0.05，FDR完整保留并按预设同类家族计算，不以FDR作为筛选必需门槛。
不要求客户新实验。先文字统计，无新增图；报告成果包含候选筛选、身份证据、患者资格、逐步全部结果和必要反向证据。
'''
(P/'PLAN_AND_SCOPE.txt').write_text(plan)
rows=[]
for acc in ['GSE210038','GSE202813','GSE222703','GSE224630']:
 base=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}'
 for suffix,out in [('/miniml/'+acc+'_family.xml.tgz',acc+'_family.xml.tgz'),('/suppl/filelist.txt',acc+'_filelist.txt')]:
  url=base+suffix
  try:
   res=S.get(url,timeout=(20,80));res.raise_for_status();data=res.content
   (P/'metadata'/out).write_bytes(data)
   rows.append(dict(accession=acc,url=url,path=str(P/'metadata'/out),bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),status='OK'))
   if out.endswith('.tgz'):
    with tarfile.open(fileobj=io.BytesIO(data)) as tar:
     f=next(x for x in tar.getmembers() if x.name.endswith('.xml'));raw=tar.extractfile(f).read()
    (P/'metadata'/f'{acc}.xml').write_bytes(raw)
    root=ET.fromstring(raw);ns={'m':'http://www.ncbi.nlm.nih.gov/geo/info/MINiML'}
    samples=[]
    for sample in root.findall('m:Sample',ns):
     samples.append(dict(accession=sample.attrib.get('iid'),title=sample.findtext('m:Title',default='',namespaces=ns),characteristics=[dict(tag=x.attrib.get('tag'),value=x.text) for x in sample.findall('.//m:Characteristics',ns)],source=[x.text for x in sample.findall('.//m:Source',ns)],processing=[x.text for x in sample.findall('.//m:Data-Processing',ns)],files=[x.text for x in sample.findall('.//m:Supplementary-Data',ns)]))
    (P/'metadata'/f'{acc}_samples.json').write_text(json.dumps(samples,ensure_ascii=False,indent=2))
    print(acc,'samples',len(samples),[(x['accession'],x['title']) for x in samples][:30],flush=True)
   else:print(acc,res.text[:1800],flush=True)
  except Exception as e:
   rows.append(dict(accession=acc,url=url,status='ERROR',error=repr(e)));print('ERROR',acc,repr(e),flush=True)
(P/'metadata/ACQUISITION.json').write_text(json.dumps(rows,indent=2))
print('METADATA_DONE',flush=True)
