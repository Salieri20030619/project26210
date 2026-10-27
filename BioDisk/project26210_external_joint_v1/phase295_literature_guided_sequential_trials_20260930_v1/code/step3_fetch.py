import subprocess,json,hashlib
from pathlib import Path
P=Path('phase295_literature_guided_sequential_trials_20260930_v1')
j=json.loads((P/'resources/TabulaTIME_tree.json').read_text());a=[]
names=['6_Application/SCINA_markergene_based_annotation.R','6_Application/Selina_Celltype_Prediction.sh','6_Application/TabulaTIME_pretrain_SELINA_meta.pkl','5_Characterization_of_cell_subtypes/NichNet.R','README.md']
for n in names:
 f=next(x for x in j['tree'] if x['path']==n);dest=P/'resources/TabulaTIME'/n;dest.parent.mkdir(exist_ok=True,parents=True)
 url='https://raw.githubusercontent.com/wanglabtongji/TabulaTiME/'+j['sha']+'/'+n
 c=subprocess.run(['curl','-L','--connect-timeout','10','--max-time','40','-fsS',url,'-o',str(dest)],capture_output=True,text=True)
 ok=False
 if c.returncode==0:
  b=dest.read_bytes();ok=hashlib.sha1(('blob '+str(len(b))+'\0').encode()+b).hexdigest()==f['sha'];assert ok
 a.append(dict(file=n,verified=ok,url=url,error=c.stderr))
for pmc in ['PMC6113240','PMC12643937']:
 dest=P/'resources'/ (pmc+'.xml');url='https://www.ebi.ac.uk/europepmc/webservices/rest/'+pmc+'/fullTextXML'
 c=subprocess.run(['curl','-L','--connect-timeout','10','--max-time','45','-fsS',url,'-o',str(dest)],capture_output=True,text=True)
 a.append(dict(file=dest.name,verified=c.returncode==0,url=url,error=c.stderr))
(P/'resources/STEP3_RESOURCE_AUDIT.json').write_text(json.dumps(a,indent=2));print(json.dumps(a,indent=2))
