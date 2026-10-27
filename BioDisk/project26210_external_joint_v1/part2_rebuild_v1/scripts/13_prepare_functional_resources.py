from pathlib import Path
import urllib.request,json,hashlib,datetime
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');O=R/'06_functional_programs/run_v1';O.mkdir(exist_ok=False);D=O/'resources';D.mkdir()
plan={'status':'FROZEN_BEFORE_RESULTS','libraries':['MSigDB_Hallmark_2020','KEGG_2021_Human','Reactome_2022'],'gene_set_size':[10,500],'primary':'C4 patient-paired mean log2(CPM+1) gene-set difference; per-library BH; complete and no TYMP/SPP1 versions','secondary':'All12 macrophage candidates same analysis; ORA on eligible nominal P<0.05 abs median gene effect>=0.25; up/down directions corrected together per cluster/library/version. Raw counts gene universe restricted to eligible measured genes.','n_unit':'patient; n<3 no test, n<10 exploratory; Wilcoxon approx continuity correction; n>=10 nominal P<0.05 acceptable exploratory support, full FDR retained','interpretation':'Program transcript expression, not metabolic flux or protein activity. Same-data discovery conditional evidence.','figures':False,'focus_categories':'lipid/cholesterol/fatty acid; lysosome/phagosome; inflammatory/TNF/interferon; antigen presentation; hypoxia; angiogenesis/vascular'}
(O/'PLAN.json').write_text(json.dumps(plan,indent=2))
reg=[]
for lib in plan['libraries']:
 url='https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName='+lib
 data=urllib.request.urlopen(url,timeout=90).read();txt=data.decode();assert len(txt.splitlines())>40 and '<html' not in txt.lower()
 (D/(lib+'.gmt')).write_bytes(data);reg.append({'library':lib,'url':url,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()});print(lib,len(txt.splitlines()),flush=True)
(O/'RESOURCE_REGISTRY.json').write_text(json.dumps(reg,indent=2))
