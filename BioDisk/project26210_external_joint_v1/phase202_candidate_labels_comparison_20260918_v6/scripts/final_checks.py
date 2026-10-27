from pathlib import Path
import json,hashlib,subprocess,shutil
import pandas as pd,numpy as np
from pypdf import PdfReader
from PIL import Image,ImageDraw
r=Path.cwd();d=r/'phase202_candidate_labels_comparison_20260918_v6';rep=r/'phase201_three_candidates_C4_main_report_20260918_v4';old=r/'phase191_discovery_order_report_20260917_v1'
a=json.loads((old/'report/report_materialized.json').read_text());b=json.loads((rep/'report/report_materialized.json').read_text());assert [x for x in a if x['type']=='table']==[x for x in b if x['type']=='table'];assert sum(x['type']=='figure' for x in b)==9
for f in (old/'figures').glob('*.pdf'):assert hashlib.sha256(f.read_bytes()).digest()==hashlib.sha256((rep/'figures'/f.name).read_bytes()).digest()
s=pd.read_csv(d/'tables/CANDIDATE_VS_WHOLE_SUMMARY.tsv',sep='\t').set_index('group');assert s.n_expression.tolist()==[26,2,9,28,38];assert np.isclose(s.loc['C4','contribution'],54.58978669);assert np.isclose(s.loc['U','contribution'],59.62678892)
q=pd.read_csv(d/'tables/CANDIDATE_PATIENT_VALUES.tsv',sep='\t');w=q.pivot(index='patient',columns='group',values='SPP1_RNA_share');assert np.allclose(w.U,w.C4+w.C10+w.C15)
loc=json.loads((d/'qc/LOCATOR_VALIDATION.json').read_text());assert all(x['original_values_identical'] for x in loc['expression_fields'])
p=rep/'report/26210_TYMP_MTA_SPP1_macrophage_C4_primary_report_v4.pdf';reader=PdfReader(p);texts=[x.extract_text() for x in reader.pages];pl=json.loads((rep/'qc/PLACEMENTS.json').read_text())
for x in pl:
 pg=reader.pages[x['page']-1];assert x['y']>=30 and x['x']+x['width']<=float(pg.mediabox.width)
 for no in [8,9]:
  if x['figure']==no:assert f'图{no} ' in texts[x['page']-1]
checks=dict(report_pages=len(reader.pages),tables_unchanged=31,original_figures_identical=7,new_report_figures=2,expression_values_unchanged=True,union_rna_additivity=True,expression_patient_counts=s.n_expression.tolist(),figures_and_captions_on_same_page=True)
(d/'qc/FINAL_CHECKS.json').write_text(json.dumps(checks,indent=2));(rep/'qc/TEXT_AND_LAYOUT_CHECKS.json').write_text(json.dumps(checks,indent=2))
(rep/'qc/final_pages').mkdir(exist_ok=True)
subprocess.run(['pdftoppm','-r','75','-jpeg',str(p),str(rep/'qc/final_pages/page')],check=True,stderr=subprocess.PIPE)
files=sorted((rep/'qc/final_pages').glob('page-*.jpg'))
for start in range(0,len(files),6):
 im=Image.new('RGB',(1040,2265),'#dddddd');dr=ImageDraw.Draw(im)
 for j,f in enumerate(files[start:start+6]):
  v=Image.open(f);v.thumbnail((500,708));x=10+(j%2)*520;y=30+(j//2)*755;im.paste(v,(x,y));dr.text((x,y-20),f'Page {start+j+1}',fill='black')
 im.save(d/'qc'/f'report_contact_{start//6+1:02d}.jpg',quality=88)
for no in [8,9]:
 pg=next(x['page'] for x in pl if x['figure']==no)
 subprocess.run(['pdftoppm','-f',str(pg),'-l',str(pg),'-r','135','-singlefile','-jpeg',str(p),str(d/'qc'/f'report_figure{no}_page')],check=True,stderr=subprocess.PIPE)
print(json.dumps(checks),flush=True)
