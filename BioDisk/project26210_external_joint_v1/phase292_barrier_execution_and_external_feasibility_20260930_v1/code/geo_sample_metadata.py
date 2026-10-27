from pathlib import Path
import requests,json,re,datetime
P=Path(__file__).resolve().parents[1]
u='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204809&targ=all&form=text&view=full'
r=requests.get(u,timeout=20);r.raise_for_status();text=r.text
(P/'resources/GSE204809_FULL_METADATA.soft.txt').write_text(text)
rows=[]
for block in re.split(r'\^SAMPLE = ',text)[1:]:
 row=dict(sample=block.splitlines()[0])
 for key in ['title','source_name_ch1','characteristics_ch1','treatment_protocol_ch1','data_processing']:
  row[key]=[s.split(' = ',1)[1] for s in block.splitlines() if s.startswith('!Sample_'+key+' = ')]
 rows.append(row)
(P/'resources/GSE204809_SAMPLE_DESIGN.json').write_text(json.dumps(rows,indent=2))
print(json.dumps([dict(sample=x['sample'],title=x['title'],characteristics=x['characteristics_ch1']) for x in rows],indent=2))
print('METADATA_PASS',len(rows))
