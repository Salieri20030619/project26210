from pathlib import Path
import json,openpyxl,zipfile,subprocess
P=Path(__file__).resolve().parents[1]
subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True)
rows=[]
for f in sorted((P/'downloads').glob('*.xlsx')):
 if f.name=='41467_2022_29366_MOESM9_ESM.xlsx':continue
 with zipfile.ZipFile(f) as z:assert z.testzip() is None
 w=openpyxl.load_workbook(f,read_only=True,data_only=True)
 row={'file':f.name,'sheets':[]}
 for ws in w:
  sample=list(ws.iter_rows(min_row=1,max_row=3,values_only=True));row['sheets'].append({'name':ws.title,'rows':ws.max_row,'cols':ws.max_column,'first3':sample})
 w.close();rows.append(row)
print(json.dumps(rows,indent=2,default=str),flush=True)
(P/'qc/INPUT_SCHEMA_v1.json').write_text(json.dumps(rows,indent=2,default=str))
