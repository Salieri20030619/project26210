from pathlib import Path
import openpyxl,pandas as pd,subprocess,json
P=Path(__file__).resolve().parents[1]
subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True)
w=pd.ExcelFile(P/'resources/author_subtype_markers.xlsx');print('SHEETS',w.sheet_names)
for n in w.sheet_names:
 d=pd.read_excel(w,sheet_name=n);print('SHEET',n,'DIM',d.shape,'COLUMNS',list(d.columns))
 for c in d.select_dtypes('object'):
  v=d[c].dropna().astype(str).unique()
  if len(v)<100:print(c,list(v))
 print(d.head(2).to_json(orient='records'))
for n in ['author_stromal.R','author_myeloid.R','author_function_score.R']:
 lines=(P/'resources'/n).read_text().splitlines();print('CODE',n)
 for i,l in enumerate(lines):
  if any(k in l for k in ['top_n','AddModuleScore','FindAllMarkers','avg_log','25','SPP1','FAP']):print(i+1,l)
