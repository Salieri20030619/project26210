from pathlib import Path
import fitz,json,collections
r=Path.cwd();d=r/'phase206_all_mainfigures_A4_20260918_v10';rows=[]
for n in [1,2,3]:
 p=fitz.open(d/'figures'/f'Main_Figure_{n}_A4_v10.pdf')[0];texts=[]
 for b in p.get_text('dict')['blocks']:
  for l in b.get('lines',[]):
   for s in l['spans']:
    if s['text'].strip():texts.append(dict(text=s['text'],size=s['size'],bbox=s['bbox']))
 small=[s for s in texts if s['size']<5.95];outside=[s for s in texts if s['bbox'][0]<25 or s['bbox'][2]>570 or s['bbox'][1]<20 or s['bbox'][3]>822]
 rec=dict(figure=n,page_mm=[p.rect.width/72*25.4,p.rect.height/72*25.4],min_font=min(z['size'] for z in texts),small=small,outside_margin=outside);rows.append(rec);print(n,'minfont',rec['min_font'],'small',len(small),'outside',len(outside),flush=True)
(d/'qc/FONT_PAGE_AUDIT.json').write_text(json.dumps(rows,indent=2))
