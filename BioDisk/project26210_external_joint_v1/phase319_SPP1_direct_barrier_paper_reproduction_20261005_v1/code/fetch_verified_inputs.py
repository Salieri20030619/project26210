from pathlib import Path
import requests,json,hashlib
P=Path(__file__).resolve().parents[1]
for a in json.loads((P/'resources/DOWNLOAD_MANIFEST.json').read_text()):
 f=P/'downloads'/a['file']
 if f.exists():
  assert hashlib.sha256(f.read_bytes()).hexdigest()==a['sha256'];continue
 q=requests.get(a['url'],timeout=(20,90));q.raise_for_status();assert hashlib.sha256(q.content).hexdigest()==a['sha256'];f.write_bytes(q.content)
