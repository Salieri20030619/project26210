from pathlib import Path
import json,csv,hashlib,datetime
P=Path(__file__).resolve().parents[1];R=P.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
st=json.loads((P/'qc/STATUS.json').read_text())
items=[x for x in csv.DictReader((R/'RESULTS_MANIFEST.tsv').open(),delimiter='\t') if x['stage']=='Phase295' and x['result_id']!='phase295_sequential_plan']
checks={}
checks['plan_unchanged']=sha(P/'CURRENT_STATE_AND_SEQUENTIAL_PLAN.md')==st['plan_sha256']
for n in ['step1_network','step1_validate','step3_alignment']:checks[n+'_exit0']=(P/'logs'/ (n+'.exit')).read_text().strip()=='0'
checks['step2_no_prediction']=json.loads((P/'step2_immune_phenotype/COMPLETION.json').read_text())['predictions']==0
with (P/'step3_reference_alignment/REFERENCE_GENE_CONCORDANCE.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
checks['27_genes_and7patients']=len(rows)==27 and all(int(x['patient_n'])==7 for x in rows)
checks['manifest7_unique']=len(items)==7 and len(set(x['result_id'] for x in items))==7
checks['manifest_hashes_match']=all(x['notes']=='SHA256='+sha(Path(x['absolute_path'])) for x in items)
names=['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv']
checks['all_root_entries']=all((R/n).read_text().startswith('## Phase295终态') for n in names if n.endswith('.md'))
checks['histories_preserved']=all((R/n).read_bytes().endswith((P/'root_record_backups/before_final'/n).read_bytes()) for n in names if n.endswith('.md'))
checks['manifest_history_preserved']=(R/'RESULTS_MANIFEST.tsv').read_bytes().startswith((P/'root_record_backups/before_final/RESULTS_MANIFEST.tsv').read_bytes())
checks['partial_explicit']=(P/'resources/RCC2021_supplement.pdf.partial').exists() and not (P/'resources/RCC2021_supplement.pdf').exists()
assert all(checks.values()),checks
out=dict(status='PASS',at=datetime.datetime.now().astimezone().isoformat(),checks={k:bool(v) for k,v in checks.items()},checks_passed=len(checks),recovery='Initial finalization finished outputs/root updates and all assertions but JSON serialization rejected numpy.bool. Recovery verifies actual saved files and writes plain bool. Initial finalize.exit1 preserved; no statistics rerun.',artifacts=[dict(path=x['absolute_path'],sha256=sha(Path(x['absolute_path']))) for x in items],root_hashes={n:sha(R/n) for n in names})
with (P/'qc/FINAL_VALIDATION.json').open('x') as f:json.dump(out,f,indent=2)
st['final_validation']='PASS';st['finalization_recovered']=True;(P/'qc/STATUS.json').write_text(json.dumps(st,indent=2))
print(json.dumps(dict(status='PASS',checks_passed=len(checks),report=items[0]['absolute_path']),ensure_ascii=False))
