from pathlib import Path
import json,hashlib,datetime,subprocess
P=Path(__file__).resolve().parents[1];R=P.parent
def sha(f):
 h=hashlib.sha256()
 with open(f,'rb') as s:
  for b in iter(lambda:s.read(1048576),b''):h.update(b)
 return h.hexdigest()
def dump(n,o):(P/n).write_text(json.dumps(o,ensure_ascii=False,indent=2))
identity=R/'phase286_CAF_identity_boundary_audit_20260929_v1/resources/IDENTITY_PANELS_FROZEN.json'
old=R/'phase289_state_source_boundary_process_20260929_v1/PANELS_FROZEN.json'
a=json.loads(identity.read_text());b=json.loads(old.read_text())
excluded={'TYMP','SPP1'}|set(sum(a.values(),[]))|set(sum([b[k] for k in ['CAF_activation','transition_RCC','perivascular','fibro_identity']],[]))
c=json.loads((P/'PANELS_CANDIDATE.json').read_text())
f={k:[g for g in v if g not in excluded] for k,v in c.items()}
assert all(len(v)>=4 for v in f.values())
dump('PANELS_FROZEN.json',f)
dump('resources/PANEL_EXCLUSIONS.json',{k:dict(removed=[g for g in c[k] if g in excluded],retained=f[k],overlap_phase289={o:sorted(set(f[k])&set(v)) for o,v in b.items() if set(f[k])&set(v)}) for k in f})
paths=[identity,old,R/'phase289_state_source_boundary_process_20260929_v1/tables/STROMAL_STATE_MEMBERSHIP.tsv.gz',R/'phase287_stromal_refinement_image_morphology_20260929_v1/inputs/STROMAL_COUNTS.npz',R/'phase287_stromal_refinement_image_morphology_20260929_v1/inputs/GENES.npy',R/'phase287_stromal_refinement_image_morphology_20260929_v1/inputs/STROMAL_INPUT_METADATA.tsv.gz',R/'phase289_state_source_boundary_process_20260929_v1/tables/MATCHED_SOURCE_RECEIVER_INPUTS.tsv',R/'phase157_ordered_bridge_reassessment_20260915_v1/02_inputs/GROUP_INDEX.tsv',R/'phase157_ordered_bridge_reassessment_20260915_v1/02_inputs/GROUP_PSEUDOBULK.npz']
assert all(x.is_file() for x in paths)
dump('qc/INPUT_MANIFEST.json',[dict(path=str(x),bytes=x.stat().st_size,sha256=sha(x)) for x in paths])
dump('FREEZE.json',dict(at=datetime.datetime.now().astimezone().isoformat(),protocol_sha256=sha(P/'PLAN_AND_PROTOCOL.md'),panels_sha256=sha(P/'PANELS_FROZEN.json'),outcomes_examined=False))
print(json.dumps(f));print('PREFLIGHT_FROZEN')
