import numpy as np,pandas as pd,json,hashlib,subprocess,urllib.request
from pathlib import Path
from scipy import sparse
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');P=R/'phase295_literature_guided_sequential_trials_20260930_v1';O=P/'step1_receiver_network'
def gate():
 with open(P/'logs/memory_gate.log','a') as f:subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],stdout=f,stderr=f,check=True)
protocol=dict(contrast='Frozen activated versus nonactivated cells within original Pericyte_SMC broad identity, purified stromal pool, Tumor only; descriptive state contrast, not treatment effect',min_cells_each_state=20,library_selection='largest minimum state count per patient, then total cells, then lexical library; no expression outcomes',min_macro_cells=100,expression_fraction=.05,patient_presence_fraction=.5,target_median_log2_cpm_delta=.25,target_positive_patient_fraction=.7,target_inference='n<10: no P/FDR; descriptive consistent effect target set; n>=10 would require paired patient test raw P<.05',ranking='NicheNet official all receiver-receptor-eligible ligands; focused subset expressed mature macrophage; Pearson primary, AUPR secondary; C4 expression separate annotation only',exclusions='all Phase286 identity/activation and all Phase289 frozen panels excluded from background and target, avoids direct marker reuse; clustering/state selection dependence remains',prior='official human Zenodo7074291 v2, verify existing resource checksum if API accessible',claim='candidate ligand-to-observed-stromal-state target network; neither SPP1 direct causality nor TYMP mediated immune exclusion')
(O/'PROTOCOL_FROZEN.json').write_text(json.dumps(protocol,ensure_ascii=False,indent=2))
gate()
m=pd.read_csv(R/'phase289_state_source_boundary_process_20260929_v1/tables/STROMAL_STATE_MEMBERSHIP.tsv.gz',sep='\t')
print('Broad labels',m.cell_type_broad_v1.value_counts().to_dict(),flush=True)
# Original mural identity is fixed, not relabeled as CAF.
mask=m.stromal_pool & m.source_tissue.eq('Tumor') & m.cell_type_broad_v1.str.contains('Pericyte|SMC',case=False,regex=True)
rows=[]
for (pat,lib),g in m[mask].groupby(['curated_patient_id','sample_id']):
 a=int(g.activated_candidate.sum());b=len(g)-a
 rows.append(dict(patient=pat,library=lib,study=g.study.iloc[0],activated=a,reference=b,min_state=min(a,b),total=len(g),eligible=a>=20 and b>=20))
q=pd.DataFrame(rows);q['selected']=False
pick=q[q.eligible].sort_values(['patient','min_state','total','library'],ascending=[True,False,False,True]).drop_duplicates('patient').index
q.loc[pick,'selected']=True;q.to_csv(O/'PATIENT_LIBRARY_QUALIFICATION.tsv',sep='\t',index=False)
sel=q[q.selected].sort_values('patient')
assert len(sel)>1
assert len(sel)<10,'Protocol requires patient DE branch for n>=10; stop before calculation'
X=sparse.load_npz(R/'phase287_stromal_refinement_image_morphology_20260929_v1/inputs/STROMAL_COUNTS.npz').astype(np.float64)
genes=np.load(R/'phase287_stromal_refinement_image_morphology_20260929_v1/inputs/GENES.npy',allow_pickle=True).astype(str)
im=pd.read_csv(R/'phase287_stromal_refinement_image_morphology_20260929_v1/inputs/STROMAL_INPUT_METADATA.tsv.gz',sep='\t')
assert np.array_equal(m.cell_id,im.cell_id) and X.shape==(len(m),len(genes)) and len(set(genes))==len(genes)
ds=[];expr=[];pb=[]
for r in sel.itertuples():
 gate()
 base=mask&m.curated_patient_id.eq(r.patient)&m.sample_id.eq(r.library)
 vals=[]
 for state,keep in [('activated',base&m.activated_candidate),('reference',base&~m.activated_candidate)]:
  z=X[np.flatnonzero(keep)];counts=np.asarray(z.sum(axis=0)).ravel();cpm=1e6*counts/counts.sum();vals.append(np.log2(cpm+1));expr.append((np.asarray((z>0).sum(axis=0)).ravel()/z.shape[0])>=.05)
  pb.append(pd.DataFrame(dict(patient=r.patient,library=r.library,state=state,gene=genes,counts=counts,cpm=cpm)))
 ds.append(vals[0]-vals[1])
D=np.array(ds);ex=np.array(expr).reshape(len(sel),2,-1).any(axis=1).mean(axis=0)>=.5
panels=json.loads((R/'phase289_state_source_boundary_process_20260929_v1/PANELS_FROZEN.json').read_text())
identity=json.loads((R/'phase286_CAF_identity_boundary_audit_20260929_v1/resources/IDENTITY_PANELS_FROZEN.json').read_text())
def flatten(x):
 if isinstance(x,dict):return sum([flatten(v) for v in x.values()],[])
 if isinstance(x,list):return sum([flatten(v) for v in x],[])
 return [x] if isinstance(x,str) else []
excluded=set(flatten(panels)+flatten(identity))
out=pd.DataFrame(dict(gene=genes,median_log2_cpm_delta=np.median(D,axis=0),positive_patient_fraction=(D>0).mean(axis=0),receiver_expressed=ex,excluded_definition_gene=np.isin(genes,list(excluded))))
out['background']=out.receiver_expressed&~out.excluded_definition_gene
out['target']=out.background&(out.median_log2_cpm_delta>=.25)&(out.positive_patient_fraction>=.7)
out.to_csv(O/'RECEIVER_GENE_INPUT.tsv',sep='\t',index=False)
pd.concat(pb).to_csv(O/'RECEIVER_PSEUDOBULK.tsv.gz',sep='\t',index=False)
pd.DataFrame(D,index=sel.patient,columns=genes).rename_axis('patient').to_csv(O/'PATIENT_CONTRASTS.tsv.gz',sep='\t')
gate()
mm=pd.read_csv(R/'phase200_union_C4_C10_C15_reassessment_20260918_v1/inputs/MACRO_METADATA.tsv.gz',sep='\t')
MX=sparse.load_npz(R/'phase200_union_C4_C10_C15_reassessment_20260918_v1/inputs/MACRO_COUNTS.npz')
mg=pd.read_csv(R/'phase200_union_C4_C10_C15_reassessment_20260918_v1/inputs/GENES.tsv',sep='\t').gene.astype(str).values
assert MX.shape==(len(mm),len(mg))
exp=[];c4=[];senderqual=[]
for r in sel.itertuples():
 idx=np.flatnonzero(mm.curated_patient_id.eq(r.patient)&mm.sample_id.eq(r.library)&mm.source_tissue.eq('Tumor'))
 senderqual.append(dict(patient=r.patient,library=r.library,macrophages=len(idx),eligible=len(idx)>=100))
 if len(idx)>=100:exp.append(np.asarray((MX[idx]>0).sum(axis=0)).ravel()/len(idx)>=.05)
 ci=idx[mm.iloc[idx].candidate.astype(str).eq('4').values]
 if len(ci)>=20:c4.append(np.asarray((MX[ci]>0).sum(axis=0)).ravel()/len(ci)>=.05)
pd.DataFrame(senderqual).to_csv(O/'SENDER_QUALIFICATION.tsv',sep='\t',index=False)
assert exp
pd.DataFrame(dict(gene=mg,macrophage_expressed=np.mean(exp,axis=0)>=.5,c4_expressed=np.mean(c4,axis=0)>=.5 if c4 else False)).to_csv(O/'SENDER_GENE_INPUT.tsv',sep='\t',index=False)
resource=Path('/mnt/BioDisk/ludingyu/project26210/external/nichenet')
audit=[]
try:
 raw=urllib.request.urlopen('https://zenodo.org/api/records/7074291',timeout=45).read();(P/'resources/zenodo7074291.json').write_bytes(raw);official={f['key']:f for f in json.loads(raw)['files']}
except Exception as e:
 official={};(P/'resources/zenodo7074291.error.txt').write_text(repr(e))
for name in ['ligand_target_matrix_nsga2r_final.rds','lr_network_human_21122021.rds','weighted_networks_nsga2r_final.rds']:
 f=resource/name;data=f.read_bytes();md5=hashlib.md5(data).hexdigest();ref=official.get(name,{})
 audit.append(dict(path=str(f),bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),md5=md5,official_checksum=ref.get('checksum',''),checksum_match=ref.get('checksum')=='md5:'+md5,source='https://zenodo.org/records/7074291'))
pd.DataFrame(audit).to_csv(O/'PRIOR_RESOURCE_AUDIT.tsv',sep='\t',index=False)
summary=dict(patients=len(sel),patients_list=sel.patient.tolist(),background=int(out.background.sum()),targets=int(out.target.sum()),sender_patients=len(exp),c4_sender_patients=len(c4),resource_verified=all(x['checksum_match'] for x in audit),P_FDR='not calculated: n<10')
(O/'INPUT_COMPLETE.json').write_text(json.dumps(summary,indent=2));print(summary,flush=True)
