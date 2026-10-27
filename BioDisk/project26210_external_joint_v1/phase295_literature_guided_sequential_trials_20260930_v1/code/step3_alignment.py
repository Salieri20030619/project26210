from pathlib import Path
import xml.etree.ElementTree as E
import pandas as pd,numpy as np,json,subprocess,hashlib
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase295_literature_guided_sequential_trials_20260930_v1');O=P/'step3_reference_alignment'
protocol=dict(reference='Tomko2018 Table1, PMID30154546, DOI10.1038/s41598-018-31126-w',endpoint='gene-level direction concordance, no composite claimed to reproduce DDR1 ssGSEA; no inferred collagen geometry',reference_scope='27 normal-versus-aligned-IDC differential/trending proteins; different from 24 alignment-correlated proteins',direction='LAMB2 and OGN down, remaining25 up per original footnote',analysis='reuse frozen 7 patient within-library state pseudobulks; report median log2(CPM+1) contrast, sign consistency, all27 and nondefinition subset; no P/FDR with n7',selection='all27 authentic genes, retain original source P including .06/.10; no external-significance relabel',annotation='SELINA downloaded metadata has45 broad classes, lacks Fibro_CTHRC1/Macro_SLPI; cannot transfer those states with this model')
(O/'PROTOCOL_FROZEN.json').write_text(json.dumps(protocol,indent=2))
t=E.parse(P/'resources/PMC6113240.xml');tab=next(x for x in t.iter('table-wrap') if x.attrib.get('id')=='Tab1')
rows=[]
for tr in tab.iter('tr'):
 v=[''.join(td.itertext()).strip() for td in tr if td.tag=='td']
 if len(v)==4:rows.append(dict(protein=v[0],gene=v[1],source_fold_change_text=v[2],source_P=float(v[3]),source_direction=-1 if v[1] in ['LAMB2','OGN'] else 1))
ref=pd.DataFrame(rows);assert len(ref)==27 and ref.gene.is_unique
ref.to_csv(O/'TOMKO2018_TABLE1_27GENES_FROZEN.tsv',sep='\t',index=False)
(O/'TOMKO_TABLE1_CAPTION_AND_FOOTNOTE.txt').write_text(' '.join(tab.itertext()))
with open(P/'logs/memory_gate.log','a') as f:subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],stdout=f,check=True)
d=pd.read_csv(P/'step1_receiver_network/PATIENT_CONTRASTS.tsv.gz',sep='\t',index_col=0)
g=pd.read_csv(P/'step1_receiver_network/RECEIVER_GENE_INPUT.tsv',sep='\t').set_index('gene')
result=[];long=[]
for r in ref.itertuples():
 present=r.gene in d.columns
 vals=d[r.gene].values if present else np.array([])
 result.append(dict(gene=r.gene,source_direction=r.source_direction,source_P=r.source_P,covered=present,patient_n=len(vals),median_log2_cpm_delta=float(np.median(vals)) if present else None,concordant_patients=int((vals*r.source_direction>0).sum()),discordant_patients=int((vals*r.source_direction<0).sum()),zero_patients=int((vals==0).sum()),receiver_expressed=bool(g.loc[r.gene,'receiver_expressed']) if present else False,definition_gene=bool(g.loc[r.gene,'excluded_definition_gene']) if present else False))
 if present:
  for pat,v in zip(d.index,vals):long.append(dict(patient=pat,gene=r.gene,delta=v,source_direction=r.source_direction,direction_concordant=v*r.source_direction>0))
res=pd.DataFrame(result);res['median_concordant']=res.median_log2_cpm_delta*res.source_direction>0
res.to_csv(O/'REFERENCE_GENE_CONCORDANCE.tsv',sep='\t',index=False)
lng=pd.DataFrame(long);lng.to_csv(O/'PATIENT_GENE_CONCORDANCE.tsv',sep='\t',index=False)
outs=[]
for name,genes in [('all27',ref.gene.tolist()),('nondefinition',[x for x in ref.gene if not bool(g.loc[x,'excluded_definition_gene'])])]:
 for pat,z in lng[lng.gene.isin(genes)].groupby('patient'):
  outs.append(dict(set=name,patient=pat,genes=len(z),direction_concordant=int(z.direction_concordant.sum()),fraction=float(z.direction_concordant.mean())))
pd.DataFrame(outs).to_csv(O/'PATIENT_DIRECTION_SUMMARY.tsv',sep='\t',index=False)
# Independent sign/direction checks, no selected positive examples.
assert set(res.gene)==set(ref.gene)
assert all(int((d[r.gene]*r.source_direction>0).sum())==r.concordant_patients for r in res.itertuples())
assert int((ref.source_direction==1).sum())==25
summary=dict(status='COMPUTED_DESCRIPTIVE_REFERENCE_CONCORDANCE',genes=27,covered=int(res.covered.sum()),patients=len(d),median_direction_concordant=int(res.median_concordant.sum()),nondefinition_genes=int((~res.definition_gene).sum()),nondefinition_direction_concordant=int(res.loc[~res.definition_gene,'median_concordant'].sum()),at_least_5_of_7_concordant=int((res.concordant_patients>=5).sum()),source_nominal_P_not_below_05=ref[ref.source_P>=.05].gene.tolist(),validation='27 rows/25up2down/symbol coverage/perpatient signs recomputed',inference='descriptive RNA concordance with source proteomic contrast; not measured alignment or DDR1 score reproduction')
(O/'CONCORDANCE_COMPLETE.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2));print(res.to_string(index=False))
