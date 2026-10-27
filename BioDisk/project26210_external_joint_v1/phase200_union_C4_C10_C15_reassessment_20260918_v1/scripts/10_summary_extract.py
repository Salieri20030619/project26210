from common import *
import json
read=lambda p:pd.read_csv(D/p,sep='\t')
out={}
met=read('bridge/SAME_LIBRARY_GROUP_METRICS.tsv')
out['same16_composition']=met.groupby('group')[['fraction','SPP1_RNA_share']].median().to_dict('index')
a=read('core/KEY_MARKER_COMPARISONS.tsv')
out['SPP1']=a[(a.gene=='SPP1')&a.contrast.isin(['all_U','all_C4','all_U_common','all_C4_common','all_U_without_C10_dominant_patient','all_C4_without_C10_dominant_patient','tumor_U','tumor_C4'])][['contrast','n','median_delta','positive','P','FDR']].to_dict('records')
p=read('functional/ALL_PROGRAM_RESULTS.tsv.gz')
terms=['Cholesterol metabolism','Lysosome','Phagosome','Antigen processing and presentation','Glycolysis','E2F Targets','G2-M Checkpoint']
out['functional']=p[p.contrast.isin(['all_U','all_C4'])&(p.version=='without_SPP1_TYMP')&p.term.isin(terms)][['contrast','term','n','median','positive','P','FDR']].to_dict('records')
g=read('gsea/FOCUS_GSEA_RESULTS.tsv')
out['gsea']=g[g.model.isin(['all_U_common','all_C4_common'])&(g.version=='complete')][['model','ID','NES','pvalue','p.adjust']].to_dict('records')
b=read('bridge/ALL_FIXED_CANDIDATE_BRIDGES.tsv')
sel=['ligand:IGFBP3','program:Cholesterol Homeostasis','program:Glycolysis','gene:SLC16A3','gene:LDHA','direct_MTA_response_RNA_proxy','polyamine_synthesis_core_RNA','MTAP_clearance_RNA']
out['candidate_bridges']=b[b.candidate.isin(sel)&(b.adjustment=='study_adjusted')][['group','candidate','a_rho','a_FDR','XY_rho','XY_P','XY_FDR','b_partial_rho','b_P','b_FDR','ab_low','ab_high']].to_dict('records')
d=read('bridge/TYMP_DIRECT_RESULTS.tsv')
out['TYMP_direct']=d[(d.group.isin(['U','C4']))&(d.exposure=='TYMP_pseudobulk')&(d.target=='fraction')].to_dict('records')
out['TYMP_group_SPP1']=read('bridge/MATCHED_GROUP_SPP1_DIRECT_TYMP.tsv').query("population == 'common_eligible'").to_dict('records')
out['MTA_response']=read('receiver/MTA_RESPONSE_SUMMARY.tsv').query("scope == 'tumor'").to_dict('records')
r=read('receiver/WITHIN_CELL_SUMMARY.tsv')
out['receiver_key']=r[r.group.isin(['U','C4','WholeMac'])&(r.adjustment=='depth')&(r.module.isin(['same_library_lactate','STAT3_frozen_collectri','RNA']))][['module','group','edge','n','median','positive','P','FDR']].to_dict('records')
out['stat3_activity']=read('receiver/STAT3_ACTIVITY_SUMMARY.tsv').to_dict('records')
t=read('regulons/ALL_GROUP_AUC_SUMMARY.tsv')
out['threeTF']=t[t.group.isin(['U','C4'])&(t.scope=='tumor')&t.TF.isin(['JUN','JUNB','CEBPD'])][['group','TF','n','median','positive','P','FDR']].to_dict('records')
out['coverage']=read('trajectory/COVERAGE_SUMMARY.tsv').query("group in ['U','C4'] and pseudotime == 'pt_primary'").to_dict('records') if (D/'trajectory/COVERAGE_SUMMARY.tsv').exists() else []
(D/'qc/RESULTS_FOR_REVIEW.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
for key in ['same16_composition','SPP1','functional','TYMP_direct','MTA_response','stat3_activity','coverage']:
 print(key,json.dumps(out[key],ensure_ascii=False))
