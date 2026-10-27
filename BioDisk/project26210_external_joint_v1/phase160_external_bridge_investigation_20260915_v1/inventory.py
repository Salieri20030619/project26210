from pathlib import Path
import pandas as pd,json,subprocess,datetime
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');Q=Path('/mnt/BioDisk/wangmingyuan/project26210_phase98_public_paired_metabolomics_v1');O=P/'phase160_external_bridge_investigation_20260915_v1';O.mkdir(exist_ok=True)
gate=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True);(O/'MEMORY_GATE.log').write_text(gate.stdout+gate.stderr);assert gate.returncode==0
agreement='\n\n## 2026-09-15 用户确认C4与整体巨噬的双层分工\n用户明确同意：C4保留用于定位SPP1信号集中的胆固醇处理状态；整体巨噬细胞继续承载SPP1代谢联系分析，并授权继续调查。后续不把C4设为所有TYMP—SPP1关系的强制终点，也不因尚未证明TYMP特异作用而删除已确立C4表型。继续独立公共资料核查，保留MTA、原始正负结果、患者单位和P/FDR；不绘图、不使用客户湿实验作为公共证据。\n'
p=P/'AGENTS.md'
if '用户确认C4与整体巨噬的双层分工' not in p.read_text():p.write_text(p.read_text()+agreement)
mapping=pd.read_csv(Q/'MasterMapping_ccRCC.csv',dtype=str).fillna('')
inv=[]
for cohort in ['RC18','RC20','CPTAC','CPTAC_val']:
 rfile=Q/'downloads_v1'/('matched_tpm_'+cohort+'.csv');mfile=Q/'downloads_v1'/('matched_Harmonized_RawData_'+cohort+'.csv' if cohort.startswith('RC') else 'matched_Harmonized_Met_'+cohort+'.csv')
 rh=pd.read_csv(rfile,nrows=0).columns.tolist();mh=pd.read_csv(mfile,nrows=0).columns.tolist();mm=mapping[mapping.Dataset==cohort]
 hits=[x for x in mh if any(y in x.lower() for y in ['lactate','lactic','cholesterol','hydroxychol','desmosterol','lanosterol','pyruvate','methylthioadenosine'])]
 inv.append(dict(cohort=cohort,rna_file=str(rfile),met_file=str(mfile),metabolite_hits=hits,tumor=int((mm.TN=='Tumor').sum()),rna_genes=len(rh)-1,genes_missing=[x for x in ['TYMP','SPP1','MSMO1','HK2','LDHA','LDHB','SLC16A3','CYP27A1','LST1','AIF1','CSF1R','CD68','FCER1G','EPCAM','PTPRC'] if x not in rh],sample_mapping=mm[['RNAID','CommonID','ITHID','Identifier','TN']].head(5).to_dict('records')))
(O/'INPUT_INVENTORY.json').write_text(json.dumps(inv,indent=2))
z=json.loads((Q/'zenodo_record_11286535_v1.json').read_text())
print('ZENODO',json.dumps(z.get('metadata',{}),ensure_ascii=False)[:6500])
print('INVENTORY',json.dumps(inv,ensure_ascii=False,indent=2))
tum=mapping[mapping.TN=='Tumor']
print('ALL TUMOR CROSS DATASET DUPLICATE RNAIDS',tum[tum.RNAID.duplicated(False)][['Dataset','RNAID']].to_string(index=False))
print('MAPPING DATASETS',mapping.groupby(['Dataset','TN']).size().to_dict())

