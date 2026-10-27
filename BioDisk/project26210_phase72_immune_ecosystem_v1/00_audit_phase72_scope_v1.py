#!/usr/bin/env python3
import argparse,json
from pathlib import Path
import pandas as pd

META67=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase67_direct_program_replication_v1/inputs_v1/metadata_v2')
META207=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase64_large_immune_validation_v1/results_v1/GSE207493/gse207493_cell_metadata_scores_v1.tsv.gz')

def external(dataset):
    m=pd.read_csv(META67/f'{dataset}_phase67_metadata_v1.tsv.gz',sep='\t',compression='gzip',dtype=str)
    units=m[['patient_key','timepoint']].drop_duplicates()
    per_patient=units.groupby('patient_key').size()
    return m,{
        'cohort':dataset,'cells':len(m),'patients':m.patient_key.nunique(),'patient_timepoints':len(units),
        'patients_with_multiple_timepoints':int((per_patient>1).sum()),
        'timepoint_counts':m[['patient_key','timepoint']].drop_duplicates().timepoint.value_counts().to_dict(),
        'treatment_counts':m[['patient_key','timepoint','treatment']].drop_duplicates().treatment.value_counts().to_dict(),
        'response_counts':m[['patient_key','response']].drop_duplicates().response.value_counts().to_dict(),
        'broad_cell_counts':m.phase7_broad_annotation.value_counts().to_dict(),
        'fine_cell_counts':m.phase7_fine_annotation.value_counts().to_dict(),
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--outdir',required=True);a=ap.parse_args();out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
    frames=[];summaries=[]
    for d in ('GSE264586','GSE294109'):
        m,s=external(d);summaries.append(s)
        frames.append(m[['dataset','patient_key','patient_id','sample_id','timepoint','treatment','response','phase7_broad_annotation','phase7_fine_annotation']])
    g=pd.read_csv(META207,sep='\t',compression='gzip',dtype=str)
    summaries.append({'cohort':'GSE207493','cells':len(g),'patients':g['sample'].nunique(),'patient_timepoints':g['sample'].nunique(),'patients_with_multiple_timepoints':0,'timepoint_counts':{'Single':g['sample'].nunique()},'treatment_counts':{'Unknown':g['sample'].nunique()},'response_counts':{'Unknown':g['sample'].nunique()},'broad_cell_counts':g.cell_type.value_counts().to_dict(),'fine_cell_counts':{}})
    (out/'phase72_scope_audit_v1.json').write_text(json.dumps(summaries,ensure_ascii=False,indent=2),encoding='utf-8')
    flat=[]
    for s in summaries:
        flat.append({k:v for k,v in s.items() if not isinstance(v,dict)})
    pd.DataFrame(flat).to_csv(out/'phase72_scope_summary_v1.tsv',sep='\t',index=False,lineterminator='\n')
    pd.concat(frames,ignore_index=True).drop_duplicates().to_csv(out/'phase72_external_unit_annotation_counts_v1.tsv.gz',sep='\t',index=False,compression='gzip',lineterminator='\n')
    print(json.dumps({'status':'PASS','cohorts':3,'summaries':summaries},ensure_ascii=False))
if __name__=='__main__':main()
