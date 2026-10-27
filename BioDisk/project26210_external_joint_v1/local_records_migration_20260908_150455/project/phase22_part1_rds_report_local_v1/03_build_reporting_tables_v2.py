#!/usr/bin/env python3
from __future__ import annotations
import json,os
from datetime import datetime
from pathlib import Path
import numpy as np,pandas as pd

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/"phase14_patient_program_validation_local_v1"/"patient_analysis_final_cnv_v6"
OUT=HERE/"reporting_tables_v2"
if OUT.exists():raise FileExistsError(OUT)
OUT.mkdir()

def atomic_tsv(frame,path):
 part=Path(str(path)+".part");frame.to_csv(part,sep="\t",index=False);os.replace(part,path)
def sign_label(value):
 if not np.isfinite(value):return "not_available"
 return "positive" if value>0 else "negative" if value<0 else "zero"

assoc=pd.read_csv(SOURCE/"within_dataset_patient_associations_v6.tsv",sep="\t")
assoc["n_independent_units"]=pd.to_numeric(assoc["n_patients"],errors="coerce")
assoc["reporting_tier"]=np.select(
 [assoc["status"].ne("TESTED"),assoc["n_independent_units"].lt(10)],
 ["NOT_TESTED","DIRECTION_ONLY_LT10"],default="FULL_STATISTICS_GE10")
assoc["observed_direction"]=[sign_label(x) for x in pd.to_numeric(assoc["spearman_rho"],errors="coerce")]
assoc["direction_summary"]=np.where(
 assoc["reporting_tier"].eq("NOT_TESTED"),"not_tested",
 np.where(assoc["direction_matches_expected"].astype(str).str.lower().eq("true"),"matches_expected","opposes_expected"))
fdr=pd.to_numeric(assoc["p_fdr_bh_within_dataset_method"],errors="coerce")
assoc["formal_result"]=np.where(
 assoc["reporting_tier"].eq("FULL_STATISTICS_GE10"),
 np.where(fdr.lt(.05),"FDR_SUPPORTED","NOT_FDR_SUPPORTED"),assoc["direction_summary"])
atomic_tsv(assoc,OUT/"association_reporting_v2.tsv")

summary=assoc.groupby(["relation","method"],as_index=False).agg(
 grid_rows=("dataset","size"),
 full_statistics_rows=("reporting_tier",lambda x:int((x=="FULL_STATISTICS_GE10").sum())),
 direction_only_rows=("reporting_tier",lambda x:int((x=="DIRECTION_ONLY_LT10").sum())),
 not_tested_rows=("reporting_tier",lambda x:int((x=="NOT_TESTED").sum())),
 expected_direction_rows=("direction_summary",lambda x:int((x=="matches_expected").sum())),
 opposite_direction_rows=("direction_summary",lambda x:int((x=="opposes_expected").sum())),
 fdr_supported_rows=("formal_result",lambda x:int((x=="FDR_SUPPORTED").sum())))
atomic_tsv(summary,OUT/"association_relation_summary_v2.tsv")

longitudinal=pd.read_csv(SOURCE/"GSE264586_paired_longitudinal_v6.tsv",sep="\t")
longitudinal["n_independent_units"]=pd.to_numeric(longitudinal["n_pairs"],errors="coerce")
longitudinal["reporting_tier"]=np.where(longitudinal["n_independent_units"].ge(10),"FULL_STATISTICS_GE10","DIRECTION_ONLY_LT10")
longitudinal["observed_direction"]=[sign_label(x) for x in pd.to_numeric(longitudinal["median_change"],errors="coerce")]
atomic_tsv(longitudinal,OUT/"longitudinal_reporting_v2.tsv")

response=pd.read_csv(SOURCE/"GSE294109_response_stratified_v6.tsv",sep="\t")
response["minimum_group_units"]=response[["n_benefit_patients","n_no_effect_patients"]].apply(pd.to_numeric,errors="coerce").min(axis=1)
response["reporting_tier"]=np.where(response["minimum_group_units"].ge(10),"FULL_STATISTICS_GE10","DIRECTION_ONLY_LT10")
response["observed_direction"]=[sign_label(x) for x in pd.to_numeric(response["median_difference_benefit_minus_no_effect"],errors="coerce")]
atomic_tsv(response,OUT/"response_reporting_v2.tsv")

f1=pd.read_csv(SOURCE/"malignant_F1_patient_reproducibility_v6.tsv",sep="\t")
f1["n_independent_units"]=pd.to_numeric(f1["paired_patients"],errors="coerce")
f1["reporting_tier"]=np.where(f1["n_independent_units"].ge(10),"FULL_STATISTICS_GE10","DIRECTION_ONLY_LT10")
f1["observed_direction"]=[sign_label(x) for x in pd.to_numeric(f1["median_malignant_minus_normal"],errors="coerce")]
atomic_tsv(f1,OUT/"f1_reproducibility_reporting_v2.tsv")

status={
 "schema_version":"part1_reporting_rules_v2",
 "status":"PASS_REPORTING_TABLES_COMPLETE",
 "completed_at":datetime.now().astimezone().isoformat(timespec="seconds"),
 "independent_unit":"patient",
 "minimum_units_for_inferential_reporting":10,
 "association_grid_rows":len(assoc),
 "full_statistics_rows":int(assoc["reporting_tier"].eq("FULL_STATISTICS_GE10").sum()),
 "direction_only_rows":int(assoc["reporting_tier"].eq("DIRECTION_ONLY_LT10").sum()),
 "not_tested_rows":int(assoc["reporting_tier"].eq("NOT_TESTED").sum()),
 "fdr_supported_rows_ge10":int(assoc["formal_result"].eq("FDR_SUPPORTED").sum()),
 "figures_generated":False
}
(OUT/"REPORTING_TABLES_STATUS_v2.json").write_text(json.dumps(status,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
print(json.dumps(status,ensure_ascii=False,indent=2))
