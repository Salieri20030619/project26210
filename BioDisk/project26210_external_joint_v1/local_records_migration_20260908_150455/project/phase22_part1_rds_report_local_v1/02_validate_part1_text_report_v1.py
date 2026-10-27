#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json
from datetime import datetime
from pathlib import Path
import pandas as pd

HERE=Path(__file__).resolve().parent
DATA=HERE/"extracted_data_v1"
REPORT=HERE/"PART1_EXTERNAL_SINGLE_CELL_REPORT_TEXT_DRAFT_v1.md"
SOURCE=HERE.parent/"phase14_patient_program_validation_local_v1"/"patient_analysis_final_cnv_v6"
OUT=HERE/"validation_v1"
if OUT.exists():raise FileExistsError(OUT)
OUT.mkdir()
checks=[]
def ck(name,value,detail):checks.append({"check":name,"status":"PASS" if bool(value) else "FAIL","detail":detail})
def sha(path):
 d=hashlib.sha256()
 with Path(path).open("rb") as h:
  for block in iter(lambda:h.read(8*1024*1024),b""):d.update(block)
 return d.hexdigest()

status=pd.read_csv(DATA/"PART1_RDS_REPORT_EXTRACTION_STATUS_v1.tsv",sep="\t")
ck("rds_extraction_status",len(status)==1 and status.loc[0,"status"]=="PASS_RDS_EXTRACTION" and int(status.loc[0,"datasets"])==5 and int(status.loc[0,"rds_rows"])==505639,"five RDS; 505639 rows")
audit=pd.read_csv(DATA/"rds_object_audit_v1.tsv",sep="\t")
ck("rds_identity",len(audit)==5 and audit["dataset"].nunique()==5 and audit["rows"].sum()==505639 and audit["duplicate_cell_ids"].sum()==0 and audit["status"].eq("PASS").all(),"5 unique objects; no duplicate cell IDs")
overview=pd.read_csv(DATA/"dataset_overview_from_rds_v1.tsv",sep="\t")
ck("patient_sample_coverage",overview["patients"].sum()==63 and overview["samples"].sum()==95 and overview["missing_patient_key_cells"].sum()==0 and overview["missing_sample_id_cells"].sum()==0,"63 patients; 95 samples; zero missing keys")
comp=pd.read_csv(DATA/"analysis_compartment_counts_from_rds_v1.tsv",sep="\t")
ck("compartment_conservation",comp["n_cells"].sum()==505639 and set(comp["analysis_compartment"])=={"t_nk","myeloid","normal_epithelial","endothelial","malignant_provisional","fibroblast"},"six target compartments; row conservation")
cnv=pd.read_csv(DATA/"final_cnv_dataset_summary_v1.tsv",sep="\t")
ck("final_cnv",len(cnv)==5 and cnv["compact_cells"].sum()==505639 and cnv["tumor_like_cells"].sum()==40302 and cnv["final_malignant_cells"].sum()==16645 and cnv["status"].eq("PASS").all(),"40302 candidates; 16645 final malignant")
relations=pd.read_csv(DATA/"relation_result_summary_v1.tsv",sep="\t")
ck("relation_grid",len(relations)==12 and relations["grid_rows"].sum()==120 and relations["tested_rows"].sum()==66 and relations["nominal_p_lt_0_05"].sum()==5 and relations["fdr_lt_0_05"].sum()==0,"120 grid; 66 tested; 5 nominal; zero FDR")
validation=pd.read_csv(SOURCE/"INDEPENDENT_VALIDATION_STATUS_v6.tsv",sep="\t")
ck("source_v6_validation",len(validation)==19 and validation.iloc[-1]["check"]=="OVERALL" and validation.iloc[-1]["result"]=="PASS","Phase14 v6 19/19 PASS")
text=REPORT.read_text(encoding="utf-8")
required=["505,639","40,302","16,645","120行","66行","没有任何关系通过","GSE207493","GSE294109","不依赖此前口头总结"]
ck("report_required_content",all(x in text for x in required),"all key counts and boundaries present")
ck("no_figures","![" not in text and not any(ext in text.lower() for ext in [".png",".jpg",".jpeg",".svg",".pdf"]),"no image links or figure files")
ck("draft_boundary","供内部审阅" in text and "尚未配置图件" in text and "不能表述为外部验证阳性" in text,"review and claim boundaries present")

frame=pd.DataFrame(checks);frame.to_csv(OUT/"PART1_TEXT_REPORT_VALIDATION_STATUS_v1.tsv",sep="\t",index=False)
passed=frame["status"].eq("PASS").all()
summary={"schema_version":"part1_text_report_validation_v1","status":"PASS_TECHNICAL_VALIDATION" if passed else "FAIL","validated_at":datetime.now().astimezone().isoformat(timespec="seconds"),"checks_passed":int(frame["status"].eq("PASS").sum()),"checks_total":len(frame),"report_sha256":sha(REPORT),"figures_generated":False,"content_review_status":"AWAITING_USER_REVIEW"}
(OUT/"PART1_TEXT_REPORT_VALIDATION_SUMMARY_v1.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
files=[REPORT,*sorted(DATA.glob("*")),Path(__file__),HERE/"01_extract_part1_rds_report_data_v1.R"]
manifest=pd.DataFrame([{"file":str(f.relative_to(HERE)).replace("\\","/"),"bytes":f.stat().st_size,"sha256":sha(f)} for f in files])
manifest.to_csv(OUT/"OUTPUT_MANIFEST_SHA256_v1.tsv",sep="\t",index=False)
print(json.dumps(summary,ensure_ascii=False,indent=2))
raise SystemExit(0 if passed else 1)
