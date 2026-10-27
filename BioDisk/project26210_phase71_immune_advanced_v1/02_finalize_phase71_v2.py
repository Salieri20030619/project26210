#!/usr/bin/env python3
import hashlib
import json
import os
from datetime import datetime
from pathlib import Path

import pandas as pd

ROOT = Path("/mnt/BioDisk/wangmingyuan/project26210_phase71_immune_advanced_v1")
RUNS = ROOT / "dataset_runs_v1"
PART = ROOT / "results_v2.part"
FINAL = ROOT / "results_v2"
DATASET_DIR = {"GSE207493":"GSE207493_v2","GSE264586":"GSE264586","GSE294109":"GSE294109"}
EXTERNAL = ("GSE264586","GSE294109")


def read(dataset, suffix):
    return pd.read_csv(RUNS / DATASET_DIR[dataset] / f"{dataset}_{suffix}", sep="\t")


def write(frame, path):
    temp = path.with_name(path.name + ".part")
    frame.to_csv(temp, sep="\t", index=False, na_rep="NA", lineterminator="\n")
    os.replace(temp, path)


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


if FINAL.exists() or PART.exists():
    raise SystemExit("Refusing existing final output")
PART.mkdir(parents=True)

summaries=[]; quartiles=[]; coverage=[]; statuses=[]
for dataset in DATASET_DIR:
    summaries.append(read(dataset,"relation_summary_v1.tsv"))
    quartiles.append(read(dataset,"quartile_summary_v1.tsv"))
    coverage.append(read(dataset,"program_coverage_v1.tsv"))
    statuses.append(read(dataset,"status_v1.tsv"))
summary=pd.concat(summaries,ignore_index=True)
quartile=pd.concat(quartiles,ignore_index=True)
coverage=pd.concat(coverage,ignore_index=True)
status=pd.concat(statuses,ignore_index=True)
if len(summary)!=39 or len(summary.query("tier=='PRIMARY'"))!=18:
    raise RuntimeError(f"Unexpected summary dimensions {summary.shape}")
if not status.status.eq("PASS_PHASE71_DATASET_V1").all():
    raise RuntimeError("Dataset status failure")

cross=[]
for relation_id, group in summary.groupby("relation_id",sort=True):
    first=group.iloc[0]
    by={r.dataset:r for r in group.itertuples(index=False)}
    supported=sum(by[d].decision=="COHORT_SUPPORTED" for d in EXTERNAL)
    associated=sum(by[d].decision=="ASSOCIATED_TWO_SIDED" for d in EXTERNAL)
    if first.expected_direction=="two_sided":
        decision="TWO_EXTERNAL_ASSOCIATED" if associated==2 else ("ONE_EXTERNAL_ASSOCIATED" if associated==1 else "NOT_SUPPORTED")
    else:
        decision="TWO_EXTERNAL_SUPPORTED" if supported==2 else ("ONE_EXTERNAL_SUPPORTED" if supported==1 else "NOT_SUPPORTED")
    row={"relation_id":relation_id,"tier":first.tier,"predictor":first.predictor,"outcome":first.outcome,"expected_direction":first.expected_direction,"external_decision":decision,"GSE207493_decision":by["GSE207493"].decision}
    for d in DATASET_DIR:
        row[f"{d}_n_patients"]=by[d].n_patients
        row[f"{d}_median_adjusted_rho"]=by[d].median_adjusted_rho
        row[f"{d}_q_value"]=by[d].q_value
        row[f"{d}_positive_fraction"]=by[d].positive_fraction
        row[f"{d}_decision"]=by[d].decision
    cross.append(row)
cross=pd.DataFrame(cross).sort_values(["tier","relation_id"]).reset_index(drop=True)
primary=cross.loc[cross.tier.eq("PRIMARY")].copy()
bridge_edges=primary.loc[primary.relation_id.isin(["P01","P02","P03"])]
bridge_decision="TRANSCRIPTIONAL_BRIDGE_SUPPORTED" if bridge_edges.external_decision.eq("TWO_EXTERNAL_SUPPORTED").all() else "TRANSCRIPTIONAL_BRIDGE_NOT_FULLY_SUPPORTED"
bridge=pd.DataFrame([{"bridge":"TYMP_to_MTA_response_to_SPP1_state","required_relations":"P01;P02;P03","relations_two_external_supported":int(bridge_edges.external_decision.eq("TWO_EXTERNAL_SUPPORTED").sum()),"required_count":3,"decision":bridge_decision,"boundary":"expression proxy; not measured MTA or causality"}])

write(summary.sort_values(["tier","relation_id","dataset"]),PART/"phase71_relation_summary_all_v1.tsv")
write(quartile.sort_values(["tier","relation_id","dataset"]),PART/"phase71_quartile_summary_all_v1.tsv")
write(coverage.sort_values(["program","dataset"]),PART/"phase71_program_coverage_all_v1.tsv")
write(status.sort_values("dataset"),PART/"phase71_dataset_status_all_v1.tsv")
write(cross,PART/"phase71_cross_cohort_decisions_v1.tsv")
write(bridge,PART/"phase71_primary_bridge_decision_v1.tsv")

def relation_line(row):
    return f"- {row.relation_id} {row.predictor} -> {row.outcome}: external={row.external_decision}; discovery={row.GSE207493_decision}; GSE264586 rho={row.GSE264586_median_adjusted_rho:.4f}, q={row.GSE264586_q_value:.4g}; GSE294109 rho={row.GSE294109_median_adjusted_rho:.4f}, q={row.GSE294109_q_value:.4g}."

lines=[
    "# Phase71 immune advanced analysis: MTA-response / SPP1 myeloid bridge", "",
    "## Scope", "",
    f"- Three serially analyzed cohorts; {int(status.selected_cells.sum()):,} frozen Myeloid/DC cells and {int(status.eligible_patients.sum())} patient-level units across cohort-specific analyses.",
    "- GSE207493 is discovery/context; GSE264586 and GSE294109 are the two external decision cohorts.",
    "- TYMP, the literature-derived MTA/A2B/STAT3 response proxy, and each outcome were evaluated within patient/timepoint after S/G2M residualization. Cells were not treated as replicates.", "",
    "## Primary results", "",
]
lines += [relation_line(r) for r in primary.itertuples(index=False)]
lines += ["", "## Bridge decision", "", f"- {bridge_decision}: {int(bridge.iloc[0].relations_two_external_supported)}/3 required external edges met the two-cohort rule.", "", "## Secondary results", ""]
lines += [relation_line(r) for r in cross.loc[cross.tier.eq("SECONDARY")].itertuples(index=False)]
lines += ["", "## Interpretation boundary", "", "Phase71 tests transcriptional convergence with a literature-derived MTA macrophage-response proxy. It does not measure 5'-MTA, establish TYMP as the biochemical cause of MTA accumulation, prove ADORA2B/STAT3 activation, demonstrate lineage conversion, or change the frozen response/PFS/OS grade.", ""]
(PART/"PHASE71_IMMUNE_ADVANCED_SUMMARY_v1.md").write_text("\n".join(lines),encoding="utf-8")
metadata={"status":"PASS_PHASE71_FINALIZATION_VALIDATION_PENDING","completed_at":datetime.now().astimezone().isoformat(timespec="seconds"),"datasets":list(DATASET_DIR),"selected_cells":int(status.selected_cells.sum()),"patient_units_sum":int(status.eligible_patients.sum()),"summary_rows":len(summary),"cross_rows":len(cross),"bridge_decision":bridge_decision}
(PART/"run_metadata_v1.json").write_text(json.dumps(metadata,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
manifest=[]
for path in sorted(PART.iterdir()):
    if path.is_file(): manifest.append({"file":path.name,"size_bytes":path.stat().st_size,"sha256":sha256(path)})
write(pd.DataFrame(manifest),PART/"FILE_MANIFEST_SHA256_v1.tsv")
(PART/"PHASE71_FINALIZATION_PASS_v1").write_text("PASS\n",encoding="ascii")
os.replace(PART,FINAL)
print(json.dumps(metadata,ensure_ascii=False))
