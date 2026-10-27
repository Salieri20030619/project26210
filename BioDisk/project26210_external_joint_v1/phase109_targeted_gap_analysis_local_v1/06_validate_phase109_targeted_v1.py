from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, wilcoxon


STATES = ["SPP1-state Macrophage", "C1QC-state Macrophage", "Other Macrophage"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bh(p):
    p = np.asarray(p, float); order = np.argsort(p); ranked = p[order]
    adj = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(adj); out[order] = np.clip(adj, 0, 1); return out


def within_rank(values, studies):
    out = np.zeros(len(values), float)
    for study in sorted(set(studies)):
        idx = np.flatnonzero(studies == study); ranks = rankdata(values[idx]); sd = ranks.std(ddof=0)
        out[idx] = (ranks - ranks.mean()) / sd if sd > 0 else 0
    return out


def fixed_estimate(frame, exposure, outcome):
    studies = frame.study_v1.to_numpy(str); levels = sorted(set(studies))
    x = within_rank(frame[exposure].to_numpy(float), studies); y = within_rank(frame[outcome].to_numpy(float), studies)
    design = np.column_stack([np.ones(len(frame))] + [(studies == level).astype(float) for level in levels[1:]])
    x = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    y = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    return float(np.corrcoef(x, y)[0, 1])


parser = argparse.ArgumentParser()
parser.add_argument("--assignments", type=Path, required=True)
parser.add_argument("--axis-scores", type=Path, required=True)
parser.add_argument("--cellchat-interactions", type=Path, required=True)
parser.add_argument("--cellchat-qc", type=Path, required=True)
parser.add_argument("--spatial-decisions", type=Path, required=True)
parser.add_argument("--analysis", type=Path, required=True)
parser.add_argument("--failure-v1", type=Path, required=True)
parser.add_argument("--failure-v2", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)

checks = []
def check(name, passed, detail=""):
    checks.append({"check": name, "passed": bool(passed), "detail": str(detail)})

required = [
    "phase109_patient_summary_v1.tsv", "abundance_primary_results_v1.tsv",
    "abundance_count_adjusted_sensitivity_v1.tsv", "abundance_leave_one_study_out_v1.tsv",
    "abundance_study_specific_v1.tsv", "abundance_exposure_sensitivity_v1.tsv",
    "cellchat_patient_fine_state_specificity_v1.tsv", "cellchat_fine_state_specificity_summary_v1.tsv",
    "cellchat_leave_one_patient_out_v1.tsv", "spatial_frozen_edge_synthesis_v1.tsv",
    "PHASE109_TARGETED_ANALYSIS_SUMMARY_v1.md", "PHASE109_ANALYSIS_STATUS_v1.json",
]
for name in required:
    check(f"required_{name}", (args.analysis / name).is_file())
check("failure_v1_preserved", args.failure_v1.is_file())
check("failure_v2_preserved", args.failure_v2.is_file())

status = json.loads((args.analysis / "PHASE109_ANALYSIS_STATUS_v1.json").read_text(encoding="utf-8"))
check("formal_status", status.get("status") == "PASS_ANALYSIS_COMPLETE")
check("formal_counts", status.get("eligible_abundance_patients") == 20 and status.get("cellchat_patients") == 14)
check("formal_support_counts", status.get("primary_abundance_supported") == 0 and status.get("fine_state_cellchat_supported") == 2)
check("no_direct_mta_claim", status.get("direct_mta_tested") is False and status.get("causal_mta_to_spp1_tested") is False)

assign = pd.read_csv(args.assignments, sep="\t")
axis = pd.read_csv(args.axis_scores, sep="\t")
patient_reported = pd.read_csv(args.analysis / "phase109_patient_summary_v1.tsv", sep="\t")
counts = assign.loc[assign.state_v2.isin(STATES)].groupby(["study_v1", "patient", "state_v2"]).size().unstack(fill_value=0).reset_index()
for state in STATES:
    if state not in counts: counts[state] = 0
counts["macrophage_total"] = counts[STATES].sum(axis=1)
counts["spp1_fraction"] = counts[STATES[0]] / counts.macrophage_total
counts["spp1_vs_c1qc_logratio"] = np.log((counts[STATES[0]] + .5) / (counts[STATES[1]] + .5))
axis_patient = axis.groupby(["study", "curated_patient_id"])[["TYMP_axis_primary", "z_TYMP_RNA", "TYMP_axis_broad_sensitivity", "F1"]].median().reset_index().rename(columns={"study":"study_v1","curated_patient_id":"patient"})
recomputed = counts.merge(axis_patient, on=["study_v1","patient"])
reported_sorted = patient_reported.sort_values(["study_v1","patient"]).reset_index(drop=True)
recomputed_sorted = recomputed.sort_values(["study_v1","patient"]).reset_index(drop=True)
check("patient_keys_recomputed", reported_sorted[["study_v1","patient"]].equals(recomputed_sorted[["study_v1","patient"]]))
for column in [*STATES, "macrophage_total", "spp1_fraction", "spp1_vs_c1qc_logratio", "TYMP_axis_primary", "z_TYMP_RNA", "TYMP_axis_broad_sensitivity", "F1"]:
    check(f"patient_column_{column}", np.allclose(reported_sorted[column], recomputed_sorted[column], equal_nan=True), np.nanmax(np.abs(reported_sorted[column]-recomputed_sorted[column])))

eligible = reported_sorted.loc[reported_sorted.eligible_abundance_primary.astype(bool)]
primary = pd.read_csv(args.analysis / "abundance_primary_results_v1.tsv", sep="\t")
check("primary_rows_4", len(primary) == 4)
for _, row in primary.iterrows():
    estimate = fixed_estimate(eligible, row.exposure, row.outcome)
    check(f"primary_estimate_{row.exposure}_{row.outcome}", abs(estimate - row.estimate) < 1e-12, estimate)
check("primary_zero_supported", int(primary.support.astype(bool).sum()) == 0)
check("primary_all_q_nonsignificant", bool((primary.bh_q_4 >= 0.05).all()))

interactions = pd.read_csv(args.cellchat_interactions, sep="\t")
qc = pd.read_csv(args.cellchat_qc, sep="\t")
pair = interactions.groupby(["study","patient","source","target"], as_index=False).prob.sum()
def val(z,s,t):
    x=z.loc[(z.source==s)&(z.target==t),"prob"]; return float(x.iloc[0]) if len(x) else 0.0
chat_rows=[]
for _,q in qc.iterrows():
    z=pair.loc[pair.patient==q.patient]
    h={s:val(z,"TYMP-high Malignant",s) for s in STATES};l={s:val(z,"TYMP-low Malignant",s) for s in STATES}
    th={s:val(z,s,"TYMP-high Malignant") for s in STATES};tl={s:val(z,s,"TYMP-low Malignant") for s in STATES}
    chat_rows.append({"study":q.study,"patient":q.patient,
      "malignant_to_spp1_specificity":(h[STATES[0]]-l[STATES[0]])-np.mean([h[s]-l[s] for s in STATES[1:]]),
      "spp1_to_malignant_specificity":(th[STATES[0]]-tl[STATES[0]])-np.mean([th[s]-tl[s] for s in STATES[1:]]),
      "spp1_to_endothelial_specificity":val(z,STATES[0],"Endothelial")-np.mean([val(z,s,"Endothelial") for s in STATES[1:]]),
      "endothelial_to_spp1_specificity":val(z,"Endothelial",STATES[0])-np.mean([val(z,"Endothelial",s) for s in STATES[1:]])})
chat_recomputed=pd.DataFrame(chat_rows).sort_values("patient").reset_index(drop=True)
chat_reported=pd.read_csv(args.analysis/"cellchat_patient_fine_state_specificity_v1.tsv",sep="\t").sort_values("patient").reset_index(drop=True)
check("chat_keys",chat_reported[["study","patient"]].equals(chat_recomputed[["study","patient"]]))
for column in chat_recomputed.columns[2:]:check(f"chat_patient_{column}",np.allclose(chat_reported[column],chat_recomputed[column]),np.max(np.abs(chat_reported[column]-chat_recomputed[column])))
summary=pd.read_csv(args.analysis/"cellchat_fine_state_specificity_summary_v1.tsv",sep="\t")
pvals=[]
for _,row in summary.iterrows():
    values=chat_recomputed[row.endpoint].to_numpy(float);test=wilcoxon(values,zero_method="wilcox",alternative="two-sided",method="approx")
    check(f"chat_median_{row.endpoint}",abs(np.median(values)-row["median"])<1e-12)
    check(f"chat_p_{row.endpoint}",abs(test.pvalue-row.wilcoxon_p)<1e-12);pvals.append(test.pvalue)
check("chat_bh",np.allclose(bh(pvals),summary.bh_q_4))
supported=set(summary.loc[summary.support.astype(bool),"endpoint"])
check("chat_supported_exact",supported=={"spp1_to_endothelial_specificity","endothelial_to_spp1_specificity"},supported)

spatial=pd.read_csv(args.spatial_decisions,sep="\t");expected={"B01_TYMP_to_SPP1":"NOT_SUPPORTED","B08_SPP1_to_Tdysfunction":"ONE_PREPARATION_SUPPORTED","B09_SPP1_to_Vascular":"TWO_PREPARATION_SUPPORTED"}
for relation,decision in expected.items():check(f"spatial_{relation}",spatial.loc[spatial.relation_id==relation,"decision"].iloc[0]==decision)

failed=[x for x in checks if not x["passed"]]
payload={"status":"PASS_PHASE109_INDEPENDENT_VALIDATION" if not failed else "FAIL_PHASE109_INDEPENDENT_VALIDATION","checks_passed":len(checks)-len(failed),"checks_total":len(checks),"failed_checks":failed,"formal_output_hashes":{name:sha256(args.analysis/name) for name in required}}
(args.output/"PHASE109_INDEPENDENT_VALIDATION_STATUS_v1.json").write_text(json.dumps(payload,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
pd.DataFrame(checks).to_csv(args.output/"phase109_independent_validation_checks_v1.tsv",sep="\t",index=False)
print(json.dumps(payload,ensure_ascii=False,indent=2));raise SystemExit(0 if not failed else 1)
