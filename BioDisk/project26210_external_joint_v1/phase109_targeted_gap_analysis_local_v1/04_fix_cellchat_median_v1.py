from pathlib import Path


path = Path(r"C:\Users\31937\project26210_phase109_targeted_gap_analysis_v1\02_run_phase109_targeted_analysis_v1.py")
text = path.read_text(encoding="utf-8")
old = "row.n_patients >= 10 and row.median > 0 and row.positive_fraction >= 2 / 3 and row.bh_q_4 < 0.05 and all_loo"
new = "row.n_patients >= 10 and row[\"median\"] > 0 and row.positive_fraction >= 2 / 3 and row.bh_q_4 < 0.05 and all_loo"
if text.count(old) != 1:
    raise RuntimeError(f"expected failure expression once, observed {text.count(old)}")
path.write_text(text.replace(old, new), encoding="utf-8", newline="\n")
print("PHASE109_CELLCHAT_MEDIAN_FIX_PASS")
