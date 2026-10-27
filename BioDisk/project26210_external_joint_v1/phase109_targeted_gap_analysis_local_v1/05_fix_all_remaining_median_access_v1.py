from pathlib import Path


path = Path(r"C:\Users\31937\project26210_phase109_targeted_gap_analysis_v1\02_run_phase109_targeted_analysis_v1.py")
text = path.read_text(encoding="utf-8")
count = text.count("row.median")
if count < 1:
    raise RuntimeError("no remaining row.median token found")
path.write_text(text.replace("row.median", 'row["median"]'), encoding="utf-8", newline="\n")
print(f"PHASE109_REMAINING_MEDIAN_FIX_PASS count={count}")
