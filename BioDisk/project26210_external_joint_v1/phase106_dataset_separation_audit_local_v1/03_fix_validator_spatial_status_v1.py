from pathlib import Path


path = Path(r"C:\Users\31937\project26210_phase106_dataset_separation_audit_v1\01_validate_phase106_dataset_separation_v1.py")
text = path.read_text(encoding="utf-8")
old = '{"SINGLE_UNIVERSE", "SINGLE_SPATIAL_UNIVERSE", "TWO_SPATIAL_COHORTS_KEEP_SEPARATE", "CROSS_SOURCE_ANALYSIS_EXPLICIT", "PREPARATION_STRATIFIED_SINGLE_UNIVERSE"}'
new = '{"SINGLE_UNIVERSE", "SINGLE_SPATIAL_UNIVERSE", "TWO_SPATIAL_COHORTS_KEEP_SEPARATE", "CROSS_SOURCE_ANALYSIS_EXPLICIT", "PREPARATION_STRATIFIED_SINGLE_UNIVERSE", "MULTI_OUTPUT_KEEP_SEPARATE"}'
if text.count(old) != 1:
    raise RuntimeError(f"expected spatial-status set once, observed {text.count(old)}")
path.write_text(text.replace(old, new), encoding="utf-8", newline="\n")
print("PHASE106_VALIDATOR_SPATIAL_STATUS_FIX_PASS")
