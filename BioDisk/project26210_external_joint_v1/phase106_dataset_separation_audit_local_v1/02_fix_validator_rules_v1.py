from pathlib import Path


path = Path(r"C:\Users\31937\project26210_phase106_dataset_separation_audit_v1\01_validate_phase106_dataset_separation_v1.py")
text = path.read_text(encoding="utf-8")

old_route = 'check("route_spatial_separate", all(any(row["primary_universe_ids"] == key for row in routes) for key in ["D04", "D05", "D06"]))'
new_route = '''check(
    "route_spatial_separate",
    all(any(key in row["primary_universe_ids"].split(";") for row in routes) for key in ["D04", "D05", "D06"])
    and all(
        row["cross_universe_status"] in {"SINGLE_UNIVERSE", "SINGLE_SPATIAL_UNIVERSE", "TWO_SPATIAL_COHORTS_KEEP_SEPARATE", "CROSS_SOURCE_ANALYSIS_EXPLICIT", "PREPARATION_STRATIFIED_SINGLE_UNIVERSE"}
        for row in routes
        if any(key in row["primary_universe_ids"].split(";") for key in ["D04", "D05", "D06"])
    ),
)'''
old_conflict = 'check("conflict_spatial_edges", any(row["issue_id"] == "N06" and "different edges" in row["correct_interpretation"] for row in conflicts))'
new_conflict = 'check("conflict_spatial_edges", any(row["issue_id"] == "N06" and row["classification"] == "DIFFERENT_EDGES_WITHIN_SAME_UNIVERSE" for row in conflicts))'

if text.count(old_route) != 1 or text.count(old_conflict) != 1:
    raise RuntimeError("validator source contract not matched exactly")
path.write_text(text.replace(old_route, new_route).replace(old_conflict, new_conflict), encoding="utf-8", newline="\n")
print("PHASE106_VALIDATOR_RULE_FIX_PASS")
