from __future__ import annotations

from pathlib import Path


source_path = Path(__file__).resolve().parent / "07_validate_main_figure_v1.py"
source = source_path.read_text(encoding="utf-8")
old_phase87 = '''    record("phase87_source_validation", f"{sum(phase87_checks.status == 'PASS')}/{len(phase87_checks)}", phase87_checks["status"].eq("PASS").all())'''
new_phase87 = '''    phase87_pass = phase87_checks["pass"].astype(str).str.lower().eq("true")
    record("phase87_source_validation", f"{int(phase87_pass.sum())}/{len(phase87_checks)}", phase87_pass.all())'''
if old_phase87 not in source:
    raise SystemExit("Expected Phase87 validator line not found")
source = source.replace(old_phase87, new_phase87)
source = source.replace("PROJECT26210_IMMUNE_MAIN_FIGURE_v6", "PROJECT26210_IMMUNE_MAIN_FIGURE_v11")
source = source.replace("BUILD_SUMMARY_v6.json", "BUILD_SUMMARY_v11.json")
namespace = {"__name__": "__main__", "__file__": str(source_path)}
exec(compile(source, str(source_path), "exec"), namespace)
