from pathlib import Path


path = Path(r"C:\Users\31937\project26210_phase109_targeted_gap_analysis_v1\02_run_phase109_targeted_analysis_v1.py")
text = path.read_text(encoding="utf-8")
text = text.replace(
    "def fixed_study_test(frame: pd.DataFrame, exposure: str, outcome: str, count_adjust: bool = False) -> dict:",
    "def fixed_study_test(frame: pd.DataFrame, exposure: str, outcome: str, count_adjust: bool = False, permutations: int = PERMUTATIONS) -> dict:",
)
old = '''    rng = np.random.default_rng(SEED)
    exceed = 0
    for _ in range(PERMUTATIONS):
        perm = yr.copy()
        for level in levels:
            idx = np.flatnonzero(studies == level)
            perm[idx] = rng.permutation(perm[idx])
        value = correlation(xr, perm)
        if np.isfinite(value) and abs(value) >= abs(observed) - 1e-15:
            exceed += 1
    p = (exceed + 1) / (PERMUTATIONS + 1)
'''
new = '''    p = float("nan")
    if permutations > 0:
        rng = np.random.default_rng(SEED)
        exceed = 0
        for _ in range(permutations):
            perm = yr.copy()
            for level in levels:
                idx = np.flatnonzero(studies == level)
                perm[idx] = rng.permutation(perm[idx])
            value = correlation(xr, perm)
            if np.isfinite(value) and abs(value) >= abs(observed) - 1e-15:
                exceed += 1
        p = (exceed + 1) / (permutations + 1)
'''
if old not in text:
    raise RuntimeError("permutation block not found")
text = text.replace(old, new)
text = text.replace(
    "loo = fixed_study_test(subset, exposure, outcome, count_adjust=False)",
    "loo = fixed_study_test(subset, exposure, outcome, count_adjust=False, permutations=0)",
)
path.write_text(text, encoding="utf-8", newline="\n")
print("PHASE109_LOO_OPTIMIZATION_PASS")
