from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("phase90_v8_for_v10", SCRIPT_DIR / "12_build_main_figure_v8.py")
v8 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(v8)
base = v8.base
original_bridge = v8.draw_replication_bridge_v8


def draw_replication_bridge_v10(ax, selected) -> None:
    original_bridge(ax, selected)
    for text in list(ax.texts):
        if text.get_text().startswith("Fill = rho"):
            text.remove()
    ax.text(
        1.58, -0.48,
        "Fill: rho (teal − / magenta +)  •  area: −log10(q)  •  coloured line: same direction",
        ha="right", va="center", fontsize=4.45, color=base.COLORS["muted"],
    )


def draw_panel_d_v10(fig, spec, retained) -> None:
    outer = fig.add_subplot(spec)
    base.clean_axis(outer)
    base.panel_heading(outer, "D", "Cross-cohort replication of TYMP-associated immune programs", None)
    gs = spec.subgridspec(1, 2, width_ratios=[1.78, 0.82], wspace=0.12)
    selected = base.select_phase67_rows(retained)
    ax = fig.add_subplot(gs[0, 0])
    draw_replication_bridge_v10(ax, selected)
    ax_net = fig.add_subplot(gs[0, 1])
    base.draw_program_constellation(ax_net, selected)
    ax_net.text(-1.18, 1.19, f"{sum(selected.endpoint == 'myeloid_dc_tymp')} two-cohort Myeloid/DC signatures", fontsize=5.3, color=base.COLORS["tymp"], fontweight="bold")
    ax_net.text(-1.18, 1.08, f"{sum(selected.endpoint == 'treg_tymp')} single-external Treg signatures", fontsize=5.3, color=base.COLORS["c1qc"], fontweight="bold")


def rehash_outputs_v10(output_dir: Path) -> None:
    summary_path = output_dir / "BUILD_SUMMARY_v1.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    renamed = {}
    for ext in ["png", "pdf", "svg"]:
        old = output_dir / f"PROJECT26210_IMMUNE_MAIN_FIGURE_v1.{ext}"
        new = output_dir / f"PROJECT26210_IMMUNE_MAIN_FIGURE_v10.{ext}"
        old.replace(new)
        renamed[ext] = {"path": str(new), "bytes": new.stat().st_size, "sha256": base.sha256(new)}
    summary["status"] = "GENERATED_PENDING_FINAL_QA_V10"
    summary["outputs"] = renamed
    (output_dir / "BUILD_SUMMARY_v10.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    summary_path.unlink()


v8.draw_replication_bridge_v8 = draw_replication_bridge_v10
base.draw_panel_d = draw_panel_d_v10


if __name__ == "__main__":
    base.main()
    rehash_outputs_v10(Path(sys.argv[3]))
