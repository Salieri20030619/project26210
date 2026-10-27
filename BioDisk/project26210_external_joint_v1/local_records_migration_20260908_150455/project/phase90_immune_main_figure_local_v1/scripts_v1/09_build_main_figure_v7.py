from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("phase90_v2_for_v7", SCRIPT_DIR / "02_build_main_figure_v2.py")
v2 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(v2)
base = v2.base
original_draw_replication_bridge = v2.draw_replication_bridge
original_figure_text = base.mpl.figure.Figure.text


def panel_heading_v7(ax, letter: str, title: str, subtitle: str | None = None) -> None:
    ax.text(-0.015, 1.045, letter, transform=ax.transAxes, ha="right", va="bottom", fontsize=13.2, fontweight="bold", color=base.COLORS["ink"])
    ax.text(0.0, 1.045, title, transform=ax.transAxes, ha="left", va="bottom", fontsize=10.2, fontweight="bold", color=base.COLORS["ink"])
    if subtitle and letter in {"A", "B"}:
        ax.text(1.0, 1.010, subtitle, transform=ax.transAxes, ha="right", va="bottom", fontsize=5.8, color=base.COLORS["muted"])


def draw_replication_bridge_v7(ax, selected) -> None:
    original_draw_replication_bridge(ax, selected)
    for child in list(ax.child_axes):
        child.remove()
    ax.text(
        0.98, 0.012,
        "rho colour: teal −0.4  •  white 0  •  magenta +0.4   |   bubble area = −log10(q)",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=4.6, color=base.COLORS["muted"],
    )


def figure_text_v7(self, x, y, s, *args, **kwargs):
    if isinstance(s, str) and s.startswith("External immune evidence after cohort expansion"):
        s = ""
    return original_figure_text(self, x, y, s, *args, **kwargs)


def rehash_outputs_v7(output_dir: Path) -> None:
    summary_path = output_dir / "BUILD_SUMMARY_v1.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    renamed = {}
    for ext in ["png", "pdf", "svg"]:
        old = output_dir / f"PROJECT26210_IMMUNE_MAIN_FIGURE_v1.{ext}"
        new = output_dir / f"PROJECT26210_IMMUNE_MAIN_FIGURE_v7.{ext}"
        old.replace(new)
        renamed[ext] = {"path": str(new), "bytes": new.stat().st_size, "sha256": base.sha256(new)}
    summary["status"] = "GENERATED_PENDING_FINAL_QA_V7"
    summary["outputs"] = renamed
    (output_dir / "BUILD_SUMMARY_v7.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    summary_path.unlink()


base.panel_heading = panel_heading_v7
v2.draw_replication_bridge = draw_replication_bridge_v7
base.mpl.figure.Figure.text = figure_text_v7


if __name__ == "__main__":
    base.main()
    rehash_outputs_v7(Path(sys.argv[3]))
