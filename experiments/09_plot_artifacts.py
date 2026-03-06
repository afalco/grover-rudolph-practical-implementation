# experiments/09_plot_artifacts.py
#
# Plot artifacts from gr_log_*.jsonl produced by experiments/08_full_gr_log_to_csv.py
#
# - No pandas dependency (uses stdlib + numpy + matplotlib)
# - Produces clean plots (not empty):
#     * Boxplots by stage for TV/L2/Fidelity (SIM in blue)
#     * Optional histograms per stage
#     * FULL-stage summary table (CSV + LaTeX)
#
# Usage:
#   python -m experiments.09_plot_artifacts --input artifacts/gr_log_YYYYMMDD-HHMMSS.jsonl
#   python -m experiments.09_plot_artifacts --input ... --include-avg
#   python -m experiments.09_plot_artifacts --input ... --hist
#
# Output (default):
#   artifacts/gr_plot_<tag>_tv_box.png/.pdf
#   artifacts/gr_plot_<tag>_l2_box.png/.pdf
#   artifacts/gr_plot_<tag>_fid_box.png/.pdf
#   artifacts/gr_summary_<tag>.csv
#   artifacts/gr_summary_<tag>.tex

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt

# -------------------------- helpers --------------------------

STAGES_ORDER = ["L0", "L01", "FULL"]

def _norm_stage(stage: str) -> str:
    s = str(stage).strip().upper()
    if s == "FULL" or s == "FULL ":
        return "FULL"
    if s in ("L0", "L01"):
        return s
    # accept "full"
    if s.lower() == "full":
        return "FULL"
    return s

def _safe_float(x: Any, default: float = float("nan")) -> float:
    try:
        return float(x)
    except Exception:
        return default

def _ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def _infer_tag(input_path: str) -> str:
    # Prefer timestamp inside filename gr_log_YYYYMMDD-HHMMSS.jsonl
    m = re.search(r"gr_log_(\d{8}-\d{6})", os.path.basename(input_path))
    return m.group(1) if m else "unknown"

def _write_summary_csv(path: str, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    _ensure_dir(os.path.dirname(path) or ".")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)

def _latex_escape(s: str) -> str:
    return (
        s.replace("\\", "\\textbackslash{}")
         .replace("_", "\\_")
         .replace("%", "\\%")
         .replace("&", "\\&")
         .replace("#", "\\#")
    )

def _write_summary_tex(path: str, caption: str, label: str, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    _ensure_dir(os.path.dirname(path) or ".")
    lines: List[str] = []
    lines.append("\\begin{table}[t]")
    lines.append("\\centering")
    lines.append("\\small")
    lines.append("\\begin{tabular}{lccc}")
    lines.append("\\hline")
    lines.append("Comparison & TV & L2 & Fidelity \\\\")
    lines.append("\\hline")
    for r in rows:
        comp = _latex_escape(str(r["comparison"]))
        tv = f"{float(r['tv']):.6g}"
        l2 = f"{float(r['l2']):.6g}"
        fid = f"{float(r['fidelity']):.6g}"
        lines.append(f"{comp} & {tv} & {l2} & {fid} \\\\")
    lines.append("\\hline")
    lines.append("\\end{tabular}")
    lines.append(f"\\caption{{{_latex_escape(caption)}}}")
    lines.append(f"\\label{{{_latex_escape(label)}}}")
    lines.append("\\end{table}")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

# -------------------------- parsing --------------------------

@dataclass
class Rec:
    backend: str          # SIM or NMR
    stage: str            # L0/L01/FULL
    ladder: str
    repeats: int          # 1 for per-run, >1 for avg rows
    readout_mitigation: bool

    tv_vs_target: float
    l2_vs_target: float
    fid_vs_target: float

    tv_vs_sim: float
    l2_vs_sim: float
    fid_vs_sim: float

    note: str = ""

def load_jsonl(path: str) -> List[Rec]:
    recs: List[Rec] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            obj = json.loads(line)
            r = obj.get("record", {})
            recs.append(
                Rec(
                    backend=str(r.get("backend", "")).upper(),
                    stage=_norm_stage(r.get("stage", "")),
                    ladder=str(r.get("ladder", "")).upper(),
                    repeats=int(r.get("repeats", 1)),
                    readout_mitigation=bool(r.get("readout_mitigation", False)),
                    tv_vs_target=_safe_float(r.get("tv_vs_target")),
                    l2_vs_target=_safe_float(r.get("l2_vs_target")),
                    fid_vs_target=_safe_float(r.get("fid_vs_target")),
                    tv_vs_sim=_safe_float(r.get("tv_vs_sim")),
                    l2_vs_sim=_safe_float(r.get("l2_vs_sim")),
                    fid_vs_sim=_safe_float(r.get("fid_vs_sim")),
                    note=str(obj.get("note", "")),
                )
            )
    return recs

# -------------------------- plotting --------------------------

def boxplot_by_stage(
    recs: List[Rec],
    metric_key: str,
    out_png: str,
    out_pdf: str,
    include_avg: bool,
    title: str,
) -> None:
    """
    metric_key in:
      tv_vs_target, l2_vs_target, fid_vs_target,
      tv_vs_sim,    l2_vs_sim,    fid_vs_sim
    """
    # Separate SIM and NMR; keep per-run only by default
    def select(rs: List[Rec], backend: str, stage: str) -> List[float]:
        vals = []
        for r in rs:
            if r.backend != backend:
                continue
            if r.stage != stage:
                continue
            if (not include_avg) and r.repeats != 1:
                continue
            vals.append(getattr(r, metric_key))
        return [v for v in vals if not (math.isnan(v) or math.isinf(v))]

    data_sim = [select(recs, "SIM", st) for st in STAGES_ORDER]
    data_nmr = [select(recs, "NMR", st) for st in STAGES_ORDER]

    # If SIM has only one point per stage, boxplot looks weird; we draw it as a scatter.
    fig = plt.figure()
    ax = plt.gca()

    positions_sim = [1, 3, 5]
    positions_nmr = [2, 4, 6]

    # NMR boxplots
    bp = ax.boxplot(
        data_nmr,
        positions=positions_nmr,
        widths=0.7,
        patch_artist=True,
        showfliers=True,
    )
    # keep default colors; just set light facecolor for readability
    for patch in bp["boxes"]:
        patch.set_alpha(0.35)

    # SIM points (blue)
    for x, vals in zip(positions_sim, data_sim):
        if len(vals) == 0:
            continue
        # plot as points with small jitter
        jitter = np.linspace(-0.08, 0.08, num=len(vals)) if len(vals) > 1 else [0.0]
        ax.scatter([x + j for j in jitter], vals, marker="o", s=35)

    ax.set_xticks([1.5, 3.5, 5.5])
    ax.set_xticklabels(STAGES_ORDER)
    ax.set_title(title)
    ax.set_xlabel("Stage")
    ax.set_ylabel(metric_key)

    # Legend proxy
    ax.scatter([], [], marker="o", label="SIM (ideal)", color="C0")
    ax.plot([], [], label="NMR (box)", color="C1")
    ax.legend(loc="best")

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_pdf)
    plt.close(fig)

def hist_by_stage(
    recs: List[Rec],
    metric_key: str,
    out_png: str,
    out_pdf: str,
    include_avg: bool,
    title: str,
) -> None:
    def select(backend: str, stage: str) -> List[float]:
        vals = []
        for r in recs:
            if r.backend != backend:
                continue
            if r.stage != stage:
                continue
            if (not include_avg) and r.repeats != 1:
                continue
            vals.append(getattr(r, metric_key))
        vals = [v for v in vals if not (math.isnan(v) or math.isinf(v))]
        return vals

    fig = plt.figure()
    ax = plt.gca()

    # concatenate per-stage NMR histograms; SIM shown as vertical line(s)
    for st in STAGES_ORDER:
        nmr_vals = select("NMR", st)
        if nmr_vals:
            ax.hist(nmr_vals, bins=15, alpha=0.35, label=f"NMR {st}")
        sim_vals = select("SIM", st)
        for v in sim_vals:
            ax.axvline(v, linestyle="--", linewidth=1.2)

    ax.set_title(title)
    ax.set_xlabel(metric_key)
    ax.set_ylabel("count")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_pdf)
    plt.close(fig)

# -------------------------- summary table --------------------------

def build_full_summary(recs: List[Rec]) -> List[Dict[str, Any]]:
    """
    Produce a Viana-style FULL summary table.
    Prefer NMR avg rows (repeats>1) if present.
    """
    rows: List[Dict[str, Any]] = []

    def pick(backend: str, want_avg: bool) -> Optional[Rec]:
        cand = [r for r in recs if r.backend == backend and r.stage == "FULL"]
        if want_avg:
            cand = [r for r in cand if r.repeats != 1]
            # prefer explicit avg note
            cand2 = [r for r in cand if "AVERAGE" in (r.note or "").upper()]
            if cand2:
                cand = cand2
        else:
            cand = [r for r in cand if r.repeats == 1]
        return cand[-1] if cand else None

    sim = pick("SIM", want_avg=True) or pick("SIM", want_avg=False)
    nmr_avg = pick("NMR", want_avg=True)

    if sim is not None:
        rows.append({
            "comparison": "Target vs SIM ideal (FULL)",
            "tv": float(sim.tv_vs_target),
            "l2": float(sim.l2_vs_target),
            "fidelity": float(sim.fid_vs_target),
        })
    if nmr_avg is not None:
        rows.append({
            "comparison": "Target vs NMR raw avg (FULL)",
            "tv": float(nmr_avg.tv_vs_target),
            "l2": float(nmr_avg.l2_vs_target),
            "fidelity": float(nmr_avg.fid_vs_target),
        })
        # If mitigation was enabled in that run, its tv_vs_target is already “final”.
        # If you log separate mitigated/mixed rows, you can add them here similarly.

    return rows

# -------------------------- main --------------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Path to artifacts/gr_log_*.jsonl")
    ap.add_argument("--outdir", default="artifacts", help="Output directory (default: artifacts)")
    ap.add_argument("--include-avg", action="store_true", help="Include repeats>1 rows in plots")
    ap.add_argument("--hist", action="store_true", help="Also output histograms per metric")
    args = ap.parse_args()

    tag = _infer_tag(args.input)
    _ensure_dir(args.outdir)

    recs = load_jsonl(args.input)
    if not recs:
        raise RuntimeError("No records found in JSONL.")

    # Boxplots: compare vs TARGET (the ones you care about most)
    for metric_key, short in [
        ("tv_vs_target", "tv"),
        ("l2_vs_target", "l2"),
        ("fid_vs_target", "fid"),
    ]:
        out_png = os.path.join(args.outdir, f"gr_plot_{tag}_{short}_box.png")
        out_pdf = os.path.join(args.outdir, f"gr_plot_{tag}_{short}_box.pdf")
        boxplot_by_stage(
            recs=recs,
            metric_key=metric_key,
            out_png=out_png,
            out_pdf=out_pdf,
            include_avg=args.include_avg,
            title=f"{metric_key} by stage (SIM points in blue; NMR boxplots)",
        )
        print("Wrote:", out_png)
        print("Wrote:", out_pdf)

        if args.hist:
            out_png_h = os.path.join(args.outdir, f"gr_plot_{tag}_{short}_hist.png")
            out_pdf_h = os.path.join(args.outdir, f"gr_plot_{tag}_{short}_hist.pdf")
            hist_by_stage(
                recs=recs,
                metric_key=metric_key,
                out_png=out_png_h,
                out_pdf=out_pdf_h,
                include_avg=args.include_avg,
                title=f"{metric_key} histograms (dashed lines = SIM)",
            )
            print("Wrote:", out_png_h)
            print("Wrote:", out_pdf_h)

    # FULL summary table (CSV + LaTeX)
    summary_rows = build_full_summary(recs)
    if summary_rows:
        summary_csv = os.path.join(args.outdir, f"gr_summary_{tag}.csv")
        summary_tex = os.path.join(args.outdir, f"gr_summary_{tag}.tex")
        _write_summary_csv(summary_csv, summary_rows)
        _write_summary_tex(
            summary_tex,
            caption="Grover–Rudolph (FULL stage) summary metrics.",
            label="tab:gr_full_summary",
            rows=summary_rows,
        )
        print("Wrote:", summary_csv)
        print("Wrote:", summary_tex)
    else:
        print("[warn] No FULL-stage rows found; no summary written.")

if __name__ == "__main__":
    main()
