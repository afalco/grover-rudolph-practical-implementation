# experiments/09_plot_artifacts.py
#
# Read artifacts/gr_log_*.jsonl produced by experiments/08_full_gr_log_to_csv.py
# and generate:
#   - artifacts/gr_metrics_<tag>.csv       (flat per-record metrics for stats)
#   - artifacts/gr_summary_<tag>.csv/.tex  ("Viana-style" summary table; FULL stage)
#   - artifacts/gr_boxplot_*_<tag>.png/.pdf  (TV/L2/Fidelity boxplots by stage)
#   - artifacts/gr_hist_*_<tag>.png/.pdf     (histograms for FULL stage metrics)
#
# Notes:
# - JSONL lines are expected to have at least: {"record": {...}, ...}
# - We use matplotlib (no seaborn).
# - SIM boxplots are colored BLUE.

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def _read_jsonl(path: str) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rows.append(json.loads(line))
    return rows


def _ensure_dir(d: str) -> None:
    os.makedirs(d, exist_ok=True)


def _tag_from_input(input_path: str) -> str:
    p = Path(input_path)
    stem = p.stem  # gr_log_....
    if stem.startswith("gr_log_"):
        return stem.replace("gr_log_", "", 1)
    return stem


def _latex_escape(s: str) -> str:
    return (
        s.replace("\\", "\\textbackslash{}")
         .replace("_", "\\_")
         .replace("%", "\\%")
         .replace("&", "\\&")
         .replace("#", "\\#")
    )


def write_summary_csv(path: str, rows: List[dict]) -> None:
    if not rows:
        return
    _ensure_dir(str(Path(path).parent))
    pd.DataFrame(rows).to_csv(path, index=False)


def write_summary_tex(path: str, caption: str, label: str, rows: List[dict]) -> None:
    if not rows:
        return
    _ensure_dir(str(Path(path).parent))
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
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def _normalize_stage(s: str) -> str:
    s2 = str(s).strip()
    if s2.lower() == "full":
        return "FULL"
    if s2.upper() in ("L0", "L01"):
        return s2.upper()
    return s2.upper()


def _extract_records(jsonl_rows: List[Dict[str, Any]]) -> pd.DataFrame:
    flat: List[Dict[str, Any]] = []
    for obj in jsonl_rows:
        rec = obj.get("record", None)
        if not isinstance(rec, dict):
            rec = obj if isinstance(obj, dict) else None
        if not isinstance(rec, dict):
            continue

        d = dict(rec)
        d["stage"] = _normalize_stage(d.get("stage", d.get("stage_label", "")))
        d["backend"] = str(d.get("backend", "")).upper()
        d["ladder"] = str(d.get("ladder", "")).upper()
        flat.append(d)

    if not flat:
        raise RuntimeError("No usable records found in JSONL. Expected a 'record' dict per line.")

    df = pd.DataFrame(flat)

    for col in [
        "shots", "repeats", "ridge", "mfull_cond",
        "tv_vs_target", "l2_vs_target", "fid_vs_target",
        "tv_vs_sim", "l2_vs_sim", "fid_vs_sim",
    ]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def _pick_metric_columns(df: pd.DataFrame, which: str) -> Tuple[str, str, str]:
    which = which.lower()
    if which == "target":
        return ("tv_vs_target", "l2_vs_target", "fid_vs_target")
    if which == "sim":
        return ("tv_vs_sim", "l2_vs_sim", "fid_vs_sim")
    raise ValueError("which must be 'target' or 'sim'")


def _boxplot_by_stage(df: pd.DataFrame, metric_col: str, title: str, out_png: str, out_pdf: str) -> None:
    stages = ["L0", "L01", "FULL"]
    backends = ["SIM", "NMR"]

    data: List[List[float]] = []
    labels: List[str] = []
    colors: List[Optional[str]] = []

    for st in stages:
        for be in backends:
            sel = df[(df["stage"] == st) & (df["backend"] == be)]
            vals = sel[metric_col].dropna().tolist() if metric_col in sel.columns else []
            if len(vals) == 0:
                continue
            data.append(vals)
            labels.append(f"{st}-{be}")
            colors.append("blue" if be == "SIM" else None)

    if not data:
        print(f"[warn] No data for boxplot {metric_col}")
        return

    plt.figure()
    bp = plt.boxplot(data, labels=labels, patch_artist=True, showmeans=False)
    for patch, col in zip(bp["boxes"], colors):
        if col is not None:
            patch.set_facecolor(col)

    plt.xticks(rotation=45, ha="right")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.savefig(out_pdf)
    plt.close()


def _hist_full(df: pd.DataFrame, metric_col: str, title: str, out_png: str, out_pdf: str) -> None:
    full = df[df["stage"] == "FULL"]
    if full.empty or metric_col not in full.columns:
        print(f"[warn] No FULL data for histogram {metric_col}")
        return

    sim = full[full["backend"] == "SIM"][metric_col].dropna().to_numpy()
    nmr = full[full["backend"] == "NMR"][metric_col].dropna().to_numpy()

    if sim.size == 0 and nmr.size == 0:
        print(f"[warn] No FULL data values for histogram {metric_col}")
        return

    plt.figure()
    if sim.size > 0:
        plt.hist(sim, bins="auto", alpha=0.7, label="SIM", color="blue")
    if nmr.size > 0:
        plt.hist(nmr, bins="auto", alpha=0.7, label="NMR")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.savefig(out_pdf)
    plt.close()


def _build_full_summary(df: pd.DataFrame, target_name: str, which: str) -> List[dict]:
    out: List[dict] = []
    full = df[df["stage"] == "FULL"].copy()
    if full.empty:
        return out

    tvc, l2c, fidc = _pick_metric_columns(df, which)

    for be in ["SIM", "NMR"]:
        sel = full[full["backend"] == be]
        if sel.empty or tvc not in sel.columns:
            continue
        tv = float(sel[tvc].dropna().mean())
        l2 = float(sel[l2c].dropna().mean())
        fid = float(sel[fidc].dropna().mean())
        out.append({
            "comparison": f"{target_name} vs {be} (FULL) — mean over records",
            "tv": tv,
            "l2": l2,
            "fidelity": fid,
        })
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Path to artifacts/gr_log_*.jsonl")
    ap.add_argument("--outdir", default="artifacts", help="Output directory (default: artifacts)")
    ap.add_argument("--metrics_against", choices=["target", "sim"], default="target",
                    help="Which metric set to plot: vs target or vs sim-ideal")
    args = ap.parse_args()

    input_path = args.input
    outdir = args.outdir
    which = args.metrics_against

    _ensure_dir(outdir)
    tag = _tag_from_input(input_path)

    rows = _read_jsonl(input_path)
    df = _extract_records(rows)

    metrics_csv = os.path.join(outdir, f"gr_metrics_{tag}.csv")
    df.to_csv(metrics_csv, index=False)
    print("Saved metrics CSV:", metrics_csv)

    tv_col, l2_col, fid_col = _pick_metric_columns(df, which)

    _boxplot_by_stage(
        df, tv_col,
        title=f"Total variation distance by stage (metrics vs {which})",
        out_png=os.path.join(outdir, f"gr_boxplot_tv_{which}_{tag}.png"),
        out_pdf=os.path.join(outdir, f"gr_boxplot_tv_{which}_{tag}.pdf"),
    )
    _boxplot_by_stage(
        df, l2_col,
        title=f"L2 distance by stage (metrics vs {which})",
        out_png=os.path.join(outdir, f"gr_boxplot_l2_{which}_{tag}.png"),
        out_pdf=os.path.join(outdir, f"gr_boxplot_l2_{which}_{tag}.pdf"),
    )
    _boxplot_by_stage(
        df, fid_col,
        title=f"Classical fidelity by stage (metrics vs {which})",
        out_png=os.path.join(outdir, f"gr_boxplot_fid_{which}_{tag}.png"),
        out_pdf=os.path.join(outdir, f"gr_boxplot_fid_{which}_{tag}.pdf"),
    )

    _hist_full(
        df, tv_col,
        title=f"FULL stage TV histogram (metrics vs {which})",
        out_png=os.path.join(outdir, f"gr_hist_tv_{which}_{tag}.png"),
        out_pdf=os.path.join(outdir, f"gr_hist_tv_{which}_{tag}.pdf"),
    )
    _hist_full(
        df, l2_col,
        title=f"FULL stage L2 histogram (metrics vs {which})",
        out_png=os.path.join(outdir, f"gr_hist_l2_{which}_{tag}.png"),
        out_pdf=os.path.join(outdir, f"gr_hist_l2_{which}_{tag}.pdf"),
    )
    _hist_full(
        df, fid_col,
        title=f"FULL stage Fidelity histogram (metrics vs {which})",
        out_png=os.path.join(outdir, f"gr_hist_fid_{which}_{tag}.png"),
        out_pdf=os.path.join(outdir, f"gr_hist_fid_{which}_{tag}.pdf"),
    )

    summary_rows = _build_full_summary(df, target_name="Target", which=which)
    if summary_rows:
        summ_csv = os.path.join(outdir, f"gr_summary_{tag}.csv")
        summ_tex = os.path.join(outdir, f"gr_summary_{tag}.tex")
        write_summary_csv(summ_csv, summary_rows)
        write_summary_tex(
            summ_tex,
            caption=f"Grover–Rudolph (FULL stage) summary metrics (mean over logged records, metrics vs {which}).",
            label="tab:gr_full_summary",
            rows=summary_rows,
        )
        print("Saved summary CSV:", summ_csv)
        print("Saved summary LaTeX:", summ_tex)
    else:
        print("[warn] No FULL summary rows produced (no FULL stage records found).")

    print("Done. Outputs written to:", outdir)


if __name__ == "__main__":
    main()
