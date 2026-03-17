# experiments/09_table_summary.py
#
# Build a compact summary table (CSV + LaTeX) from gr_log_*.jsonl
# and plot statewise probability bar charts: Target vs SIM vs NMR (optionally mitigated/mixed),
# plus delta bar charts (SIM-Target, NMR-Target).
#
# Usage:
#   python -m experiments.09_table_summary --input artifacts/gr_log_YYYYMMDD-HHMMSS.jsonl
#
# Outputs (in artifacts/ by default):
#   gr_table_summary_<tag>.csv
#   gr_table_summary_<tag>.tex
#   gr_probs_<tag>_<stage>.png/.pdf
#   gr_diff_<tag>_<stage>.png/.pdf
#
# Notes:
# - Expects JSONL lines written by experiments/08_full_gr_log_to_csv.py
# - Works even if mitigation was disabled (mitigated/mixed rows will be skipped).
# - "Histogram" here means bar chart of probabilities per basis state (3 qubits => 8 bars).

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt

# Make "gr/" importable even if launched from inside "experiments/".
import sys
try:
    import gr  # noqa: F401
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from gr import STATES_3Q, target_from_prob8  # STATES_3Q = ['000',...,'111']


PROB8_DEFAULT = [1, 2, 3, 4, 4, 3, 2, 1]


# ----------------------------- Metrics -----------------------------

def normalize(d: Dict[str, float]) -> Dict[str, float]:
    s = float(sum(d.values()))
    if s <= 0:
        return {k: 0.0 for k in STATES_3Q}
    return {k: float(d.get(k, 0.0)) / s for k in STATES_3Q}


def vec(d: Dict[str, float]) -> np.ndarray:
    d = normalize(d)
    return np.array([float(d.get(k, 0.0)) for k in STATES_3Q], dtype=float)


def tv_l2_fidelity(a: Dict[str, float], b: Dict[str, float]) -> Tuple[float, float, float]:
    p = vec(a)
    q = vec(b)
    tv = 0.5 * float(np.sum(np.abs(p - q)))
    l2 = float(np.sqrt(np.sum((p - q) ** 2)))
    fid = float((np.sum(np.sqrt(np.maximum(p, 0) * np.maximum(q, 0)))) ** 2)
    return tv, l2, fid


def mix(raw: Dict[str, float], mitig: Dict[str, float], lam: float) -> Dict[str, float]:
    raw = normalize(raw)
    mitig = normalize(mitig)
    out = {k: (1.0 - lam) * raw[k] + lam * mitig[k] for k in STATES_3Q}
    return normalize(out)


# ----------------------------- Parsing -----------------------------

@dataclass
class Row:
    stage: str                  # L0/L01/FULL
    variant: str                # raw/mitigated/mixed
    ref: str                    # vs_target / vs_sim
    n: int
    tv_mean: float
    tv_std: float
    l2_mean: float
    l2_std: float
    fid_mean: float
    fid_std: float
    tv_median: float
    l2_median: float
    fid_median: float


def _stage_norm(x: str) -> str:
    x = (x or "").strip().upper()
    if x == "FULL" or x == "FULL ":
        return "FULL"
    if x in ("L0", "L01"):
        return x
    # some logs may store 'full'
    if x.lower() == "full":
        return "FULL"
    return x


def _get_probs(obj: dict, key: str) -> Optional[Dict[str, float]]:
    v = obj.get(key, None)
    if isinstance(v, dict) and len(v) > 0:
        return normalize({k: float(v.get(k, 0.0)) for k in STATES_3Q})
    return None


def load_jsonl(path: Path) -> List[dict]:
    records = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            records.append(json.loads(line))
    return records


# ----------------------------- Aggregation -----------------------------

def summarize_metrics(values: List[Tuple[float, float, float]]) -> Tuple[float, float, float, float, float, float, float, float, float]:
    # returns mean/std/median for tv,l2,fid
    arr = np.array(values, dtype=float)  # shape (n,3)
    tv = arr[:, 0]; l2 = arr[:, 1]; fid = arr[:, 2]
    return (
        float(tv.mean()), float(tv.std(ddof=1) if len(tv) > 1 else 0.0), float(np.median(tv)),
        float(l2.mean()), float(l2.std(ddof=1) if len(l2) > 1 else 0.0), float(np.median(l2)),
        float(fid.mean()), float(fid.std(ddof=1) if len(fid) > 1 else 0.0), float(np.median(fid)),
    )


def build_summary_table(
    jsonl: List[dict],
    target: Dict[str, float],
    lam_mix: float,
) -> List[Row]:
    """
    Uses per-run NMR records (repeats==1) and SIM stage reference from the JSONL.
    Computes metrics vs Target and vs SIM(stage).
    """
    # Collect SIM reference per stage (ideal)
    sim_by_stage: Dict[str, Dict[str, float]] = {}
    for item in jsonl:
        rec = item.get("record", {})
        if rec.get("backend") == "SIM":
            st = _stage_norm(rec.get("stage", rec.get("stage_label", "")))
            pr = _get_probs(item, "probs_raw")
            if pr is not None:
                sim_by_stage[st] = pr

    # Collect per-run NMR distributions by stage
    # We only consider individual runs: backend==NMR and repeats==1
    by_group: Dict[Tuple[str, str, str], List[Tuple[float, float, float]]] = {}
    # key: (stage, variant, ref) where ref in {"vs_target","vs_sim"}

    for item in jsonl:
        rec = item.get("record", {})
        if rec.get("backend") != "NMR":
            continue
        if int(rec.get("repeats", 0)) != 1:
            continue

        st = _stage_norm(rec.get("stage", rec.get("stage_label", "")))
        sim_ref = sim_by_stage.get(st, None)

        raw = _get_probs(item, "probs_raw")
        mitig = _get_probs(item, "probs_mitig")
        mixed = None
        if mitig is not None and raw is not None:
            mixed = _get_probs(item, "probs_mixed") or mix(raw, mitig, lam=lam_mix)

        variants: List[Tuple[str, Optional[Dict[str, float]]]] = [
            ("raw", raw),
            ("mitigated", mitig),
            ("mixed", mixed),
        ]

        for variant_name, dist in variants:
            if dist is None:
                continue
            # vs target
            key_t = (st, variant_name, "vs_target")
            by_group.setdefault(key_t, []).append(tv_l2_fidelity(target, dist))
            # vs sim(stage)
            if sim_ref is not None:
                key_s = (st, variant_name, "vs_sim")
                by_group.setdefault(key_s, []).append(tv_l2_fidelity(sim_ref, dist))

    rows: List[Row] = []
    for (st, variant, ref), vals in sorted(by_group.items(), key=lambda x: (x[0][0], x[0][2], x[0][1])):
        tv_m, tv_s, tv_med, l2_m, l2_s, l2_med, f_m, f_s, f_med = summarize_metrics(vals)
        rows.append(Row(
            stage=st, variant=variant, ref=ref,
            n=len(vals),
            tv_mean=tv_m, tv_std=tv_s,
            l2_mean=l2_m, l2_std=l2_s,
            fid_mean=f_m, fid_std=f_s,
            tv_median=tv_med, l2_median=l2_med, fid_median=f_med,
        ))
    return rows


# ----------------------------- Output: CSV + LaTeX -----------------------------

def write_csv(path: Path, rows: List[Row]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "stage", "variant", "ref", "n",
            "tv_mean", "tv_std", "tv_median",
            "l2_mean", "l2_std", "l2_median",
            "fid_mean", "fid_std", "fid_median",
        ])
        for r in rows:
            w.writerow([
                r.stage, r.variant, r.ref, r.n,
                r.tv_mean, r.tv_std, r.tv_median,
                r.l2_mean, r.l2_std, r.l2_median,
                r.fid_mean, r.fid_std, r.fid_median,
            ])


def latex_escape(s: str) -> str:
    return (
        s.replace("\\", "\\textbackslash{}")
         .replace("_", "\\_")
         .replace("%", "\\%")
         .replace("&", "\\&")
         .replace("#", "\\#")
    )


def write_tex(path: Path, rows: List[Row], caption: str, label: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    def fmt_pm(m: float, s: float) -> str:
        return f"{m:.6g} $\\pm$ {s:.2g}"

    lines = []
    lines.append("\\begin{table}[t]")
    lines.append("\\centering")
    lines.append("\\small")
    lines.append("\\begin{tabular}{lllcccc}")
    lines.append("\\hline")
    lines.append("Stage & Variant & Ref & $n$ & TV & $\\ell_2$ & Fidelity \\\\")
    lines.append("\\hline")
    for r in rows:
        stage = latex_escape(r.stage)
        variant = latex_escape(r.variant)
        ref = latex_escape(r.ref.replace("vs_", "vs\\_"))
        lines.append(
            f"{stage} & {variant} & {ref} & {r.n:d} & "
            f"{fmt_pm(r.tv_mean, r.tv_std)} & {fmt_pm(r.l2_mean, r.l2_std)} & {fmt_pm(r.fid_mean, r.fid_std)} \\\\"
        )
    lines.append("\\hline")
    lines.append("\\end{tabular}")
    lines.append(f"\\caption{{{latex_escape(caption)}}}")
    lines.append(f"\\label{{{latex_escape(label)}}}")
    lines.append("\\end{table}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ----------------------------- Plots: Target vs SIM vs NMR + deltas -----------------------------

def avg_dist(dists: List[Dict[str, float]]) -> Dict[str, float]:
    if not dists:
        return {k: 0.0 for k in STATES_3Q}
    acc = {k: 0.0 for k in STATES_3Q}
    for d in dists:
        d = normalize(d)
        for k in STATES_3Q:
            acc[k] += d[k]
    for k in STATES_3Q:
        acc[k] /= len(dists)
    return normalize(acc)


def collect_for_plot(jsonl: List[dict], stage: str, lam_mix: float) -> Tuple[Optional[Dict[str, float]], List[Dict[str, float]], List[Dict[str, float]], List[Dict[str, float]]]:
    """
    Returns (sim_stage, raw_runs, mitig_runs, mixed_runs) for NMR per-run (repeats==1).
    """
    st = _stage_norm(stage)

    sim = None
    for item in jsonl:
        rec = item.get("record", {})
        if rec.get("backend") == "SIM":
            st2 = _stage_norm(rec.get("stage", rec.get("stage_label", "")))
            if st2 == st:
                sim = _get_probs(item, "probs_raw")
                break

    raw_runs: List[Dict[str, float]] = []
    mitig_runs: List[Dict[str, float]] = []
    mixed_runs: List[Dict[str, float]] = []

    for item in jsonl:
        rec = item.get("record", {})
        if rec.get("backend") != "NMR":
            continue
        if int(rec.get("repeats", 0)) != 1:
            continue
        st2 = _stage_norm(rec.get("stage", rec.get("stage_label", "")))
        if st2 != st:
            continue

        raw = _get_probs(item, "probs_raw")
        mitig = _get_probs(item, "probs_mitig")
        mixed = None
        if raw is not None and mitig is not None:
            mixed = _get_probs(item, "probs_mixed") or mix(raw, mitig, lam=lam_mix)

        if raw is not None:
            raw_runs.append(raw)
        if mitig is not None:
            mitig_runs.append(mitig)
        if mixed is not None:
            mixed_runs.append(mixed)

    return sim, raw_runs, mitig_runs, mixed_runs


def plot_probs_and_diffs(outdir: Path, tag: str, stage: str,
                         target: Dict[str, float],
                         sim: Optional[Dict[str, float]],
                         nmr_avg: Optional[Dict[str, float]]) -> None:
    """
    Two figures per stage:
      - probabilities bar chart: Target vs SIM vs NMR-avg
      - delta bar chart: (SIM-Target), (NMR-Target)
    """
    stage_u = _stage_norm(stage)
    outdir.mkdir(parents=True, exist_ok=True)

    t = vec(target)
    s = vec(sim) if sim is not None else None
    n = vec(nmr_avg) if nmr_avg is not None else None

    x = np.arange(len(STATES_3Q))
    width = 0.28

    # --- Probs plot ---
    plt.figure(figsize=(10, 3.5))
    plt.bar(x - width, t, width=width, label="Target")
    if s is not None:
        plt.bar(x, s, width=width, label="SIM ideal")
    if n is not None:
        plt.bar(x + width, n, width=width, label="NMR avg")
    plt.xticks(x, STATES_3Q)
    plt.ylabel("Probability")
    plt.title(f"Stage {stage_u}: Target vs SIM vs NMR (avg)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / f"gr_probs_{tag}_{stage_u}.png", dpi=200)
    plt.savefig(outdir / f"gr_probs_{tag}_{stage_u}.pdf")
    plt.close()

    # --- Diffs plot ---
    plt.figure(figsize=(10, 3.5))
    if s is not None:
        plt.bar(x - width/2, s - t, width=width, label="SIM − Target")
    if n is not None:
        plt.bar(x + width/2, n - t, width=width, label="NMR − Target")
    plt.axhline(0.0, linewidth=1.0)
    plt.xticks(x, STATES_3Q)
    plt.ylabel("Δ probability")
    plt.title(f"Stage {stage_u}: Differences vs Target")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / f"gr_diff_{tag}_{stage_u}.png", dpi=200)
    plt.savefig(outdir / f"gr_diff_{tag}_{stage_u}.pdf")
    plt.close()


# ----------------------------- Main -----------------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, type=str, help="Path to artifacts/gr_log_<timestamp>.jsonl")
    ap.add_argument("--outdir", default="artifacts", type=str, help="Output directory (default: artifacts)")
    ap.add_argument("--tag", default=None, type=str, help="Tag used for outputs (default: derived from input filename)")
    ap.add_argument("--lam-mix", default=0.3, type=float, help="Lambda for mixed = (1-lam)*raw + lam*mitig (default: 0.3)")
    ap.add_argument("--prob8", default=None, type=str, help="Optional PROB8 as comma-separated list of 8 numbers")
    ap.add_argument("--stages", default="L0,L01,FULL", type=str, help="Comma-separated stages to summarize/plot")
    args = ap.parse_args()

    in_path = Path(args.input)
    outdir = Path(args.outdir)

    if args.tag is None:
        # e.g., gr_log_YYYYMMDD-HHMMSS.jsonl -> YYYYMMDD-HHMMSS
        stem = in_path.stem
        tag = stem.replace("gr_log_", "")
    else:
        tag = args.tag

    if args.prob8 is None:
        prob8 = PROB8_DEFAULT
    else:
        xs = [float(x.strip()) for x in args.prob8.split(",")]
        if len(xs) != 8:
            raise ValueError("prob8 must have exactly 8 numbers")
        prob8 = xs

    target = target_from_prob8(prob8)
    jsonl = load_jsonl(in_path)

    # --- Summary table ---
    rows = build_summary_table(jsonl=jsonl, target=target, lam_mix=float(args.lam_mix))
    out_csv = outdir / f"gr_table_summary_{tag}.csv"
    out_tex = outdir / f"gr_table_summary_{tag}.tex"
    write_csv(out_csv, rows)
    write_tex(
        out_tex,
        rows,
        caption="Grover–Rudolph summary metrics per stage (mean $\\pm$ std over NMR runs).",
        label="tab:gr_summary_by_stage",
    )

    print("Saved summary CSV:", out_csv)
    print("Saved summary LaTeX:", out_tex)

    # --- Plots: Target vs SIM vs NMR(avg) + diffs ---
    stages = [s.strip() for s in args.stages.split(",") if s.strip()]
    for st in stages:
        sim, raw_runs, mitig_runs, mixed_runs = collect_for_plot(jsonl, stage=st, lam_mix=float(args.lam_mix))

        # Prefer NMR raw avg (over runs). If none, skip plot.
        nmr_avg = avg_dist(raw_runs) if raw_runs else None

        plot_probs_and_diffs(
            outdir=outdir,
            tag=tag,
            stage=st,
            target=target,
            sim=sim,
            nmr_avg=nmr_avg,
        )

    print("Saved plots: gr_probs_* and gr_diff_* in", outdir)


if __name__ == "__main__":
    main()
