#!/usr/bin/env python3
"""Serial dependence in the paired differences of the October intervention
campaign.

The main text reports five paired contrasts over 25 repetitions. The paired
analysis tolerates arbitrary dependence *within* a repetition, which is what
the shot-count question puts in doubt; it still assumes the 25 differences are
independent, or at least exchangeable, *across* repetitions. A drift slow on
the scale of one repetition but not on the scale of the session would violate
that, and the acquisition order was counterbalanced rather than randomised, so
the assumption is not guaranteed by the design.

This script tests it. For each contrast it reports:

  * the iid paired interval the paper quotes, for reference;
  * the autocorrelation of the 25 differences at lags 1-5, with Bartlett
    bands, plus Ljung-Box and Durbin-Watson;
  * a trend in repetition index (OLS slope and Spearman);
  * three dependence-robust alternatives to the iid interval --- a Newey-West
    HAC standard error, a moving-block bootstrap, and a block sign-flip
    permutation test.

A contrast whose conclusion is unchanged under all three is not resting on the
independence assumption. One whose interval widens enough to cross zero is.

Usage:
    python3 01_serial_dependence.py [--runs PATH] [--trials N] [--seed S]
                                    [--json OUT.json]
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
DEFAULT_RUNS = HERE.parent / "data" / "ladder_ab" / "D1_ladder_ab_runs.csv"

# config label in the CSV -> name used in the paper
CONTROL, REDUCED = "control", "exp1"


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------

def load_cells(path: Path) -> dict[tuple[str, str], np.ndarray]:
    """Return {(config, ladder): fidelity series ordered by repetition}."""
    rows = list(csv.DictReader(open(path)))
    cells: dict[tuple[str, str], dict[int, float]] = {}
    for r in rows:
        if r["stage"] != "FULL":
            continue
        key = (r["config"], r["ladder"])
        cells.setdefault(key, {})[int(r["repeat"])] = float(r["fid_tgt"])
    out = {}
    for key, by_rep in cells.items():
        reps = sorted(by_rep)
        if reps != list(range(min(reps), max(reps) + 1)):
            raise SystemExit(f"repetitions not contiguous for {key}: {reps}")
        out[key] = np.array([by_rep[i] for i in reps], float)
    return out


def contrasts(cells) -> dict[str, np.ndarray]:
    """The five paired difference series the main text reports."""
    cA, cB = cells[(CONTROL, "A")], cells[(CONTROL, "B")]
    rA, rB = cells[(REDUCED, "A")], cells[(REDUCED, "B")]
    return {
        "A-B, Control": cA - cB,
        "A-B, Reduced": rA - rB,
        "Reduction, ladder A": rA - cA,
        "Reduction, ladder B": rB - cB,
        "Ladder x reduction": (rB - cB) - (rA - cA),
    }


# ---------------------------------------------------------------------------
# dependence diagnostics
# ---------------------------------------------------------------------------

def acf(x: np.ndarray, nlags: int) -> np.ndarray:
    """Sample autocorrelation at lags 1..nlags (denominator n, as is standard)."""
    x = x - x.mean()
    denom = float(x @ x)
    return np.array([float(x[k:] @ x[:-k]) / denom for k in range(1, nlags + 1)])


def ljung_box(x: np.ndarray, nlags: int) -> tuple[float, float]:
    n = len(x)
    r = acf(x, nlags)
    q = n * (n + 2) * sum(r[k] ** 2 / (n - k - 1) for k in range(nlags))
    return float(q), float(stats.chi2.sf(q, nlags))


def durbin_watson(x: np.ndarray) -> float:
    e = x - x.mean()
    return float(np.sum(np.diff(e) ** 2) / np.sum(e ** 2))


def newey_west_se(x: np.ndarray, lags: int | None = None) -> tuple[float, int]:
    """HAC standard error of the sample mean (Bartlett kernel)."""
    n = len(x)
    if lags is None:                       # Newey-West rule of thumb
        lags = int(np.floor(4 * (n / 100.0) ** (2.0 / 9.0)))
    e = x - x.mean()
    gamma0 = float(e @ e) / n
    var = gamma0
    for k in range(1, lags + 1):
        gk = float(e[k:] @ e[:-k]) / n
        var += 2.0 * (1.0 - k / (lags + 1.0)) * gk
    var = max(var, 1e-300)
    return float(np.sqrt(var / n)), lags


def moving_block_bootstrap(x: np.ndarray, block: int, trials: int,
                           rng: np.random.Generator) -> np.ndarray:
    """Bootstrap distribution of the mean under a moving-block resample."""
    n = len(x)
    nblocks = int(np.ceil(n / block))
    starts_max = n - block + 1
    out = np.empty(trials)
    for t in range(trials):
        starts = rng.integers(0, starts_max, size=nblocks)
        samp = np.concatenate([x[s:s + block] for s in starts])[:n]
        out[t] = samp.mean()
    return out


def block_signflip(x: np.ndarray, block: int, trials: int,
                   rng: np.random.Generator) -> float:
    """Two-sided p-value for mean zero, flipping the sign of whole blocks.

    Sign-flipping preserves the dependence structure inside a block, so the
    reference distribution is valid under within-block dependence of any form
    provided the series is symmetric about its (null) mean.
    """
    n = len(x)
    edges = list(range(0, n, block))
    obs = abs(x.mean())
    count = 0
    for _ in range(trials):
        signs = rng.choice([-1.0, 1.0], size=len(edges))
        y = x.copy()
        for s, e0 in zip(signs, edges):
            y[e0:e0 + block] *= s
        if abs(y.mean()) >= obs - 1e-15:
            count += 1
    return (count + 1) / (trials + 1)


# ---------------------------------------------------------------------------

def analyse(name: str, d: np.ndarray, trials: int,
            rng: np.random.Generator) -> dict:
    n = len(d)
    mean = float(d.mean())
    sd = float(d.std(ddof=1))
    se = sd / np.sqrt(n)
    tcrit = float(stats.t.ppf(0.975, n - 1))
    t_p = float(stats.ttest_1samp(d, 0.0).pvalue)

    nlags = min(5, n // 3)
    r = acf(d, nlags)
    band = 1.96 / np.sqrt(n)
    lb_q, lb_p = ljung_box(d, nlags)
    dw = durbin_watson(d)

    idx = np.arange(1, n + 1, dtype=float)
    slope, intercept, rval, trend_p, stderr = stats.linregress(idx, d)
    sp_rho, sp_p = stats.spearmanr(idx, d)

    nw_se, nw_lags = newey_west_se(d)
    z = float(stats.norm.ppf(0.975))

    block = max(2, int(round(n ** (1.0 / 3.0))))       # ~3 for n = 25
    boot = moving_block_bootstrap(d, block, trials, rng)
    blo, bhi = np.percentile(boot, [2.5, 97.5])
    perm_p = block_signflip(d, block, trials, rng)

    return {
        "name": name, "n": n, "mean": mean, "sd": sd,
        "iid": {"se": se, "lo": mean - tcrit * se, "hi": mean + tcrit * se,
                "p": t_p},
        "acf": {f"lag{k + 1}": float(r[k]) for k in range(nlags)},
        "acf_band": float(band),
        "ljung_box": {"Q": lb_q, "p": lb_p, "lags": nlags},
        "durbin_watson": dw,
        "trend": {"slope_per_rep": float(slope), "p": float(trend_p),
                  "spearman_rho": float(sp_rho), "spearman_p": float(sp_p)},
        "hac": {"se": nw_se, "lags": nw_lags,
                "lo": mean - z * nw_se, "hi": mean + z * nw_se,
                "inflation": nw_se / se},
        "block_bootstrap": {"block": block, "trials": trials,
                            "lo": float(blo), "hi": float(bhi)},
        "block_signflip_p": perm_p,
    }


def crosses_zero(lo: float, hi: float) -> bool:
    return lo <= 0.0 <= hi


def report(res: list[dict]) -> None:
    w = 22
    print("\n" + "=" * 100)
    print("PAIRED DIFFERENCES: DEPENDENCE DIAGNOSTICS  (n = %d repetitions)"
          % res[0]["n"])
    print("=" * 100)
    print(f"{'contrast':{w}s} {'mean':>9s} {'acf(1)':>8s} {'DW':>6s} "
          f"{'LB p':>8s} {'trend p':>8s} {'HAC/iid':>8s}")
    print("-" * 100)
    for r in res:
        print(f"{r['name']:{w}s} {r['mean']:+9.4f} {r['acf']['lag1']:+8.3f} "
              f"{r['durbin_watson']:6.2f} {r['ljung_box']['p']:8.3f} "
              f"{r['trend']['p']:8.3f} {r['hac']['inflation']:8.2f}")
    print("-" * 100)
    print("Bartlett band for a white-noise acf at 95%%: +/- %.3f"
          % res[0]["acf_band"])

    print("\n" + "=" * 100)
    print("INTERVALS: iid PAIRED vs DEPENDENCE-ROBUST")
    print("=" * 100)
    print(f"{'contrast':{w}s} {'iid 95% CI':>24s} {'HAC 95% CI':>24s} "
          f"{'block boot 95% CI':>24s}")
    print("-" * 100)
    for r in res:
        i, h, b = r["iid"], r["hac"], r["block_bootstrap"]
        print(f"{r['name']:{w}s} "
              f"{'[%+.4f,%+.4f]' % (i['lo'], i['hi']):>24s} "
              f"{'[%+.4f,%+.4f]' % (h['lo'], h['hi']):>24s} "
              f"{'[%+.4f,%+.4f]' % (b['lo'], b['hi']):>24s}")
    print("-" * 100)
    print(f"{'contrast':{w}s} {'p (paired t)':>14s} {'p (block flip)':>16s} "
          f"{'verdict':>34s}")
    print("-" * 100)
    for r in res:
        same = (crosses_zero(r['iid']['lo'], r['iid']['hi'])
                == crosses_zero(r['hac']['lo'], r['hac']['hi'])
                == crosses_zero(r['block_bootstrap']['lo'],
                                r['block_bootstrap']['hi']))
        verdict = "unchanged" if same else "CHANGES under dependence"
        print(f"{r['name']:{w}s} {r['iid']['p']:14.5f} "
              f"{r['block_signflip_p']:16.4f} {verdict:>34s}")
    print("=" * 100 + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default=str(DEFAULT_RUNS),
                    help="run-level CSV of the October campaign")
    ap.add_argument("--trials", type=int, default=20000,
                    help="bootstrap / permutation replicates (default 20000)")
    ap.add_argument("--seed", type=int, default=20261007)
    ap.add_argument("--json", default=None, help="write full results here")
    a = ap.parse_args()

    path = Path(a.runs)
    if not path.exists():
        raise SystemExit(f"run file not found: {path}")
    cells = load_cells(path)
    missing = [k for k in ((CONTROL, "A"), (CONTROL, "B"),
                           (REDUCED, "A"), (REDUCED, "B")) if k not in cells]
    if missing:
        raise SystemExit(f"missing cells: {missing}")

    rng = np.random.default_rng(a.seed)
    res = [analyse(name, d, a.trials, rng)
           for name, d in contrasts(cells).items()]
    report(res)

    if a.json:
        Path(a.json).write_text(json.dumps(res, indent=2))
        print(f"wrote {a.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
