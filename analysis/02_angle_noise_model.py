#!/usr/bin/env python3
"""Does an independent per-gate angle error predict the measured A/B effects?

Section 8 of the manuscript entertains the simplest mechanical model of the
observed compilation sensitivity: every commanded rotation is delivered with
an angle error. An earlier draft argued that no model of this family could
produce an A/B difference at all, because the perturbation is exchangeable
over gates and the two ladders command identical angle multisets. That
argument is wrong -- a perturbation drawn at a given position propagates
through the CNOT pattern that follows it, and the two ladders differ in
exactly that pattern -- so the question is quantitative and has to be
simulated. This script simulates it.

Four configurations are evaluated, matching the October campaign:

    Control  = the compiled FULL circuit including the identity-safe tail
    Reduced  = the same circuit with the three identity-equivalent gates removed
    ladders A and B = the two Gray-code orderings of the level-2 ladder

with the tail placed as the instrument's own record has it: adjacent at the end
for ladder A, split by the final CNOT for ladder B.

Two error models are run, because the campaign's own angle measurement
distinguishes them (Section 8: the error at the level-1 gate is dominated by a
systematic offset that changes between sessions, not by independent
fluctuation):

    independent  epsilon drawn afresh for every gate of every run
    systematic   one epsilon per gate, shared by all four configurations
                 within a run, redrawn between runs

and for each the script reports the mean classical fidelity of each cell, the
A/B contrasts, and the ladder x intervention interaction, each with a Monte
Carlo standard error, against the measured values.

Usage:
    python3 02_angle_noise_model.py [--sigma 1,2,3.73,6] [--trials 20000]
                                    [--grtri PATH] [--json OUT.json]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

# measured values the model has to reproduce (October campaign, D1, FULL)
MEASURED = {
    "F control A": 0.9669, "F control B": 0.8481,
    "A-B, Control": +0.0395, "A-B, Reduced": -0.0229,
    "Reduction, ladder A": +0.0140, "Reduction, ladder B": +0.0765,
    "Ladder x reduction": +0.0625,
}
# run-to-run sample s.d. of the same paired differences over the 25 repetitions
MEASURED_SD = {
    "A-B, Control": 0.0558, "A-B, Reduced": 0.0490,
    "Reduction, ladder A": 0.0656, "Reduction, ladder B": 0.0466,
    "Ladder x reduction": 0.0773,
}
# r.m.s. angle error measured at the level-1 gate across sessions, in degrees
MEASURED_RMS_DEG = 3.73

D1 = np.array([0.05, 0.10, 0.15, 0.20, 0.20, 0.15, 0.10, 0.05])

TAIL_DEG = (0.06, 719.94)       # the identity-safe pair, summing to 720


# ---------------------------------------------------------------------------
# circuits
# ---------------------------------------------------------------------------

def import_grtri(explicit: str | None):
    """Import the verification package's exact simulator."""
    cands = []
    if explicit:
        cands.append(Path(explicit))
    cands += [HERE.parent.parent / "gr-triangulum-verification",
              HERE.parent.parent.parent / "gr-triangulum-verification"]
    for c in cands:
        if (c / "grtri" / "simulator.py").exists():
            sys.path.insert(0, str(c))
            import grtri                                   # noqa: F401
            return grtri
    raise SystemExit("grtri not found; pass --grtri /path/to/"
                     "gr-triangulum-verification")


def build_gate_list(grtri, ladder: str, reduced: bool) -> list[tuple]:
    """FULL circuit for D1, with the tail placed as the instrument records it.

    The shipped compiler emits 19 gates and appends the tail at the end for
    both ladders. The device record has it at positions 19-20 (0-indexed) for
    ladder A, and at 18 and 20 -- split by the final CNOT -- for ladder B.
    """
    angles = grtri.angles.gr_angles(D1)
    ops = list(grtri.simulator._gate_list(angles, "FULL", ladder,
                                          identity_tail=False))
    assert len(ops) == 19, f"expected 19 gates, got {len(ops)}"
    if reduced:
        # the three identity-equivalent gates: the zero-angle rotation inside
        # the level-2 ladder, and the two tail pulses (never added here)
        ops = [op for op in ops
               if not (op[0] == "ry" and abs(op[2] % 720.0) < 1e-9)]
        return ops
    t0 = ("ry", 0, TAIL_DEG[0])
    t1 = ("ry", 0, TAIL_DEG[1])
    if ladder == "A":
        return ops + [t0, t1]
    return ops[:-1] + [t0, ops[-1], t1]          # split by the final CNOT


def n_rotations(ops: list[tuple]) -> int:
    return sum(1 for op in ops if op[0] == "ry")


# ---------------------------------------------------------------------------
# noisy evaluation
# ---------------------------------------------------------------------------

def run_probs(grtri, ops: list[tuple], eps_deg: np.ndarray) -> np.ndarray:
    """Output distribution with eps_deg[k] added to the k-th Ry angle."""
    sim = grtri.simulator
    state = np.zeros(8, dtype=complex)
    state[0] = 1.0
    k = 0
    for op in ops:
        if op[0] == "ry":
            state = sim._apply_1q(state, sim._ry(op[2] + eps_deg[k]), op[1])
            k += 1
        elif op[0] == "x":
            state = sim._apply_1q(state, sim._X, op[1])
        else:
            state = sim._apply_cx(state, op[1], op[2])
    p = np.abs(state) ** 2
    return p / p.sum()


def fidelity(p: np.ndarray, q: np.ndarray) -> float:
    return float(np.sum(np.sqrt(p * q)) ** 2)


def simulate(grtri, sigma_deg: float, trials: int, model: str,
             rng: np.random.Generator) -> dict:
    """Mean fidelity of each cell, and the contrasts, under one error model."""
    circuits = {(cfg, lad): build_gate_list(grtri, lad, cfg == "Reduced")
                for cfg in ("Control", "Reduced") for lad in ("A", "B")}
    nmax = max(n_rotations(ops) for ops in circuits.values())

    acc = {k: np.empty(trials) for k in circuits}
    for t in range(trials):
        shared = rng.normal(0.0, sigma_deg, nmax) if model == "systematic" else None
        for key, ops in circuits.items():
            m = n_rotations(ops)
            eps = shared[:m] if shared is not None else rng.normal(0.0, sigma_deg, m)
            acc[key][t] = fidelity(run_probs(grtri, ops, eps), D1)

    def mse(x):                                  # Monte Carlo standard error
        return float(x.std(ddof=1) / np.sqrt(len(x)))

    cA, cB = acc[("Control", "A")], acc[("Control", "B")]
    rA, rB = acc[("Reduced", "A")], acc[("Reduced", "B")]
    series = {
        "A-B, Control": cA - cB,
        "A-B, Reduced": rA - rB,
        "Reduction, ladder A": rA - cA,
        "Reduction, ladder B": rB - cB,
        "Ladder x reduction": (rB - cB) - (rA - cA),
    }
    return {
        "sigma_deg": sigma_deg, "model": model, "trials": trials,
        "gates": {f"{c} {l}": n_rotations(o) for (c, l), o in circuits.items()},
        "cells": {f"F {c.lower()} {l}": {"mean": float(acc[(c, l)].mean()),
                                         "mcse": mse(acc[(c, l)])}
                  for c in ("Control", "Reduced") for l in ("A", "B")},
        "contrasts": {k: {"mean": float(v.mean()), "mcse": mse(v),
                          "sd": float(v.std(ddof=1))}
                      for k, v in series.items()},
    }


# ---------------------------------------------------------------------------

def report(results: list[dict]) -> None:
    keys = ["A-B, Control", "A-B, Reduced", "Reduction, ladder A",
            "Reduction, ladder B", "Ladder x reduction"]
    for model in ("independent", "systematic"):
        rs = [r for r in results if r["model"] == model]
        if not rs:
            continue
        print("\n" + "=" * 104)
        print(f"ANGLE-ERROR MODEL: {model.upper()}"
              f"     ({rs[0]['trials']} trials per cell per sigma)")
        print("=" * 104)
        head = "".join(f"{('s=%.2f' % r['sigma_deg']):>14s}" for r in rs)
        print(f"{'quantity':26s}{head}{'MEASURED':>14s}")
        print("-" * 104)
        for lab in ("F control A", "F control B"):
            row = "".join(f"{r['cells'][lab]['mean']:14.4f}" for r in rs)
            print(f"{lab:26s}{row}{MEASURED[lab]:14.4f}")
        print("-" * 104)
        for k in keys:
            row = "".join(f"{r['contrasts'][k]['mean']:+14.4f}" for r in rs)
            print(f"{k:26s}{row}{MEASURED[k]:+14.4f}")
        print(f"{'(Monte Carlo s.e.)':26s}"
              + "".join(f"{r['contrasts']['Ladder x reduction']['mcse']:14.5f}"
                        for r in rs))
        print("-" * 104)
        print("predicted run-to-run s.d. of the same differences, "
              "against the measured s.d. over 25 repetitions")
        for k in keys:
            row = "".join(f"{r['contrasts'][k]['sd']:14.5f}" for r in rs)
            print(f"{'  sd ' + k:26s}{row}{MEASURED_SD[k]:14.4f}")
        print("=" * 104)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sigma", default=f"1,2,{MEASURED_RMS_DEG},6,10",
                    help="comma-separated angle-error sigmas in degrees")
    ap.add_argument("--trials", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=20261007)
    ap.add_argument("--grtri", default=None,
                    help="path to the gr-triangulum-verification checkout")
    ap.add_argument("--json", default=None)
    a = ap.parse_args()

    grtri = import_grtri(a.grtri)
    import grtri.angles, grtri.simulator                    # noqa: F401

    # the ideal circuits must still be exact, or nothing below means anything
    for lad in ("A", "B"):
        for red in (False, True):
            ops = build_gate_list(grtri, lad, red)
            p = run_probs(grtri, ops, np.zeros(n_rotations(ops)))
            err = float(np.abs(p - D1).max())
            tag = "Reduced" if red else "Control"
            assert err < 1e-12, f"{tag} {lad} not exact: {err:.2e}"
            print(f"  ideal check  {tag:8s} {lad}  {n_rotations(ops):2d} Ry  "
                  f"max |p - target| = {err:.2e}")

    sigmas = [float(s) for s in a.sigma.split(",")]
    rng = np.random.default_rng(a.seed)
    results = [simulate(grtri, s, a.trials, m, rng)
               for m in ("independent", "systematic") for s in sigmas]
    report(results)

    if a.json:
        Path(a.json).write_text(json.dumps(results, indent=2))
        print(f"wrote {a.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
