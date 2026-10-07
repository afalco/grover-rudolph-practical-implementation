#!/usr/bin/env python3
"""
10_ladder_ab_campaign.py -- interleaved ladder A/B campaign on the SpinQ Triangulum.

Purpose
-------
Produce a run-level, fully deposited dataset for the three-configuration study of
the manuscript (Section 9), on a locally operated Triangulum unit.

It covers the two configurations that are reachable through the gate-level SDK:

  control : the compiled circuit exactly as the paper defines it, including the
            Ry(0) placeholders and the identity-safe tail, executed under the
            stock pulse calibration of the device.
  exp1    : the same logical circuit with the identity-equivalent gates removed
            (zero-angle Ry gates and the identity-safe tail), pulse calibration
            left unchanged.

Experiment 2 of the study -- the globally optimised pulse set derived from the
circuit unitary -- is NOT reachable from this SDK surface and is not attempted
here. See README notes.

What makes this different from the original 10-run study
--------------------------------------------------------
  * ladders A and B are INTERLEAVED, not run in blocks, and the within-repeat
    order alternates, so the A/B comparison is balanced in acquisition time;
  * every run is logged individually, with the shot count recorded explicitly;
  * before any hardware job is submitted, each circuit variant is checked on the
    simulator against the analytic target, so a pruning mistake cannot silently
    corrupt a hardware session;
  * the output schema is the one of runs_flat_v2.csv, so the existing
    verification scripts read it without modification.

Credentials (environment)
-------------------------
    SPINQ_IP, SPINQ_PORT, SPINQ_USER, SPINQ_PASS
    SPINQ_BITORDER=MSB->LSB      # validate first with calibrate_bit_order.py

PowerShell equivalent (quotes required on the bit-order value, or '>' redirects):

    $env:SPINQ_IP = "<ip>"; $env:SPINQ_PORT = "55444"
    $env:SPINQ_USER = "<account>"; $env:SPINQ_PASS = "<password>"
    $env:SPINQ_BITORDER = "MSB->LSB"

Typical use
-----------
    # 1. dry run: no hardware, validates the circuits and prints the schedule
    python experiments/10_ladder_ab_campaign.py --backend sim --repeats 2

    # 2. real campaign: 25 repeats x 2 configs x 2 ladders = 100 runs
    python experiments/10_ladder_ab_campaign.py \
        --backend nmr --repeats 25 --shots 2048 \
        --unit-label "triangulum-ceu" --outdir artifacts/ladder_ab

Interrupting with Ctrl-C is safe: rerunning the same command resumes.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import gr  # noqa: F401
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from spinqit import Circuit, CX, Ry, X

from gr import (
    STATES_3Q,
    target_from_prob8,
    angles_3q_asin_child1,
    build_gr_circuit_3q,
    ucry_coeffs_from_thetas,
    run_sim_probs,
    run_nmr_probs_robust,
    tv_l2_fidelity,
    per_qubit_marginals,
)

# D1 of the benchmark suite: (0.05, 0.10, 0.15, 0.20, 0.20, 0.15, 0.10, 0.05)
PROB8_D1 = [1, 2, 3, 4, 4, 3, 2, 1]

ZERO_TOL = 1e-9          # an Ry below this (radians) is an identity gate
SIM_SHOTS = 200_000      # reference simulator shots
SIM_TV_TOL = 5e-3        # sampled-simulator agreement with the analytic target


# --------------------------------------------------------------------------- #
# Circuit construction                                                         #
# --------------------------------------------------------------------------- #

def build_control(prob8: List[float], ladder: str) -> Circuit:
    """The paper's compiled FULL circuit, identity-equivalent gates retained."""
    return build_gr_circuit_3q(
        prob8, [1.0, 1.0, 1.0], None,
        depth="full", ladder=ladder, ensure_nmr_attrs=True,
    )


def build_exp1(prob8: List[float], ladder: str) -> Circuit:
    """
    The same logical circuit with identity-equivalent gates removed.

    Two kinds of gate are dropped, and only these two:
      * the identity-safe tail Ry(+eps) Ry(-eps) on q0;
      * any Ry whose commanded angle is zero to within ZERO_TOL, which includes
        the Ry(0) placeholders inside the level-2 UCRy ladder.

    The CNOT ladder is left intact: its four CNOTs carry the parity structure and
    are not identity-equivalent individually. A level-1 controlled rotation with a
    zero angle collapses entirely (its two CNOTs share a control and cancel), so
    that block is dropped as a whole.

    Logical equivalence with build_control is not assumed -- it is verified
    against the analytic target on the simulator before any hardware job runs.
    """
    ladder_u = ladder.upper()
    if ladder_u not in ("A", "B"):
        raise ValueError("ladder must be 'A' or 'B'")

    ang = angles_3q_asin_child1(prob8)
    circ = Circuit()
    q = circ.allocateQubits(3)

    # --- level 0 -----------------------------------------------------------
    th0 = float(ang[0][()])
    if abs(th0) > ZERO_TOL:
        circ << (Ry, q[0], th0)

    # --- level 1: q1 conditioned on q0 -------------------------------------
    for q0_value in (0, 1):
        theta = float(ang[1][(q0_value,)])
        if abs(theta) <= ZERO_TOL:
            continue                      # Ry(0) CX Ry(0) CX == identity
        half = 0.5 * theta
        if q0_value == 0:
            circ << (X, q[0])
        circ << (Ry, q[1], half)
        circ << (CX, (q[0], q[1]))
        circ << (Ry, q[1], -half)
        circ << (CX, (q[0], q[1]))
        if q0_value == 0:
            circ << (X, q[0])

    # --- level 2: q2 conditioned on (q0,q1) via the UCRy ladder -------------
    a0, a1, a2, a3 = ucry_coeffs_from_thetas(
        float(ang[2][(0, 0)]), float(ang[2][(0, 1)]),
        float(ang[2][(1, 0)]), float(ang[2][(1, 1)]),
        ladder=ladder_u,
    )
    controls = [q[1], q[0], q[1], q[0]] if ladder_u == "A" else [q[0], q[1], q[0], q[1]]
    for coeff, ctrl in zip((a0, a1, a2, a3), controls):
        if abs(coeff) > ZERO_TOL:
            circ << (Ry, q[2], float(coeff))
        circ << (CX, (ctrl, q[2]))

    return circ


BUILDERS = {"control": build_control, "exp1": build_exp1}


def gate_summary(circ: Circuit) -> Dict[str, Optional[int]]:
    """Best-effort gate census. Never fatal: returns None where introspection fails."""
    out: Dict[str, Optional[int]] = {"n_gates": None, "n_ry": None, "n_cx": None, "n_x": None}
    seq = None
    for attr in ("gates", "instructions", "ops", "_gates", "instruction_list"):
        cand = getattr(circ, attr, None)
        if cand is not None:
            try:
                seq = list(cand)
                break
            except Exception:
                continue
    if seq is None:
        return out
    names = []
    for g in seq:
        n = getattr(getattr(g, "gate", None), "label", None) or getattr(g, "label", None) or repr(g)
        names.append(str(n).lower())
    out["n_gates"] = len(names)
    out["n_ry"] = sum(1 for n in names if re.search(r"\bry\b", n))
    out["n_cx"] = sum(1 for n in names if re.search(r"\bcx\b|\bcnot\b", n))
    out["n_x"] = sum(1 for n in names if re.search(r"\bx\b", n))
    return out


# --------------------------------------------------------------------------- #
# Schedule                                                                     #
# --------------------------------------------------------------------------- #

def build_schedule(configs: List[str], ladders: List[str], repeats: int) -> List[Tuple[int, str, str]]:
    """
    Interleaved schedule. Within each repeat every (config, ladder) cell is visited
    once; the order is reversed on odd repeats so that no cell systematically
    occupies an early or a late slot.
    """
    cells = [(c, l) for c in configs for l in ladders]
    schedule: List[Tuple[int, str, str]] = []
    for r in range(1, repeats + 1):
        order = cells if (r % 2 == 1) else list(reversed(cells))
        schedule.extend((r, c, l) for c, l in order)
    return schedule


def run_name_of(tag: str, config: str, ladder: str, repeat: int) -> str:
    return f"{tag}_{config}_{ladder}_{repeat:03d}"


# --------------------------------------------------------------------------- #
# Pre-flight                                                                   #
# --------------------------------------------------------------------------- #

def preflight(prob8: List[float], configs: List[str], ladders: List[str]) -> Dict[Tuple[str, str], Dict[str, float]]:
    """
    Build and simulate every circuit variant, verify it reproduces the analytic
    target, and return the reference simulator distributions.

    Raises SystemExit if any variant fails -- no hardware time is spent on a
    circuit we cannot vouch for.
    """
    target = target_from_prob8(prob8)
    print("=" * 74)
    print("PRE-FLIGHT -- simulator validation of every circuit variant")
    print("=" * 74)
    refs: Dict[Tuple[str, str], Dict[str, float]] = {}
    failures = []
    for config in configs:
        for ladder in ladders:
            circ = BUILDERS[config](prob8, ladder)
            probs = run_sim_probs(circ, shots=SIM_SHOTS)
            tv, l2, fid = tv_l2_fidelity(target, probs)
            gs = gate_summary(circ)
            gc = "" if gs["n_gates"] is None else (
                f"  gates={gs['n_gates']} (Ry={gs['n_ry']}, CX={gs['n_cx']}, X={gs['n_x']})"
            )
            ok = tv <= SIM_TV_TOL
            print(f"  {config:8s} ladder {ladder}   TV(target,sim)={tv:.2e}  fid={fid:.6f}{gc}"
                  f"   {'OK' if ok else 'FAIL'}")
            if not ok:
                failures.append((config, ladder, tv))
            refs[(config, ladder)] = probs
    if failures:
        print("\nPre-flight FAILED for: " + ", ".join(f"{c}/{l} (TV={t:.2e})" for c, l, t in failures))
        print("No hardware job has been submitted. Fix the circuit construction first.")
        raise SystemExit(2)
    print("\nAll variants reproduce the target on the simulator.\n")
    return refs


# --------------------------------------------------------------------------- #
# Record assembly                                                              #
# --------------------------------------------------------------------------- #

FLAT_COLUMNS = (
    ["run_name", "dist_id", "stage", "ladder", "repeat", "duration_s",
     "fid_sim", "tv_sim", "l2_sim", "fid_tgt", "tv_tgt", "l2_tgt"]
    + [f"{src}_{s}" for s in STATES_3Q for src in ("exp", "sim")]
    + ["sim_q0_p0", "exp_q0_p0", "sim_q1_p0", "exp_q1_p0", "sim_q2_p0", "exp_q2_p0"]
    + ["config", "shots", "backend", "unit_label",
       "n_gates", "n_ry", "n_cx", "n_x", "started_utc", "finished_utc"]
)


def make_record(*, run_name, dist_id, ladder, repeat, config, shots, backend, unit_label,
                exp_probs, sim_probs, target, gs, started, finished, duration) -> Dict:
    tv_t, l2_t, fid_t = tv_l2_fidelity(target, exp_probs)
    tv_s, l2_s, fid_s = tv_l2_fidelity(sim_probs, exp_probs)
    m_exp = per_qubit_marginals(exp_probs)
    m_sim = per_qubit_marginals(sim_probs)
    rec = {
        "run_name": run_name, "dist_id": dist_id, "stage": "FULL", "ladder": ladder,
        "repeat": repeat, "duration_s": round(duration, 3),
        "fid_sim": fid_s, "tv_sim": tv_s, "l2_sim": l2_s,
        "fid_tgt": fid_t, "tv_tgt": tv_t, "l2_tgt": l2_t,
        "config": config, "shots": shots, "backend": backend, "unit_label": unit_label,
        "started_utc": started, "finished_utc": finished,
    }
    for s in STATES_3Q:
        rec[f"exp_{s}"] = float(exp_probs[s])
        rec[f"sim_{s}"] = float(sim_probs[s])
    for i in range(3):
        rec[f"exp_q{i}_p0"] = float(m_exp[i][0])
        rec[f"sim_q{i}_p0"] = float(m_sim[i][0])
    rec.update({k: gs.get(k) for k in ("n_gates", "n_ry", "n_cx", "n_x")})
    return rec


def write_csv(records: List[Dict], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FLAT_COLUMNS, extrasaction="ignore")
        w.writeheader()
        for r in records:
            w.writerow(r)


def load_done(jsonl: Path, backend: str) -> Tuple[List[Dict], set, set]:
    """
    Read an existing log for resumption.

    A run counts as already done only if it was acquired on the SAME backend.
    Simulator and hardware rows must never share a log: the run names are
    identical, so without this check a dry run would mark the corresponding
    hardware runs as complete and they would silently never be acquired.
    The caller refuses to continue when foreign backends are present.
    """
    if not jsonl.exists():
        return [], set(), set()
    records, names, foreign = [], set(), set()
    with jsonl.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            records.append(rec)
            rec_backend = rec.get("backend")
            if rec_backend == backend:
                names.add(rec.get("run_name"))
            else:
                foreign.add(str(rec_backend))
    return records, names, foreign


# --------------------------------------------------------------------------- #
# Main                                                                         #
# --------------------------------------------------------------------------- #

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--backend", choices=("sim", "nmr"), default="sim",
                    help="'sim' validates everything without touching the device (default)")
    ap.add_argument("--repeats", type=int, default=25, help="repeats per (config, ladder) cell")
    ap.add_argument("--shots", type=int, default=2048, help="shots per run; RECORDED in the output")
    ap.add_argument("--configs", default="control,exp1")
    ap.add_argument("--ladders", default="A,B")
    ap.add_argument("--dist-id", default="D1", help="label written to the dist_id column")
    ap.add_argument("--prob8", default=None,
                    help="comma-separated 8 weights; default is D1 = 1,2,3,4,4,3,2,1")
    ap.add_argument("--unit-label", default="unspecified",
                    help="which physical machine this is, e.g. 'triangulum-ceu'")
    ap.add_argument("--outdir", default="artifacts/ladder_ab")
    ap.add_argument("--tag", default=None, help="run-name prefix; defaults to the dist id")
    ap.add_argument("--cooldown", type=float, default=2.0, help="seconds between hardware jobs")
    args = ap.parse_args()

    prob8 = [float(x) for x in args.prob8.split(",")] if args.prob8 else list(PROB8_D1)
    if len(prob8) != 8:
        print("--prob8 needs exactly 8 values"); return 2
    configs = [c.strip() for c in args.configs.split(",") if c.strip()]
    ladders = [l.strip().upper() for l in args.ladders.split(",") if l.strip()]
    for c in configs:
        if c not in BUILDERS:
            print(f"unknown config '{c}'; choose from {sorted(BUILDERS)}"); return 2

    tag = args.tag or args.dist_id
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    jsonl_path = outdir / f"{tag}_ladder_ab_runs.jsonl"
    csv_path = outdir / f"{tag}_ladder_ab_runs.csv"

    target = target_from_prob8(prob8)
    refs = preflight(prob8, configs, ladders)

    schedule = build_schedule(configs, ladders, args.repeats)
    records, done, foreign = load_done(jsonl_path, args.backend)
    if foreign:
        print(f"\nREFUSING TO CONTINUE: {jsonl_path} already holds runs acquired on "
              f"{', '.join(sorted(foreign))}, and this invocation uses '{args.backend}'.")
        print("Simulator and hardware runs share the same run names, so mixing them in one")
        print("log would both corrupt the table and silently skip hardware acquisitions.")
        print("Use a different --outdir for this backend, or move the existing log aside.")
        return 2
    todo = [job for job in schedule if run_name_of(tag, job[1], job[2], job[0]) not in done]

    print("=" * 74)
    print(f"CAMPAIGN   target={args.dist_id}  backend={args.backend}  unit={args.unit_label}")
    print(f"           configs={configs}  ladders={ladders}  repeats={args.repeats}")
    print(f"           shots={args.shots}   total runs={len(schedule)}   already done={len(done)}")
    if args.backend == "nmr":
        print(f"           estimated time ~{len(todo) * 2.0 / 60.0:.1f} h at ~2 min/run")
    print(f"           output: {jsonl_path}")
    print("=" * 74)

    if not todo:
        print("Nothing to do -- every run in the schedule is already in the log.")
        write_csv(records, csv_path)
        print(f"CSV rewritten: {csv_path}")
        return 0

    if args.backend == "nmr":
        missing = [v for v in ("SPINQ_IP", "SPINQ_PORT", "SPINQ_USER", "SPINQ_PASS")
                   if not os.environ.get(v)]
        if missing:
            print("Missing environment variables: " + ", ".join(missing)); return 2
        if not os.environ.get("SPINQ_BITORDER"):
            print("WARNING: SPINQ_BITORDER is not set. Validate the bit order before trusting "
                  "these results (see calibrate_bit_order.py), then export SPINQ_BITORDER.")

    circuits = {(c, l): BUILDERS[c](prob8, l) for c in configs for l in ladders}
    gsums = {k: gate_summary(v) for k, v in circuits.items()}

    try:
        with jsonl_path.open("a", encoding="utf-8") as log:
            for i, (repeat, config, ladder) in enumerate(todo, start=1):
                name = run_name_of(tag, config, ladder, repeat)
                started = datetime.now(timezone.utc).isoformat()
                t0 = time.time()
                if args.backend == "nmr":
                    exp_probs = run_nmr_probs_robust(
                        circuits[(config, ladder)], name=name,
                        shots=args.shots, cooldown_s=args.cooldown,
                    )
                else:
                    exp_probs = run_sim_probs(circuits[(config, ladder)], shots=args.shots)
                dt = time.time() - t0
                rec = make_record(
                    run_name=name, dist_id=args.dist_id, ladder=ladder, repeat=repeat,
                    config=config, shots=args.shots, backend=args.backend,
                    unit_label=args.unit_label, exp_probs=exp_probs,
                    sim_probs=refs[(config, ladder)], target=target,
                    gs=gsums[(config, ladder)], started=started,
                    finished=datetime.now(timezone.utc).isoformat(), duration=dt,
                )
                log.write(json.dumps(rec) + "\n")
                log.flush()
                os.fsync(log.fileno())
                records.append(rec)
                write_csv(records, csv_path)
                print(f"[{i:4d}/{len(todo)}] {name:28s} fid={rec['fid_tgt']:.4f} "
                      f"TV={rec['tv_tgt']:.4f}  {dt:5.1f}s")
    except KeyboardInterrupt:
        print("\nInterrupted. Progress is saved; rerun the same command to resume.")
        write_csv(records, csv_path)
        return 130

    write_csv(records, csv_path)

    # ---- closing summary --------------------------------------------------
    print("\n" + "=" * 74)
    print("SUMMARY -- mean fidelity against the target")
    print("=" * 74)
    for config in configs:
        means = {}
        for ladder in ladders:
            vals = [r["fid_tgt"] for r in records
                    if r.get("config") == config and r.get("ladder") == ladder]
            if vals:
                means[ladder] = (float(np.mean(vals)), float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0, len(vals))
                print(f"  {config:8s} ladder {ladder}:  {means[ladder][0]:.4f} +- {means[ladder][1]:.4f}   (n={means[ladder][2]})")
        if len(means) == 2 and set(means) == {"A", "B"}:
            mean_a = {s: float(np.mean([r[f"exp_{s}"] for r in records
                                        if r.get("config") == config and r.get("ladder") == "A"]))
                      for s in STATES_3Q}
            mean_b = {s: float(np.mean([r[f"exp_{s}"] for r in records
                                        if r.get("config") == config and r.get("ladder") == "B"]))
                      for s in STATES_3Q}
            tv_ab, l2_ab, fid_ab = tv_l2_fidelity(mean_a, mean_b)
            print(f"  {config:8s} A vs B:      TV={tv_ab:.4f}  l2={l2_ab:.4f}  fid={fid_ab:.4f}")
    print(f"\nRun-level log: {jsonl_path}")
    print(f"Flat table:    {csv_path}")
    print("\nAfter the session, export the instrument job database as well, so the")
    print("campaign has a device-side record from day one.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
