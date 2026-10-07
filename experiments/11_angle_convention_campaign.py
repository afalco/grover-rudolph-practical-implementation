#!/usr/bin/env python3
"""
11_angle_convention_campaign.py -- the angle-representative experiment.

Why this exists
---------------
The instrument's own job record shows that the 700-run benchmark campaign of
1-2 April 2026 and the 100-run campaign of 5 October 2026 did NOT execute the
same compiled circuit, although both were produced from the same source for the
same target. They differ in two ways, and one of them is this:

    level-1 "undo" rotations, April  :  Ry(q1, +663.211 deg), Ry(q1, +686.789 deg)
    level-1 "undo" rotations, October:  Ry(q1,  -56.789 deg), Ry(q1,  -33.211 deg)

                 663.211 = 720 - 56.789        686.789 = 720 - 33.211

The April compiler emitted -theta as 720 deg - theta. The two commands are the
same unitary -- exactly, not approximately -- and the April convention is
systematic: it appears on every target present in the export (D0, D1, D4, D5,
D6), on both ladders, at the same two circuit positions. October never commands
a rotation larger than 90 deg.

On a liquid-state NMR processor the rotation angle sets the pulse area, so a
675 deg rotation is not the same physical operation as a -45 deg one even
though it implements the same gate. This script measures whether that matters.

The design
----------
A 2x2 factorial, convention x ladder, on one target, in one session:

    direct : the compiled FULL circuit, undo rotations commanded as -theta
    folded : the same circuit, undo rotations commanded as -theta + 720 deg

Nothing else differs. Same gate sequence, same gate count, same gate types,
same controls, same every other angle, same target, same session, interleaved
run by run and paired by repetition. This is a stronger comparison than the
ladder A/B one, because A and B are different gate sequences whereas these two
are the same sequence with one number written two ways.

What a null result would mean
-----------------------------
If the device internally wraps the commanded angle into (-180, 180], the two
cells will agree and the experiment returns a null. That null is itself the
answer: it would establish that the April convention cost nothing, and remove
the circuit difference as a confounder from the longitudinal comparison. Either
outcome settles a question the deposited data cannot.

The pre-flight prints the commanded angles of every cell and refuses to submit a
hardware job unless (i) every variant reproduces the analytic target on the
simulator and (ii) the direct and folded variants agree on the simulator to
numerical precision while differing in exactly the two intended angles.

Credentials (environment)
-------------------------
    SPINQ_IP, SPINQ_PORT, SPINQ_USER, SPINQ_PASS
    SPINQ_BITORDER=MSB->LSB      # validate first with calibrate_bit_order.py

PowerShell (quotes required on the bit-order value, or '>' redirects):

    $env:SPINQ_IP = "<ip>"; $env:SPINQ_PORT = "55444"
    $env:SPINQ_USER = "<account>"; $env:SPINQ_PASS = "<password>"
    $env:SPINQ_BITORDER = "MSB->LSB"

Typical use
-----------
    # 1. dry run: no hardware, validates the circuits and prints the angle table
    python experiments/11_angle_convention_campaign.py --backend sim --repeats 2

    # 2. the campaign: 25 repeats x 2 conventions x 2 ladders = 100 runs
    python experiments/11_angle_convention_campaign.py \
        --backend nmr --repeats 25 --shots 2048 \
        --unit-label "triangulum-ceu" --outdir artifacts/angle_convention

Interrupting with Ctrl-C is safe: rerunning the same command resumes. Export the
instrument job database afterwards, so the commanded angles have a device-side
record as well -- that record is what makes the experiment checkable.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
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

TAIL_RAD = math.radians(0.0572958)   # the identity-safe tail pair on q0
SIM_SHOTS = 200_000
SIM_TV_TOL = 5e-3                    # sampled simulator vs analytic target
PAIR_TV_TOL = 5e-3                   # direct vs folded, on the simulator

# The April convention adds two full turns. Ry has period 4*pi for a spin-1/2,
# so +4*pi is the representative that leaves even the global phase untouched;
# +2*pi gives the same unitary up to an unobservable overall sign.
TURN_RAD = 2.0 * math.pi


# --------------------------------------------------------------------------- #
# Circuit construction                                                         #
# --------------------------------------------------------------------------- #

def build_full(prob8: List[float], ladder: str, fold_turns: float = 0.0,
               with_tail: bool = True) -> Circuit:
    """
    The compiled FULL Grover-Rudolph circuit, written out explicitly so that the
    two level-1 undo rotations can be commanded in either representative.

    fold_turns = 0 reproduces the canonical compilation (verified against
    build_gr_circuit_3q in the pre-flight); fold_turns = 2 reproduces the April
    convention, -theta commanded as -theta + 720 deg.

    Everything except those two angles is identical between the two variants:
    same gates, same order, same controls, same remaining angles.
    """
    ladder_u = ladder.upper()
    if ladder_u not in ("A", "B"):
        raise ValueError("ladder must be 'A' or 'B'")
    shift = fold_turns * TURN_RAD

    ang = angles_3q_asin_child1(prob8)
    circ = Circuit()
    q = circ.allocateQubits(3)

    # --- level 0 -----------------------------------------------------------
    circ << (Ry, q[0], float(ang[0][()]))

    # --- level 1: q1 conditioned on q0 -------------------------------------
    # Each block is  Ry(+half) CX Ry(-half) CX , bracketed by X on q0 for the
    # q0 = 0 branch. The second rotation is the one whose representative the
    # experiment varies.
    for q0_value in (0, 1):
        half = 0.5 * float(ang[1][(q0_value,)])
        if q0_value == 0:
            circ << (X, q[0])
        circ << (Ry, q[1], half)
        circ << (CX, (q[0], q[1]))
        circ << (Ry, q[1], -half + shift)
        circ << (CX, (q[0], q[1]))
        if q0_value == 0:
            circ << (X, q[0])

    # --- level 2: q2 conditioned on (q0, q1) via the Gray-code ladder ------
    a0, a1, a2, a3 = ucry_coeffs_from_thetas(
        float(ang[2][(0, 0)]), float(ang[2][(0, 1)]),
        float(ang[2][(1, 0)]), float(ang[2][(1, 1)]),
        ladder=ladder_u,
    )
    controls = [q[1], q[0], q[1], q[0]] if ladder_u == "A" else [q[0], q[1], q[0], q[1]]
    for coeff, ctrl in zip((a0, a1, a2, a3), controls):
        circ << (Ry, q[2], float(coeff))
        circ << (CX, (ctrl, q[2]))

    # --- identity-safe tail on q0 ------------------------------------------
    if with_tail:
        circ << (Ry, q[0], TAIL_RAD)
        circ << (Ry, q[0], -TAIL_RAD)

    return circ


def _builder(fold_turns: float, with_tail: bool):
    def f(prob8, ladder):
        return build_full(prob8, ladder, fold_turns=fold_turns, with_tail=with_tail)
    return f


# `direct` and `folded` are the factorial. `canonical` is the library's own
# builder, carried along in the pre-flight as a reference only.
BUILDERS = {
    "direct": _builder(0.0, True),
    "folded": _builder(2.0, True),
    "direct_notail": _builder(0.0, False),
    "folded_notail": _builder(2.0, False),
}


def gate_list(circ: Circuit) -> Optional[List]:
    for attr in ("gates", "instructions", "ops", "_gates", "instruction_list"):
        cand = getattr(circ, attr, None)
        if cand is not None:
            try:
                return list(cand)
            except Exception:
                continue
    return None


def gate_summary(circ: Circuit) -> Dict[str, Optional[int]]:
    """Best-effort gate census. Never fatal: returns None where introspection fails."""
    out: Dict[str, Optional[int]] = {"n_gates": None, "n_ry": None, "n_cx": None, "n_x": None}
    seq = gate_list(circ)
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


def commanded_angles(prob8: List[float], ladder: str, fold_turns: float,
                     with_tail: bool) -> List[float]:
    """
    The Ry angles this variant commands, in degrees, in emission order.

    Derived from the same expressions build_full uses, so it reports what is
    sent rather than what the SDK object happens to expose.
    """
    shift = fold_turns * TURN_RAD
    ang = angles_3q_asin_child1(prob8)
    out = [float(ang[0][()])]
    for q0_value in (0, 1):
        half = 0.5 * float(ang[1][(q0_value,)])
        out += [half, -half + shift]
    a = ucry_coeffs_from_thetas(
        float(ang[2][(0, 0)]), float(ang[2][(0, 1)]),
        float(ang[2][(1, 0)]), float(ang[2][(1, 1)]),
        ladder=ladder.upper(),
    )
    out += [float(x) for x in a]
    if with_tail:
        out += [TAIL_RAD, -TAIL_RAD]
    return [math.degrees(x) for x in out]


# --------------------------------------------------------------------------- #
# Schedule                                                                     #
# --------------------------------------------------------------------------- #

def build_schedule(configs: List[str], ladders: List[str], repeats: int) -> List[Tuple[int, str, str]]:
    """
    Interleaved schedule. Within each repeat every (config, ladder) cell is
    visited once; the order is reversed on odd repeats so that no cell
    systematically occupies an early or a late slot.
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

def preflight(prob8: List[float], configs: List[str], ladders: List[str],
              with_tail: bool) -> Dict[Tuple[str, str], Dict[str, float]]:
    """
    Validate every variant before any hardware time is spent.

      1. each variant reproduces the analytic target on the simulator;
      2. the canonical library builder is reported alongside, as a reference;
      3. direct and folded agree on the simulator and differ in exactly the two
         intended angles -- the premise of the whole experiment.
    """
    target = target_from_prob8(prob8)
    print("=" * 78)
    print("PRE-FLIGHT -- simulator validation and the commanded-angle table")
    print("=" * 78)

    refs: Dict[Tuple[str, str], Dict[str, float]] = {}
    failures = []
    for config in configs:
        for ladder in ladders:
            circ = BUILDERS[config](prob8, ladder)
            probs = run_sim_probs(circ, shots=SIM_SHOTS)
            tv, l2, fid = tv_l2_fidelity(target, probs)
            gs = gate_summary(circ)
            gc = "" if gs["n_gates"] is None else (
                f"  gates={gs['n_gates']} (Ry={gs['n_ry']}, CX={gs['n_cx']}, X={gs['n_x']})")
            ok = tv <= SIM_TV_TOL
            print(f"  {config:14s} ladder {ladder}   TV(target,sim)={tv:.2e}  fid={fid:.6f}{gc}"
                  f"   {'OK' if ok else 'FAIL'}")
            if not ok:
                failures.append((config, ladder, tv))
            refs[(config, ladder)] = probs

    # the library's own builder, for reference only
    print()
    for ladder in ladders:
        try:
            ref = build_gr_circuit_3q(prob8, [1.0, 1.0, 1.0], None,
                                      depth="full", ladder=ladder, ensure_nmr_attrs=True)
            p = run_sim_probs(ref, shots=SIM_SHOTS)
            tv, _, _ = tv_l2_fidelity(target, p)
            gs = gate_summary(ref)
            print(f"  [reference] build_gr_circuit_3q ladder {ladder}: TV={tv:.2e}  "
                  f"gates={gs['n_gates']}")
        except Exception as exc:                                   # pragma: no cover
            print(f"  [reference] build_gr_circuit_3q ladder {ladder}: unavailable ({exc})")

    # the commanded angles
    print("\n  Commanded Ry angles (degrees, emission order)")
    for config in configs:
        fold = 2.0 if "folded" in config else 0.0
        for ladder in ladders:
            a = commanded_angles(prob8, ladder, fold, with_tail)
            print(f"    {config:14s} {ladder}: " + " ".join(f"{x:+10.4f}" for x in a))

    # the premise: identical unitary, different representative
    if "direct" in configs and "folded" in configs:
        print()
        for ladder in ladders:
            tv_pair, _, _ = tv_l2_fidelity(refs[("direct", ladder)], refs[("folded", ladder)])
            ad = commanded_angles(prob8, ladder, 0.0, with_tail)
            af = commanded_angles(prob8, ladder, 2.0, with_tail)
            differ = [i for i, (x, y) in enumerate(zip(ad, af)) if abs(x - y) > 1e-9]
            same_count = len(ad) == len(af)
            ok = tv_pair <= PAIR_TV_TOL and differ == [2, 4] and same_count
            print(f"  ladder {ladder}: direct vs folded  TV(sim,sim)={tv_pair:.2e}   "
                  f"angles differing at positions {differ}   {'OK' if ok else 'FAIL'}")
            if not ok:
                failures.append((f"direct/folded pair", ladder, tv_pair))

    if failures:
        print("\nPre-flight FAILED for: " + ", ".join(f"{c}/{l} ({t:.2e})" for c, l, t in failures))
        print("No hardware job has been submitted. Fix the circuit construction first.")
        raise SystemExit(2)
    print("\nAll variants reproduce the target, and the two conventions differ only in the\n"
          "two level-1 undo rotations.\n")
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
       "n_gates", "n_ry", "n_cx", "n_x", "undo_angle_deg", "started_utc", "finished_utc"]
)


def make_record(*, run_name, dist_id, ladder, repeat, config, shots, backend, unit_label,
                exp_probs, sim_probs, target, gs, undo_deg, started, finished, duration) -> Dict:
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
        "undo_angle_deg": undo_deg, "started_utc": started, "finished_utc": finished,
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
    Read an existing log for resumption. A run counts as done only if it was
    acquired on the SAME backend: simulator and hardware rows share run names,
    so without this check a dry run would mark hardware runs complete and they
    would silently never be acquired.
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
            if rec.get("backend") == backend:
                names.add(rec.get("run_name"))
            else:
                foreign.add(str(rec.get("backend")))
    return records, names, foreign


# --------------------------------------------------------------------------- #
# Main                                                                         #
# --------------------------------------------------------------------------- #

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--backend", choices=("sim", "nmr"), default="sim",
                    help="'sim' validates everything without touching the device (default)")
    ap.add_argument("--repeats", type=int, default=25, help="repeats per (config, ladder) cell")
    ap.add_argument("--shots", type=int, default=2048, help="shots per run; RECORDED in the output")
    ap.add_argument("--configs", default="direct,folded",
                    help="comma-separated; append '_notail' to drop the identity-safe tail")
    ap.add_argument("--ladders", default="A,B")
    ap.add_argument("--dist-id", default="D1", help="label written to the dist_id column")
    ap.add_argument("--prob8", default=None,
                    help="comma-separated 8 weights; default is D1 = 1,2,3,4,4,3,2,1")
    ap.add_argument("--unit-label", default="unspecified",
                    help="which physical machine this is, e.g. 'triangulum-ceu'")
    ap.add_argument("--outdir", default="artifacts/angle_convention")
    ap.add_argument("--tag", default=None, help="run-name prefix; defaults to '<dist>_ang'")
    ap.add_argument("--cooldown", type=float, default=2.0, help="seconds between hardware jobs")
    args = ap.parse_args()

    prob8 = [float(x) for x in args.prob8.split(",")] if args.prob8 else list(PROB8_D1)
    if len(prob8) != 8:
        print("--prob8 needs exactly 8 values")
        return 2
    configs = [c.strip() for c in args.configs.split(",") if c.strip()]
    ladders = [l.strip().upper() for l in args.ladders.split(",") if l.strip()]
    for c in configs:
        if c not in BUILDERS:
            print(f"unknown config '{c}'; choose from {sorted(BUILDERS)}")
            return 2
    tails = {not c.endswith("_notail") for c in configs}
    if len(tails) > 1:
        print("mixing tailed and tail-free configs in one campaign confounds the comparison; "
              "run them as separate campaigns")
        return 2
    with_tail = tails.pop()

    tag = args.tag or f"{args.dist_id}_ang"
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    jsonl_path = outdir / f"{tag}_runs.jsonl"
    csv_path = outdir / f"{tag}_runs.csv"

    target = target_from_prob8(prob8)
    refs = preflight(prob8, configs, ladders, with_tail)

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

    print("=" * 78)
    print(f"ANGLE-CONVENTION CAMPAIGN   target={args.dist_id}  backend={args.backend}  "
          f"unit={args.unit_label}")
    print(f"           configs={configs}  ladders={ladders}  repeats={args.repeats}")
    print(f"           shots={args.shots}   total runs={len(schedule)}   already done={len(done)}")
    if args.backend == "nmr":
        print(f"           estimated time ~{len(todo) * 2.0 / 60.0:.1f} h at ~2 min/run")
    print(f"           output: {jsonl_path}")
    print("=" * 78)

    if not todo:
        print("Nothing to do -- every run in the schedule is already in the log.")
        write_csv(records, csv_path)
        print(f"CSV rewritten: {csv_path}")
        return 0

    if args.backend == "nmr":
        missing = [v for v in ("SPINQ_IP", "SPINQ_PORT", "SPINQ_USER", "SPINQ_PASS")
                   if not os.environ.get(v)]
        if missing:
            print("Missing environment variables: " + ", ".join(missing))
            return 2
        if not os.environ.get("SPINQ_BITORDER"):
            print("WARNING: SPINQ_BITORDER is not set. Validate the bit order before trusting "
                  "these results (see calibrate_bit_order.py), then export SPINQ_BITORDER.")

    circuits = {(c, l): BUILDERS[c](prob8, l) for c in configs for l in ladders}
    gsums = {k: gate_summary(v) for k, v in circuits.items()}
    undo = {}
    for c in configs:
        fold = 2.0 if "folded" in c else 0.0
        for l in ladders:
            undo[(c, l)] = round(commanded_angles(prob8, l, fold, with_tail)[2], 4)

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
                    gs=gsums[(config, ladder)], undo_deg=undo[(config, ladder)],
                    started=started, finished=datetime.now(timezone.utc).isoformat(),
                    duration=dt,
                )
                log.write(json.dumps(rec) + "\n")
                log.flush()
                os.fsync(log.fileno())
                records.append(rec)
                write_csv(records, csv_path)
                print(f"[{i:4d}/{len(todo)}] {name:30s} fid={rec['fid_tgt']:.4f} "
                      f"TV={rec['tv_tgt']:.4f}  {dt:5.1f}s")
    except KeyboardInterrupt:
        print("\nInterrupted. Progress is saved; rerun the same command to resume.")
        write_csv(records, csv_path)
        return 130

    write_csv(records, csv_path)

    # ---- closing summary --------------------------------------------------
    print("\n" + "=" * 78)
    print("SUMMARY -- mean fidelity against the target")
    print("=" * 78)
    cell = {}
    for config in configs:
        for ladder in ladders:
            vals = [r["fid_tgt"] for r in records
                    if r.get("config") == config and r.get("ladder") == ladder]
            if vals:
                cell[(config, ladder)] = np.array(vals, float)
                print(f"  {config:14s} ladder {ladder}:  {np.mean(vals):.4f} "
                      f"+- {np.std(vals, ddof=1) if len(vals) > 1 else 0.0:.4f}   (n={len(vals)})")

    if {"direct", "folded"} <= set(configs):
        print("\n  Paired by repetition, folded - direct:")
        for ladder in ladders:
            d = cell.get(("direct", ladder))
            f = cell.get(("folded", ladder))
            if d is None or f is None or len(d) != len(f):
                continue
            diff = f - d
            m = float(diff.mean())
            se = float(diff.std(ddof=1) / np.sqrt(len(diff))) if len(diff) > 1 else 0.0
            print(f"    ladder {ladder}:  {m:+.4f}  [{m - 1.96 * se:+.4f}, {m + 1.96 * se:+.4f}]"
                  f"   (n={len(diff)} pairs)")
        da = cell.get(("direct", "A")); db = cell.get(("direct", "B"))
        fa = cell.get(("folded", "A")); fb = cell.get(("folded", "B"))
        if all(x is not None for x in (da, db, fa, fb)) and len(da) == len(db) == len(fa) == len(fb):
            inter = (fa - da) - (fb - db)
            m = float(inter.mean())
            se = float(inter.std(ddof=1) / np.sqrt(len(inter))) if len(inter) > 1 else 0.0
            print(f"\n  Convention x ladder interaction: {m:+.4f} "
                  f"[{m - 1.96 * se:+.4f}, {m + 1.96 * se:+.4f}]")

    print(f"\nRun-level log: {jsonl_path}")
    print(f"Flat table:    {csv_path}")
    print("\nAfter the session, export the instrument job database as well: the commanded")
    print("angles in that export are what make this experiment checkable by a referee.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
