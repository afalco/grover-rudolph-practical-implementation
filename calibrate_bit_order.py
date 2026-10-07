#!/usr/bin/env python3
"""
calibrate_bit_order.py -- measure the measurement bit-order convention of the backend.

Why this exists
---------------
Every probability vector in this repository, and every metric computed from one,
assumes the canonical order

    bitstring  b0 b1 b2   <->   qubits (q0, q1, q2)        [MSB -> LSB]

Until it is measured, that is an assumption. This utility measures it: it applies a
single X to each qubit in turn and records which position of the reported bitstring
carries the excitation. Three one-gate circuits settle the convention.

It is deliberately read-only with respect to configuration. It does not edit any
file; it prints the `SPINQ_BITORDER` value you should export and writes a JSON
report you can deposit alongside the campaign data.

Important: the measurement is taken on the RAW backend output. `gr.backends`
applies the `SPINQ_BITORDER` remap on the way out, which would mask exactly what we
are trying to measure, so that variable is neutralised for the duration of the run
and its prior value is recorded in the report.

Usage
-----
    python calibrate_bit_order.py --backend sim --shots 1024 --outdir artifacts

    python calibrate_bit_order.py \\
        --backend triangulum \\
        --ip <TRIANGULUM_IP> --port 55444 \\
        --account <ACCOUNT> --password <PASSWORD> \\
        --shots 1024 --outdir artifacts --unit-label triangulum-ceu

PowerShell, same call (backtick is the continuation character):

    python calibrate_bit_order.py `
        --backend triangulum `
        --ip <TRIANGULUM_IP> --port 55444 `
        --account <ACCOUNT> --password <PASSWORD> `
        --shots 1024 --outdir artifacts --unit-label triangulum-ceu

Credentials may equally be supplied through SPINQ_IP / SPINQ_PORT / SPINQ_USER /
SPINQ_PASS, in which case the corresponding flags can be omitted.

Exit status
-----------
    0  conclusive
    2  usage or connection error
    3  inconclusive -- the three tests do not define a permutation (see the report)
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import gr  # noqa: F401
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parent))

from spinqit import Circuit, Ry, X

from gr import STATES_3Q, run_sim_probs, run_nmr_probs_robust

EPS_TAIL = 1e-3
DOMINANCE_MIN = 0.30   # the flipped state must at least be the clear mode


# --------------------------------------------------------------------------- #
# Circuits                                                                     #
# --------------------------------------------------------------------------- #

def build_x_probe(qi: int, with_tail: bool = True) -> Circuit:
    """A single X on qubit qi, optionally followed by the identity-safe tail on q0.

    The tail is present by default because some NMR backends reject circuits that
    are 'too trivial'. It is Ry(+eps) Ry(-eps), an identity up to numerical noise,
    and it does not touch the excitation created by the X unless qi == 0, where its
    effect on the population is of order eps^2/4 ~ 2.5e-7.
    """
    c = Circuit()
    q = c.allocateQubits(3)
    c << (X, q[qi])
    if with_tail:
        c << (Ry, q[0], EPS_TAIL)
        c << (Ry, q[0], -EPS_TAIL)
    return c


# --------------------------------------------------------------------------- #
# Inference                                                                    #
# --------------------------------------------------------------------------- #

def dominant_state(probs: Dict[str, float]) -> Tuple[str, float]:
    k, v = max(probs.items(), key=lambda kv: kv[1])
    return k, float(v)


def flipped_position(bitstring: str) -> Optional[int]:
    """Index of the single '1' in the reported string, or None if not unique."""
    ones = [i for i, ch in enumerate(bitstring) if ch == "1"]
    return ones[0] if len(ones) == 1 else None


def infer_spec(positions: List[Optional[int]]) -> Tuple[Optional[str], str]:
    """
    positions[i] = index of the reported bit that lit up when X was applied to qi.

    The remap that carries backend order into canonical order is

        new[i] = old[positions[i]]

    which is exactly the permutation string consumed by gr.backends. Returns
    (spec, human-readable verdict); spec is None when the result is not a
    permutation.
    """
    if any(p is None for p in positions):
        return None, "at least one probe did not produce a unique excited bit"
    if sorted(positions) != [0, 1, 2]:
        return None, f"positions {positions} are not a permutation of (0,1,2)"
    spec = "".join(str(p) for p in positions)
    if spec == "012":
        return spec, "canonical: the backend already reports MSB->LSB (q0 q1 q2)"
    if spec == "210":
        return spec, "reversed: the backend reports LSB->MSB (q2 q1 q0)"
    return spec, f"non-trivial permutation: canonical[i] = backend[{spec}[i]]"


def export_value(spec: str) -> str:
    if spec == "012":
        return "MSB->LSB"
    if spec == "210":
        return "LSB->MSB"
    return spec


# --------------------------------------------------------------------------- #
# Main                                                                         #
# --------------------------------------------------------------------------- #

def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--backend", choices=("sim", "triangulum", "nmr"), default="sim")
    ap.add_argument("--shots", type=int, default=1024)
    ap.add_argument("--outdir", default="artifacts")
    ap.add_argument("--ip"); ap.add_argument("--port")
    ap.add_argument("--account"); ap.add_argument("--password")
    ap.add_argument("--unit-label", default="unspecified",
                    help="which physical machine this is, recorded in the report")
    ap.add_argument("--no-tail", action="store_true",
                    help="omit the identity-safe tail (probe is then a bare X)")
    args = ap.parse_args()

    hardware = args.backend in ("triangulum", "nmr")

    if hardware:
        for flag, var in (("ip", "SPINQ_IP"), ("port", "SPINQ_PORT"),
                          ("account", "SPINQ_USER"), ("password", "SPINQ_PASS")):
            val = getattr(args, flag)
            if val:
                os.environ[var] = str(val)
        missing = [v for v in ("SPINQ_IP", "SPINQ_PORT", "SPINQ_USER", "SPINQ_PASS")
                   if not os.environ.get(v)]
        if missing:
            print("Missing connection parameters: " + ", ".join(missing))
            print("Supply them as flags or as environment variables.")
            return 2

    # Neutralise the remap so that we observe the backend's raw convention.
    prior_bitorder = os.environ.pop("SPINQ_BITORDER", None)

    print("=" * 74)
    print("BIT-ORDER CALIBRATION")
    print("=" * 74)
    print(f"  backend        : {args.backend}")
    print(f"  unit           : {args.unit_label}")
    print(f"  shots          : {args.shots}")
    print(f"  probe          : X on q_i" + ("" if args.no_tail else " + identity-safe tail"))
    print(f"  SPINQ_BITORDER : {prior_bitorder!r} (neutralised for this measurement)")
    print()

    results = []
    positions: List[Optional[int]] = []
    try:
        for qi in (0, 1, 2):
            circ = build_x_probe(qi, with_tail=not args.no_tail)
            if hardware:
                probs = run_nmr_probs_robust(circ, name=f"BITORDER_X_q{qi}", shots=args.shots)
            else:
                probs = run_sim_probs(circ, shots=args.shots)
            dom, p = dominant_state(probs)
            pos = flipped_position(dom)
            runner_up = sorted(probs.values(), reverse=True)[1] if len(probs) > 1 else 0.0
            clear = (p >= DOMINANCE_MIN) and (p > runner_up)
            if not clear:
                pos = None
            positions.append(pos)
            results.append({
                "qubit": qi, "dominant": dom, "dominant_prob": p,
                "runner_up_prob": float(runner_up),
                "flipped_position": pos,
                "probs": {s: float(probs.get(s, 0.0)) for s in STATES_3Q},
            })
            mark = "ok" if pos is not None else "UNCLEAR"
            print(f"  X on q{qi}  ->  dominant {dom}  (p={p:.4f}, next={runner_up:.4f})  [{mark}]")
    except Exception as exc:                              # noqa: BLE001
        print(f"\nBackend error: {exc}")
        if prior_bitorder is not None:
            os.environ["SPINQ_BITORDER"] = prior_bitorder
        return 2

    spec, verdict = infer_spec(positions)

    print("\n" + "-" * 74)
    print(f"  Verdict: {verdict}")
    if spec is not None:
        val = export_value(spec)
        print("\n  bash / zsh:")
        print(f"    export SPINQ_BITORDER={val}")
        print("\n  PowerShell (the quotes are required: unquoted, '>' redirects):")
        print(f'    $env:SPINQ_BITORDER = "{val}"')
        print("\n  PowerShell, persistent for this user:")
        print(f'    [Environment]::SetEnvironmentVariable("SPINQ_BITORDER", "{val}", "User")')
        if spec != "012":
            print("\n  NOTE: this is NOT the convention assumed throughout the repository.")
            print("  Any data already acquired without this setting must be re-examined.")
    else:
        print("\n  No SPINQ_BITORDER value can be recommended from this run.")
        print("  Check readout calibration (05_readout_calibration_8x8.py) and retry:")
        print("  a strongly biased readout can blur the single-excitation probe.")
    print("-" * 74)

    report = {
        "utc": datetime.now(timezone.utc).isoformat(),
        "backend": args.backend,
        "unit_label": args.unit_label,
        "shots": args.shots,
        "probe_includes_identity_tail": not args.no_tail,
        "spinq_bitorder_before_run": prior_bitorder,
        "canonical_state_order": STATES_3Q,
        "measurements": results,
        "inferred_spec": spec,
        "recommended_SPINQ_BITORDER": export_value(spec) if spec else None,
        "verdict": verdict,
        "conclusive": spec is not None,
    }
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
    out = outdir / f"bit_order_{args.backend}_{stamp}.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"\n  Report: {out}")

    if prior_bitorder is not None:
        os.environ["SPINQ_BITORDER"] = prior_bitorder

    return 0 if spec is not None else 3


if __name__ == "__main__":
    raise SystemExit(main())
