# experiments/02_trivial_suite.py
from __future__ import annotations

import sys
from pathlib import Path

# Make "gr/" importable even if this script is launched from inside "experiments/".
try:
    import gr  # noqa: F401
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import math
from spinqit import Circuit, Ry, CX

from gr import (
    run_sim_probs,
    run_nmr_probs_robust,
    print_distribution,
    print_marginals,
)

RUN_NMR = False
SHOTS_SIM = 100_000
SHOTS_NMR = 2048


def identity_safe() -> Circuit:
    c = Circuit()
    q = c.allocateQubits(3)
    c << (Ry, q[0], 1e-3)
    c << (Ry, q[0], -1e-3)
    return c


def ry_pi_over_2(qi: int) -> Circuit:
    c = Circuit()
    q = c.allocateQubits(3)
    c << (Ry, q[qi], math.pi / 2)
    c << (Ry, q[0], 1e-3)
    c << (Ry, q[0], -1e-3)
    return c


def bell_q0_q1() -> Circuit:
    c = Circuit()
    q = c.allocateQubits(3)
    c << (Ry, q[0], math.pi / 2)  # acts like H up to phase for our purposes
    c << (CX, (q[0], q[1]))
    c << (Ry, q[0], 1e-3)
    c << (Ry, q[0], -1e-3)
    return c


def run_case(name: str, circ: Circuit) -> None:
    sim = run_sim_probs(circ, shots=SHOTS_SIM)
    print_distribution(f"SIM: {name}", sim)
    print_marginals(f"SIM: {name}", sim)

    if RUN_NMR:
        nmr = run_nmr_probs_robust(circ, name=name, shots=SHOTS_NMR)
        print_distribution(f"NMR: {name}", nmr)
        print_marginals(f"NMR: {name}", nmr)


def main() -> None:
    run_case("IDENTITY_SAFE", identity_safe())
    for qi in [0, 1, 2]:
        run_case(f"RY_PI_OVER_2_q{qi}", ry_pi_over_2(qi))
    run_case("BELL_q0_q1", bell_q0_q1())


if __name__ == "__main__":
    main()
