# experiments/01_bit_order_check.py
from __future__ import annotations

import sys
from pathlib import Path

# Make "gr/" importable even if this script is launched from inside "experiments/".
try:
    import gr  # noqa: F401
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from spinqit import Circuit, X

from gr import (
    STATES_3Q,
    run_sim_probs,
    run_nmr_probs_robust,
    print_distribution,
)

RUN_NMR = False  # set True to also run on hardware (requires SPINQ_* env vars)
SHOTS_SIM = 50_000
SHOTS_NMR = 2048


def build_x_circuit(qi: int, ensure_nmr_attrs: bool) -> Circuit:
    c = Circuit()
    q = c.allocateQubits(3)
    c << (X, q[qi])
    # identity-safe tail handled inside build_gr_circuit_3q; here we keep it minimal.
    # If your backend needs it for trivial circuits, wrap with GR circuit builder or
    # manually append Ry(eps),Ry(-eps). We typically find X alone is accepted.
    return c


def print_top(probs, title: str, top: int = 8) -> None:
    items = sorted(probs.items(), key=lambda kv: kv[1], reverse=True)[:top]
    print(f"\n{title}")
    for k, v in items:
        print(f"  {k}: {v:.6f}")
    dom = max(probs.items(), key=lambda kv: kv[1])
    print(f"Dominant state: {dom[0]}  (prob={dom[1]:.6f})")


def main() -> None:
    for qi in [0, 1, 2]:
        circ = build_x_circuit(qi, ensure_nmr_attrs=False)

        sim = run_sim_probs(circ, shots=SHOTS_SIM)
        print_top(sim, f"=== SIM: X applied to q{qi} ===")

        if RUN_NMR:
            nmr = run_nmr_probs_robust(circ, name=f"X_q{qi}", shots=SHOTS_NMR)
            print_top(nmr, f"=== NMR: X applied to q{qi} ===")


if __name__ == "__main__":
    main()
