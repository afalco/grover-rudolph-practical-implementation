# experiments/00_smoke_sim_only.py
from __future__ import annotations

from gr import (
    target_from_prob8,
    build_gr_circuit_3q,
    run_sim_probs,
    print_distribution,
    print_marginals,
    print_compare,
)

PROB8 = [1, 2, 3, 4, 4, 3, 2, 1]
SHOTS_SIM = 200_000
LADDER = "B"   # should not matter in ideal simulation if implementation is correct


def main() -> None:
    target = target_from_prob8(PROB8)
    print_distribution("Target (normalized)", target)
    print_marginals("Target", target)

    circ = build_gr_circuit_3q(
        PROB8, k_gain=[1, 1, 1], clip_cmd=None, depth="full",
        ladder=LADDER, ensure_nmr_attrs=False,
    )
    sim = run_sim_probs(circ, shots=SHOTS_SIM)
    print_distribution(f"SIM FULL (ladder {LADDER})", sim)
    print_marginals("SIM FULL", sim)

    print_compare(target, sim, "Target", "SIM FULL")


if __name__ == "__main__":
    main()
