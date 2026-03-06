# experiments/03_ladder_full_A_vs_B.py
from __future__ import annotations

import sys
from pathlib import Path

# Make "gr/" importable even if this script is launched from inside "experiments/".
try:
    import gr  # noqa: F401
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from gr import (
    target_from_prob8,
    build_gr_circuit_3q,
    run_sim_probs,
    run_nmr_repeated_avg,
    print_distribution,
    print_compare,
)

PROB8 = [1, 2, 3, 4, 4, 3, 2, 1]

RUN_SIM = True
RUN_NMR = False

SHOTS_SIM = 200_000
SHOTS_NMR = 2048
REPEATS_NMR = 5

CLIP_THETA_CMD = None  # keep None first


def main() -> None:
    target = target_from_prob8(PROB8)
    print_distribution("Target", target)

    # --- SIM ---
    if RUN_SIM:
        circA_sim = build_gr_circuit_3q(PROB8, [1, 1, 1], None, depth="full", ladder="A", ensure_nmr_attrs=False)
        circB_sim = build_gr_circuit_3q(PROB8, [1, 1, 1], None, depth="full", ladder="B", ensure_nmr_attrs=False)
        simA = run_sim_probs(circA_sim, shots=SHOTS_SIM)
        simB = run_sim_probs(circB_sim, shots=SHOTS_SIM)

        print_distribution("SIM FULL (ladder A)", simA)
        print_distribution("SIM FULL (ladder B)", simB)
        print_compare(simA, simB, "SIM FULL A", "SIM FULL B")
        print_compare(target, simA, "Target", "SIM FULL A")
        print_compare(target, simB, "Target", "SIM FULL B")

    # --- NMR ---
    if RUN_NMR:
        circA = build_gr_circuit_3q(PROB8, [1, 1, 1], CLIP_THETA_CMD, depth="full", ladder="A", ensure_nmr_attrs=True)
        circB = build_gr_circuit_3q(PROB8, [1, 1, 1], CLIP_THETA_CMD, depth="full", ladder="B", ensure_nmr_attrs=True)

        nmrA = run_nmr_repeated_avg(circA, base_name="GR_FULL_A", shots=SHOTS_NMR, repeats=REPEATS_NMR)
        nmrB = run_nmr_repeated_avg(circB, base_name="GR_FULL_B", shots=SHOTS_NMR, repeats=REPEATS_NMR)

        print_distribution("NMR FULL avg (ladder A)", nmrA)
        print_distribution("NMR FULL avg (ladder B)", nmrB)
        print_compare(nmrA, nmrB, "NMR FULL A avg", "NMR FULL B avg")
        print_compare(target, nmrA, "Target", "NMR FULL A avg")
        print_compare(target, nmrB, "Target", "NMR FULL B avg")


if __name__ == "__main__":
    main()
