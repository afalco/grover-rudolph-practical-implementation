# experiments/06_readout_mitigation_apply.py
from __future__ import annotations

import sys
from pathlib import Path

# Make "gr/" importable even if this script is launched from inside "experiments/".
try:
    import gr  # noqa: F401
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from gr import (
    target_from_prob8,
    build_gr_circuit_3q,
    run_nmr_repeated_avg,
    run_sim_probs,
    print_distribution,
    print_compare,
    mitigate_readout,
)

from .utils_io import load_matrix

PROB8 = [1, 2, 3, 4, 4, 3, 2, 1]
LADDER = "B"

RUN_SIM = True
RUN_NMR = True

SHOTS_SIM = 200_000
SHOTS_NMR = 2048
REPEATS_NMR = 5

RIDGE = 1e-3
LAM_MIX = 0.3

# Path to a previously saved Mfull_*.npy
MFULL_PATH = "artifacts/Mfull_latest.npy"  # replace with your file


def mix(raw: dict, mitig: dict, lam: float) -> dict:
    keys = raw.keys()
    out = {k: (1.0 - lam) * raw[k] + lam * mitig[k] for k in keys}
    s = sum(out.values())
    return {k: out[k] / s for k in keys}


def main() -> None:
    target = target_from_prob8(PROB8)
    print_distribution("Target", target)

    M = load_matrix(MFULL_PATH)
    print("Loaded Mfull. Condition number:", float(np.linalg.cond(M)))

    circ_sim = build_gr_circuit_3q(PROB8, [1, 1, 1], None, depth="full", ladder=LADDER, ensure_nmr_attrs=False)
    circ_nmr = build_gr_circuit_3q(PROB8, [1, 1, 1], None, depth="full", ladder=LADDER, ensure_nmr_attrs=True)

    if RUN_SIM:
        sim = run_sim_probs(circ_sim, shots=SHOTS_SIM)
        print_distribution("SIM ideal — FULL", sim)

    nmr = run_nmr_repeated_avg(circ_nmr, base_name=f"GR_FULL_{LADDER}", shots=SHOTS_NMR, repeats=REPEATS_NMR)
    print_distribution("NMR raw avg — FULL", nmr)

    mitig = mitigate_readout(nmr, M, ridge=RIDGE)
    print_distribution(f"NMR mitigated — FULL (ridge={RIDGE})", mitig)

    mixed = mix(nmr, mitig, lam=LAM_MIX)
    print_distribution(f"NMR mixed — FULL (lam={LAM_MIX})", mixed)

    print_compare(target, nmr, "Target", "NMR raw avg")
    print_compare(target, mitig, "Target", "NMR mitigated")
    print_compare(target, mixed, "Target", "NMR mixed")


if __name__ == "__main__":
    main()
