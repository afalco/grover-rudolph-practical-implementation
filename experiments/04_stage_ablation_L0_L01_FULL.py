# experiments/04_stage_ablation_L0_L01_FULL.py
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
LADDER = "B"

RUN_SIM = True
RUN_NMR = False

SHOTS_SIM = 200_000
SHOTS_NMR = 2048
REPEATS_NMR = 5

CLIP_THETA_CMD = None  # set e.g. 0.95*pi if needed


def run_stage(stage: str, target: dict) -> None:
    stage_u = stage.upper()
    print("\n" + "=" * 90)
    print(f"STAGE {stage_u} (ladder {LADDER})")
    print("=" * 90)

    circ_sim = build_gr_circuit_3q(PROB8, [1, 1, 1], None, depth=stage, ladder=LADDER, ensure_nmr_attrs=False)
    sim = run_sim_probs(circ_sim, shots=SHOTS_SIM) if RUN_SIM else None
    if sim is not None:
        print_distribution(f"SIM ideal — {stage_u}", sim)

    if RUN_NMR:
        circ_nmr = build_gr_circuit_3q(PROB8, [1, 1, 1], CLIP_THETA_CMD, depth=stage, ladder=LADDER, ensure_nmr_attrs=True)
        nmr = run_nmr_repeated_avg(circ_nmr, base_name=f"GR_{stage_u}_{LADDER}", shots=SHOTS_NMR, repeats=REPEATS_NMR)
        print_distribution(f"NMR raw avg — {stage_u}", nmr)

        if sim is not None:
            print_compare(sim, nmr, f"SIM {stage_u}", f"NMR {stage_u} avg")
        if stage.lower() == "full":
            print_compare(target, nmr, "Target", "NMR FULL avg")


def main() -> None:
    target = target_from_prob8(PROB8)
    print_distribution("Target", target)

    for stage in ["L0", "L01", "full"]:
        run_stage(stage, target)

    if RUN_SIM:
        circ_full = build_gr_circuit_3q(PROB8, [1, 1, 1], None, depth="full", ladder=LADDER, ensure_nmr_attrs=False)
        sim_full = run_sim_probs(circ_full, shots=SHOTS_SIM)
        print("\n" + "=" * 90)
        print("SANITY CHECK (Target vs SIM FULL)")
        print("=" * 90)
        print_compare(target, sim_full, "Target", "SIM FULL")


if __name__ == "__main__":
    main()
