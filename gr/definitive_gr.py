# src/definitive_gr.py (ANONYMIZED, package-based runner)
#
# Grover–Rudolph state preparation (3 qubits) using the `gr/` package:
#   - Ideal simulation verification (sanity check vs target)
#   - Optional execution on a real NMR device via SpinQ backend (credentials via env vars)
#   - Ancilla-free transpilation into {Ry, X, CNOT} only
#   - Staged runs: L0 / L01 / FULL
#   - Metrics: TV / L2 / classical fidelity + per-qubit marginals
#   - Optional 8×8 readout calibration + mitigation
#
# NMR mode requires environment variables:
#   SPINQ_IP, SPINQ_PORT, SPINQ_USER, SPINQ_PASS
#
# Two-line meaning:
# RUN_NMR = False  -> simulation-only mode (no hardware connection, no SPINQ_* needed)
# RUN_NMR = True   -> hardware mode (runs on SpinQ Triangulum NMR; requires SPINQ_IP/PORT/USER/PASS)

from __future__ import annotations

import numpy as np

from gr import (
    STATES_3Q,
    target_from_prob8,
    build_gr_circuit_3q,
    run_sim_probs,
    run_nmr_repeated_avg,
    print_distribution,
    print_marginals,
    print_compare,
    calibrate_readout_matrix_8x8,
    mitigate_readout,
)

# ===================== USER CONFIG =====================

# Target weights (8 entries for 3 qubits): order is |q0 q1 q2> = 000..111 (MSB->LSB)
PROB8 = [1, 2, 3, 4, 4, 3, 2, 1]

# Execution knobs
RUN_SIM = True
RUN_NMR = False  # <-- set True to run on hardware (requires SPINQ_* env vars)

SHOTS_SIM = 200_000
SHOTS_NMR = 2048
REPEATS_NMR = 5

# Robustness/cooldown knobs for NMR are configured inside gr/backends.py defaults.

# Optional: minimal Ry gain calibration is NOT included in this runner (keep it minimal/clean).
# You can add it later if needed by extending gr/ with a calibration helper.

# Optional: hardware command clipping (try None first; enable if you see collapse/instability)
CLIP_THETA_CMD = None  # e.g. 0.95 * 3.141592653589793

# Choose level-2 CX ladder ordering: "A" or "B"
UCRY_LADDER = "B"

# Optional: readout mitigation (8×8). Needs NMR; takes extra time.
DO_READOUT_MITIGATION = False
SHOTS_RO = 4096
RIDGE = 1e-3

# ===================== INTERNALS =====================


def run_stage(depth: str, target_ref: dict, Mfull: np.ndarray | None) -> None:
    depth_u = depth.upper()

    # Ideal circuit for simulation (no scaling, no clipping, no tail)
    circ_sim = build_gr_circuit_3q(
        PROB8,
        k_gain=[1.0, 1.0, 1.0],
        clip_cmd=None,
        depth=depth,
        ladder=UCRY_LADDER,
        ensure_nmr_attrs=False,
    )

    # Hardware circuit for NMR (optional clipping + identity-safe tail)
    circ_nmr = build_gr_circuit_3q(
        PROB8,
        k_gain=[1.0, 1.0, 1.0],  # keep [1,1,1] for a clean baseline
        clip_cmd=CLIP_THETA_CMD,
        depth=depth,
        ladder=UCRY_LADDER,
        ensure_nmr_attrs=True,
    )

    print("\n" + "=" * 90)
    print(f"STAGE {depth_u}  (ladder {UCRY_LADDER})")
    print("=" * 90)

    sim_probs = None
    if RUN_SIM:
        sim_probs = run_sim_probs(circ_sim, shots=SHOTS_SIM)
        print_distribution(f"SIM ideal — {depth_u}", sim_probs)
        print_marginals(f"SIM — {depth_u}", sim_probs)

    nmr_probs = None
    nmr_mitig = None
    if RUN_NMR:
        nmr_probs = run_nmr_repeated_avg(
            circ_nmr,
            base_name=f"GR_{depth_u}_{UCRY_LADDER}",
            shots=SHOTS_NMR,
            repeats=REPEATS_NMR,
        )
        print_distribution(f"NMR raw avg — {depth_u}", nmr_probs)
        print_marginals(f"NMR raw — {depth_u}", nmr_probs)

        if Mfull is not None:
            nmr_mitig = mitigate_readout(nmr_probs, Mfull, ridge=RIDGE)
            print_distribution(f"NMR mitigated — {depth_u} (ridge={RIDGE})", nmr_mitig)
            print_marginals(f"NMR mitigated — {depth_u}", nmr_mitig)

    # Comparisons
    if RUN_SIM and RUN_NMR and sim_probs is not None and nmr_probs is not None:
        print_compare(sim_probs, nmr_probs, name_a=f"SIM {depth_u}", name_b=f"NMR {depth_u} raw avg")
        if nmr_mitig is not None:
            print_compare(sim_probs, nmr_mitig, name_a=f"SIM {depth_u}", name_b=f"NMR {depth_u} mitigated")

    if depth.lower() == "full" and RUN_NMR and nmr_probs is not None:
        print_compare(target_ref, nmr_probs, name_a="Target", name_b="NMR FULL raw avg")
        if nmr_mitig is not None:
            print_compare(target_ref, nmr_mitig, name_a="Target", name_b="NMR FULL mitigated")


def main() -> None:
    # Target distribution
    target = target_from_prob8(PROB8)
    print_distribution("Target distribution (normalized)", target)
    print_marginals("Target", target)

    # Optional: readout mitigation matrix
    if RUN_NMR and DO_READOUT_MITIGATION:
        print("\n=== Calibrating 8x8 readout matrix (Mfull) ===")
        Mfull = calibrate_readout_matrix_8x8(
            shots=SHOTS_RO,
            base_name="RO_FULL",
        )
        print("Done. Condition number:", float(np.linalg.cond(Mfull)))
    else:
        Mfull = None

    # Stage runs
    for depth in ["L0", "L01", "full"]:
        run_stage(depth=depth, target_ref=target, Mfull=Mfull)

    # Final sanity check: SIM FULL vs Target should be ~0 TV and ~1 fidelity
    if RUN_SIM:
        circ_check = build_gr_circuit_3q(
            PROB8,
            k_gain=[1.0, 1.0, 1.0],
            clip_cmd=None,
            depth="full",
            ladder=UCRY_LADDER,
            ensure_nmr_attrs=False,
        )
        sim_full = run_sim_probs(circ_check, shots=SHOTS_SIM)
        print("\n" + "=" * 90)
        print("SANITY CHECK (SIM FULL vs Target)")
        print("=" * 90)
        print_compare(target, sim_full, name_a="Target", name_b="SIM FULL")


if __name__ == "__main__":
    main()
