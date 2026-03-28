# experiments/06_readout_mitigation_apply.py
from __future__ import annotations

import sys
from pathlib import Path
import argparse
import numpy as np

# Make "gr/" importable even if this script is launched from inside "experiments/".
try:
    import gr  # noqa: F401
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

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

SHOTS_SIM = 200000
SHOTS_NMR = 2048
REPEATS_NMR = 5

RIDGE = 1e-3
LAM_MIX = 0.3

ARTIFACTS_DIR = Path("artifacts")

# Convention note:
# Raw, simulated, and mitigated distributions are all interpreted in the same
# canonical state order used throughout this repository:
#   ['000', '001', '010', '011', '100', '101', '110', '111']
# In this workflow, the effective comparison convention is MSB for both
# simulator and hardware.



def mix(raw: dict, mitig: dict, lam: float) -> dict:
    keys = raw.keys()
    out = {k: (1.0 - lam) * raw[k] + lam * mitig[k] for k in keys}
    s = sum(out.values())
    return {k: out[k] / s for k in keys}


def find_latest_mfull(artifacts_dir: Path) -> Path:
    """
    Returns the most recently modified artifacts/Mfull_*.npy file.
    Raises FileNotFoundError if none exist.
    """
    candidates = sorted(
        artifacts_dir.glob("Mfull_*.npy"),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    if not candidates:
        raise FileNotFoundError(
            f"No Mfull_*.npy found under '{artifacts_dir}'. "
            "Run: python -m experiments.05_readout_calibration_8x8 first."
        )
    return candidates[0]


def maybe_update_latest_copy(src: Path, artifacts_dir: Path) -> None:
    """
    Optionally maintain a convenience copy artifacts/Mfull_latest.npy.
    This avoids editing scripts when you regenerate Mfull.
    """
    dst = artifacts_dir / "Mfull_latest.npy"
    try:
        # overwrite copy
        data = np.load(str(src))
        np.save(str(dst), data)
    except Exception:
        # If anything goes wrong, don't block the experiment.
        pass


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Apply readout mitigation to FULL GR run using latest Mfull.")
    ap.add_argument(
        "--mfull",
        type=str,
        default=None,
        help="Path to Mfull_*.npy. If omitted, picks latest artifacts/Mfull_*.npy.",
    )
    ap.add_argument(
        "--no-update-latest",
        action="store_true",
        help="Do not write/update artifacts/Mfull_latest.npy convenience copy.",
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    ARTIFACTS_DIR.mkdir(parents=True, exist_ok=True)

    target = target_from_prob8(PROB8)
    print_distribution("Target", target)

    # Resolve Mfull path
    if args.mfull is not None:
        mfull_path = Path(args.mfull)
        if not mfull_path.exists():
            raise FileNotFoundError(f"--mfull file not found: {mfull_path}")
    else:
        mfull_path = find_latest_mfull(ARTIFACTS_DIR)

    if not args.no_update_latest:
        maybe_update_latest_copy(mfull_path, ARTIFACTS_DIR)

    M = load_matrix(str(mfull_path))
    print(f"Loaded Mfull: {mfull_path}")
    print("Condition number:", float(np.linalg.cond(M)))

    circ_sim = build_gr_circuit_3q(PROB8, [1, 1, 1], None, depth="full", ladder=LADDER, ensure_nmr_attrs=False)
    circ_nmr = build_gr_circuit_3q(PROB8, [1, 1, 1], None, depth="full", ladder=LADDER, ensure_nmr_attrs=True)

    if RUN_SIM:
        sim = run_sim_probs(circ_sim, shots=SHOTS_SIM)
        print_distribution("SIM ideal — FULL", sim)

    if RUN_NMR:
        nmr = run_nmr_repeated_avg(circ_nmr, base_name=f"GR_FULL_{LADDER}", shots=SHOTS_NMR, repeats=REPEATS_NMR)
        print_distribution("NMR raw avg — FULL", nmr)

        mitig = mitigate_readout(nmr, M, ridge=RIDGE)
        print_distribution(f"NMR mitigated — FULL (ridge={RIDGE})", mitig)

        mixed = mix(nmr, mitig, lam=LAM_MIX)
        print_distribution(f"NMR mixed — FULL (lam={LAM_MIX})", mixed)

        print_compare(target, nmr, "Target", "NMR raw avg")
        print_compare(target, mitig, "Target", "NMR mitigated")
        print_compare(target, mixed, "Target", "NMR mixed")
    else:
        print("\nRUN_NMR=False: nothing else to do (mitigation needs NMR).")


if __name__ == "__main__":
    main()
