# experiments/05_readout_calibration_8x8.py
from __future__ import annotations

import numpy as np

from gr import calibrate_readout_matrix_8x8

from .utils_io import ensure_dir, save_matrix, timestamp

RUN_NMR = True  # this script is hardware-only
SHOTS_RO = 4096

OUTDIR = "artifacts"


def main() -> None:
    ensure_dir(OUTDIR)

    M = calibrate_readout_matrix_8x8(shots=SHOTS_RO, base_name="RO_FULL")
    cond = float(np.linalg.cond(M))
    print("Done. Condition number:", cond)

    t = timestamp()
    path_npy = f"{OUTDIR}/Mfull_{t}.npy"
    path_csv = f"{OUTDIR}/Mfull_{t}.csv"
    save_matrix(path_npy, M, path_csv=path_csv)
    print("Saved:", path_npy)
    print("Saved:", path_csv)


if __name__ == "__main__":
    main()
