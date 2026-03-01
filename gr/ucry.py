# gr/ucry.py
from __future__ import annotations

import math
from typing import Optional, Tuple

import numpy as np
from spinqit import Circuit, CX, Ry


def ladder_sign_matrix(ladder: str) -> np.ndarray:
    """
    Build the 4x4 sign matrix S such that S a = theta, where:
      a = (a0,a1,a2,a3) are the Ry angles in the ladder,
      theta = (theta00,theta01,theta10,theta11) are conditional rotations,
    and ladder defines the CNOT control order.

    Ladder A: controls [c1, c0, c1, c0]
    Ladder B: controls [c0, c1, c0, c1]
    """
    ladder = ladder.upper()
    if ladder not in ("A", "B"):
        raise ValueError("ladder must be 'A' or 'B'")

    # For the 3 CNOTs BETWEEN the 4 Ry's, we track which control toggles parity.
    # For ladder A: between a0-a1 uses c1, a1-a2 uses c0, a2-a3 uses c1.
    # For ladder B: between a0-a1 uses c0, a1-a2 uses c1, a2-a3 uses c0.
    cx_controls = [1, 0, 1, 0] if ladder == "A" else [0, 1, 0, 1]  # for the 4 CNOTs in total
    # We derive the 4 signs for (a0..a3) per control pattern by parity tracking.

    rows = []
    for c0, c1 in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        bits = [c0, c1]
        parity = 0
        signs = []
        signs.append(+1 if parity == 0 else -1)  # a0 always positive initially
        for k in range(3):  # toggles before a1, a2, a3
            ctrl = cx_controls[k]  # which control drives the CNOT
            if bits[ctrl] == 1:
                parity ^= 1
            signs.append(+1 if parity == 0 else -1)
        rows.append(signs)
    return np.array(rows, dtype=float)


def ucry_coeffs_from_thetas(
    theta00: float, theta01: float, theta10: float, theta11: float, ladder: str
) -> Tuple[float, float, float, float]:
    """
    Solve for (a0,a1,a2,a3) in the ladder so that it implements the conditional thetas.
    """
    S = ladder_sign_matrix(ladder)
    theta_vec = np.array([theta00, theta01, theta10, theta11], dtype=float)
    a = np.linalg.solve(S, theta_vec)
    return float(a[0]), float(a[1]), float(a[2]), float(a[3])


def _clip(x: float, clip_cmd: Optional[float]) -> float:
    if clip_cmd is None:
        return x
    if abs(x) > clip_cmd:
        return math.copysign(clip_cmd, x)
    return x


def apply_ucry_2ctrl(
    circ: Circuit,
    c0,
    c1,
    t,
    theta00_cmd: float,
    theta01_cmd: float,
    theta10_cmd: float,
    theta11_cmd: float,
    ladder: str,
    clip_cmd: Optional[float] = None,
) -> None:
    """
    Apply a 2-control uniformly-controlled Ry on target t controlled by (c0,c1),
    using only Ry + CX ladder. The theta*_cmd must already be "hardware commands"
    (e.g., scaled/clipped as desired).
    """
    a0, a1, a2, a3 = ucry_coeffs_from_thetas(
        theta00_cmd, theta01_cmd, theta10_cmd, theta11_cmd, ladder=ladder
    )

    ladder = ladder.upper()
    circ << (Ry, t, _clip(a0, clip_cmd))

    if ladder == "A":
        circ << (CX, (c1, t))
        circ << (Ry, t, _clip(a1, clip_cmd))
        circ << (CX, (c0, t))
        circ << (Ry, t, _clip(a2, clip_cmd))
        circ << (CX, (c1, t))
        circ << (Ry, t, _clip(a3, clip_cmd))
        circ << (CX, (c0, t))
    elif ladder == "B":
        circ << (CX, (c0, t))
        circ << (Ry, t, _clip(a1, clip_cmd))
        circ << (CX, (c1, t))
        circ << (Ry, t, _clip(a2, clip_cmd))
        circ << (CX, (c0, t))
        circ << (Ry, t, _clip(a3, clip_cmd))
        circ << (CX, (c1, t))
    else:
        raise ValueError("ladder must be 'A' or 'B'")
