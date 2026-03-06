# gr/ucry.py
from __future__ import annotations

import math
from typing import Optional, Tuple

import numpy as np
from spinqit import Circuit, CX, Ry


def ladder_sign_matrix(ladder: str) -> np.ndarray:
    """
    Build the 4×4 sign matrix S such that:

        S · a = θ,

    where:
      - a = (a0, a1, a2, a3) are the Ry angles used in the UCRy ladder,
      - θ = (θ00, θ01, θ10, θ11) are the desired conditional rotations on the target,
        indexed by (c0, c1) ∈ {0,1}².

    The ladder choice defines the order of CNOT controls in the standard ancilla-free
    2-control uniformly-controlled Ry decomposition:

      Ladder A (controls):  c1, c0, c1, c0
      Ladder B (controls):  c0, c1, c0, c1

    Note on counting CNOTs:
      - The ladder contains 4 CNOTs total (it “opens” and “closes”).
      - Only the first 3 CNOTs occur *between* successive Ry gates and therefore
        determine the sign pattern (parity toggles) for a1, a2, a3 relative to a0.
      - The final (4th) CNOT closes the ladder but does not precede an additional Ry,
        so it does not introduce a new sign in S.
    """
    ladder = ladder.upper()
    if ladder not in ("A", "B"):
        raise ValueError("ladder must be 'A' or 'B'")

    # Controls of the CNOTs, in the order they appear in the ladder:
    #   A: [c1, c0, c1, c0]
    #   B: [c0, c1, c0, c1]
    #
    # For the sign system we only need the *three* toggles that happen
    # before applying a1, a2, a3 (i.e., between the four Ry gates).
    cx_controls = [1, 0, 1, 0] if ladder == "A" else [0, 1, 0, 1]

    rows = []
    for c0, c1 in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        bits = [c0, c1]
        parity = 0
        signs = []

        # a0 contributes with the initial sign (+) before any toggle
        signs.append(+1 if parity == 0 else -1)

        # Update parity through the three toggles that occur before a1, a2, a3
        for k in range(3):
            ctrl = cx_controls[k]  # which control is used by this CNOT
            if bits[ctrl] == 1:
                parity ^= 1
            signs.append(+1 if parity == 0 else -1)

        rows.append(signs)

    return np.array(rows, dtype=float)


def ucry_coeffs_from_thetas(
    theta00: float, theta01: float, theta10: float, theta11: float, ladder: str
) -> Tuple[float, float, float, float]:
    """
    Solve for (a0, a1, a2, a3) such that the chosen ladder implements the desired
    conditional rotations (θ00, θ01, θ10, θ11).

    This is a purely classical step: build the sign matrix S(ladder) and solve:

        S(ladder) · a = θ.
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
    Apply an ancilla-free 2-control uniformly-controlled Ry (UCRy) on target t
    with controls (c0, c1), using only the native gate set {Ry, CX}.

    Inputs:
      - theta**_cmd are *hardware command angles* (already scaled and/or clipped
        upstream if you use gain calibration or pulse-length limits).
      - ladder ∈ {"A","B"} selects the CNOT control order.

    The circuit is:

      Ry(a0) on t
      CNOT(control=?, target=t)
      Ry(a1) on t
      CNOT(control=?, target=t)
      Ry(a2) on t
      CNOT(control=?, target=t)
      Ry(a3) on t
      CNOT(control=?, target=t)

    where (a0..a3) are computed by solving S(ladder)·a = θ.
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
