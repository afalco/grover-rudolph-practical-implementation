# gr/circuit.py
from __future__ import annotations

import math
from typing import List, Optional

from spinqit import Circuit, CX, Ry, X

from .angles import angles_3q_asin_child1
from .ucry import apply_ucry_2ctrl


def add_identity_safe_tail(circ: Circuit, q, eps: float = 1e-3) -> None:
    """
    Append a tiny identity on q0: Ry(eps) then Ry(-eps).

    Rationale:
      Some NMR/IR backends exhibit edge cases for circuits that are "too trivial"
      (e.g., missing attributes / overly optimized DAG). This tail keeps the
      circuit non-empty without changing the logical unitary (up to numerical noise).
    """
    circ << (Ry, q[0], eps)
    circ << (Ry, q[0], -eps)


def scale_clip(theta: float, k: float, clip_cmd: Optional[float]) -> float:
    """
    Convert an *ideal* rotation angle theta into a *hardware command* angle theta_cmd.

    Model:
      hardware implements approximately: theta_eff ≈ k * theta_cmd
      => command: theta_cmd = theta / k

    Notes:
      - This is a minimal gain correction (amplitude scaling) applied on the *target qubit*.
      - If clip_cmd is provided, clipping is applied to the *command* magnitude |theta_cmd|
        (useful to avoid excessively long pulses when k<1).
    """
    theta_cmd = theta / (k if k > 1e-6 else 1.0)
    if clip_cmd is not None and abs(theta_cmd) > clip_cmd:
        theta_cmd = math.copysign(clip_cmd, theta_cmd)
    return theta_cmd


def apply_cry_1ctrl_decomp(circ: Circuit, c, t, theta_cmd: float, clip_cmd: Optional[float]) -> None:
    """
    Decompose a controlled-Ry(c -> t, theta_cmd) using only {Ry, CX}:

        Ry(theta/2) on t
        CX(c,t)
        Ry(-theta/2) on t
        CX(c,t)

    This avoids relying on nested parametric controlled gates in the compiler.
    """
    half = 0.5 * theta_cmd
    if clip_cmd is not None and abs(half) > clip_cmd:
        half = math.copysign(clip_cmd, half)

    circ << (Ry, t, half)
    circ << (CX, (c, t))
    circ << (Ry, t, -half)
    circ << (CX, (c, t))


def apply_patterned_cry_on_q0(
    circ: Circuit,
    q,
    target: int,
    theta_cmd: float,
    q0_value: int,
    clip_cmd: Optional[float],
) -> None:
    """
    Apply a controlled-Ry on q[target] conditioned on (q0 == q0_value).

    Implementation:
      - q0_value == 1: direct CRy with control q0
      - q0_value == 0: X on q0, then CRy, then X on q0 (0-control via X-sandwich)

    The CRy itself is decomposed into {Ry, CX} via apply_cry_1ctrl_decomp.
    """
    if q0_value not in (0, 1):
        raise ValueError("q0_value must be 0 or 1")

    if q0_value == 0:
        circ << (X, q[0])
    apply_cry_1ctrl_decomp(circ, c=q[0], t=q[target], theta_cmd=theta_cmd, clip_cmd=clip_cmd)
    if q0_value == 0:
        circ << (X, q[0])


def build_gr_circuit_3q(
    prob8: List[float],
    k_gain: List[float],
    clip_cmd: Optional[float],
    depth: str = "full",
    ladder: str = "B",
    ensure_nmr_attrs: bool = True,
    eps_tail: float = 1e-3,
) -> Circuit:
    """
    Build Grover–Rudolph state preparation for 3 qubits (q0,q1,q2), ancilla-free,
    using only the native gate set {Ry, X, CX}.

    Conventions:
      - Bitstrings are MSB→LSB: |q0 q1 q2⟩ corresponds to '000'..'111'.
      - Angle tree uses the asin / "child=1" convention:
          θ(prefix) = 2 asin( sqrt(P(prefix+1) / P(prefix)) )
        so that applying Ry(θ) on |0⟩ yields P(1|prefix) = sin²(θ/2).

    Parameters:
      - k_gain = [k0,k1,k2]: minimal per-qubit Ry gain correction (command scaling).
        Use [1,1,1] for ideal/simulation reference.
      - clip_cmd: optional cap on |theta_cmd| (hardware safety against very long pulses).
      - depth: "L0", "L01", or "full"
          * L0   : prepares q0 only
          * L01  : prepares q0 and q1|q0
          * full : prepares q0, q1|q0, and q2|(q0,q1)
      - ladder: "A" or "B" selects the UCRy CNOT-control ordering at level 2
        (see gr/ucry.py).
      - ensure_nmr_attrs: append an identity-safe tail (recommended for NMR).
    """
    if len(k_gain) != 3:
        raise ValueError("k_gain must be [k0,k1,k2] for 3 qubits.")

    ladder_u = ladder.upper()
    if ladder_u not in ("A", "B"):
        raise ValueError("ladder must be 'A' or 'B'.")

    depth_l = depth.lower()
    if depth_l not in ("l0", "l01", "full"):
        raise ValueError("depth must be one of: 'L0', 'L01', 'full'.")

    ang = angles_3q_asin_child1(prob8)

    circ = Circuit()
    q = circ.allocateQubits(3)

    # Level 0: prepare q0
    theta0_cmd = scale_clip(ang[0][()], k_gain[0], clip_cmd)
    circ << (Ry, q[0], theta0_cmd)
    if depth_l == "l0":
        if ensure_nmr_attrs:
            add_identity_safe_tail(circ, q, eps=eps_tail)
        return circ

    # Level 1: prepare q1 conditioned on q0
    th_q1_q0eq0_cmd = scale_clip(ang[1][(0,)], k_gain[1], clip_cmd)
    th_q1_q0eq1_cmd = scale_clip(ang[1][(1,)], k_gain[1], clip_cmd)
    apply_patterned_cry_on_q0(circ, q, target=1, theta_cmd=th_q1_q0eq0_cmd, q0_value=0, clip_cmd=clip_cmd)
    apply_patterned_cry_on_q0(circ, q, target=1, theta_cmd=th_q1_q0eq1_cmd, q0_value=1, clip_cmd=clip_cmd)

    if depth_l == "l01":
        if ensure_nmr_attrs:
            add_identity_safe_tail(circ, q, eps=eps_tail)
        return circ

    # Level 2: prepare q2 conditioned on (q0,q1) via UCRy (Ry+CX only)
    th00_cmd = scale_clip(ang[2][(0, 0)], k_gain[2], clip_cmd)
    th01_cmd = scale_clip(ang[2][(0, 1)], k_gain[2], clip_cmd)
    th10_cmd = scale_clip(ang[2][(1, 0)], k_gain[2], clip_cmd)
    th11_cmd = scale_clip(ang[2][(1, 1)], k_gain[2], clip_cmd)

    apply_ucry_2ctrl(
        circ,
        c0=q[0],
        c1=q[1],
        t=q[2],
        theta00_cmd=th00_cmd,
        theta01_cmd=th01_cmd,
        theta10_cmd=th10_cmd,
        theta11_cmd=th11_cmd,
        ladder=ladder_u,
        clip_cmd=clip_cmd,
    )

    if ensure_nmr_attrs:
        add_identity_safe_tail(circ, q, eps=eps_tail)
    return circ
