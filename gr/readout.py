# gr/readout.py
from __future__ import annotations

import numpy as np
from typing import Dict

from spinqit import Circuit, X

from .constants import STATES_3Q
from .backends import run_nmr_probs_robust
from .circuit import add_identity_safe_tail


def _dict_to_vec(d: Dict[str, float]) -> np.ndarray:
    return np.array([float(d.get(s, 0.0)) for s in STATES_3Q], dtype=float)


def _vec_to_dict(v: np.ndarray) -> Dict[str, float]:
    return {STATES_3Q[i]: float(v[i]) for i in range(8)}


def circ_prep_basis(bits: str, eps_tail: float = 1e-3) -> Circuit:
    c = Circuit()
    q = c.allocateQubits(3)
    if bits[0] == "1":
        c << (X, q[0])
    if bits[1] == "1":
        c << (X, q[1])
    if bits[2] == "1":
        c << (X, q[2])
    add_identity_safe_tail(c, q, eps=eps_tail)
    return c


def calibrate_readout_matrix_8x8(
    *,
    shots: int,
    eps_tail: float = 1e-3,
    base_name: str = "RO",
    max_tries: int = 6,
    base_sleep: float = 2.0,
    jitter: float = 0.35,
    cooldown_s: float = 2.0,
) -> np.ndarray:
    """
    Calibrate the full 8x8 confusion matrix M where:
      M[i,j] ~ P(meas=i | prep=j)
    with i,j in STATES_3Q order.
    """
    M = np.zeros((8, 8), dtype=float)
    for j, prep in enumerate(STATES_3Q):
        probs = run_nmr_probs_robust(
            circ_prep_basis(prep, eps_tail=eps_tail),
            name=f"{base_name}_prep_{prep}",
            shots=shots,
            max_tries=max_tries,
            base_sleep=base_sleep,
            jitter=jitter,
            cooldown_s=cooldown_s,
        )
        for i, meas in enumerate(STATES_3Q):
            M[i, j] = float(probs.get(meas, 0.0))
    return M


def mitigate_readout(probs_meas: Dict[str, float], M: np.ndarray, ridge: float = 1e-3) -> Dict[str, float]:
    """
    Ridge-regularized inversion:
      p_true ≈ (M + ridge I)^{-1} p_meas
    then project to simplex (clip negative, renormalize).
    """
    y = _dict_to_vec(probs_meas)
    A = M + ridge * np.eye(8)
    p = np.linalg.solve(A, y)
    p = np.maximum(p, 0.0)
    s = float(p.sum())
    if s > 0:
        p /= s
    return _vec_to_dict(p)
