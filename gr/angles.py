# gr/angles.py
from __future__ import annotations

import math
from typing import Dict, List, Tuple

from .constants import STATES_3Q


def normalize_prob8(prob8: List[float]) -> List[float]:
    """Normalize a length-8 nonnegative vector into probabilities."""
    if len(prob8) != 8:
        raise ValueError("prob8 must have length 8 (3 qubits).")
    q = [max(0.0, float(x)) for x in prob8]
    s = float(sum(q))
    if s <= 0.0:
        raise ValueError("prob8 must have positive sum.")
    return [x / s for x in q]


def target_from_prob8(prob8: List[float]) -> Dict[str, float]:
    """Return target distribution dict keyed by '000'..'111' (MSB->LSB)."""
    p = normalize_prob8(prob8)
    return {format(i, "03b"): p[i] for i in range(8)}


def _sum_indices_prob(p: List[float], prefix_bits: List[int]) -> float:
    """
    Sum probabilities p[i] whose MSB prefix matches prefix_bits.
    Index bits are (b0 b1 b2) with b0=MSB corresponding to q0.
    """
    total = 0.0
    for i in range(8):
        b0 = (i >> 2) & 1
        b1 = (i >> 1) & 1
        b2 = (i >> 0) & 1
        bits = [b0, b1, b2]
        if bits[: len(prefix_bits)] == prefix_bits:
            total += p[i]
    return total


def angles_3q_asin_child1(prob8: List[float]) -> Dict[int, Dict[Tuple[int, ...], float]]:
    """
    Grover–Rudolph angles for 3 qubits using the convention:

        theta(prefix) = 2*asin( sqrt( P(prefix+[1]) / P(prefix) ) )

    so that applying Ry(theta) on |0> yields:
        P(1 | prefix) = sin^2(theta/2).

    Returns a dict with:
      level 0: {(): theta0}
      level 1: {(0,): theta10, (1,): theta11}
      level 2: {(0,0): t200, (0,1): t201, (1,0): t210, (1,1): t211}
    """
    p = normalize_prob8(prob8)

    def theta(P_parent: float, P_child1: float) -> float:
        if P_parent <= 0.0:
            return 0.0
        x = min(1.0, max(0.0, P_child1 / P_parent))
        return 2.0 * math.asin(math.sqrt(x))

    ang: Dict[int, Dict[Tuple[int, ...], float]] = {}

    # Level 0: q0
    P = _sum_indices_prob(p, [])
    P1 = _sum_indices_prob(p, [1])
    ang[0] = {(): theta(P, P1)}

    # Level 1: q1 conditioned on q0
    P0 = _sum_indices_prob(p, [0])
    P01 = _sum_indices_prob(p, [0, 1])
    P1p = _sum_indices_prob(p, [1])
    P11 = _sum_indices_prob(p, [1, 1])
    ang[1] = {(0,): theta(P0, P01), (1,): theta(P1p, P11)}

    # Level 2: q2 conditioned on (q0,q1)
    def Ppref(a: int, b: int) -> float:
        return _sum_indices_prob(p, [a, b])

    def Ppref1(a: int, b: int) -> float:
        return _sum_indices_prob(p, [a, b, 1])

    ang[2] = {
        (0, 0): theta(Ppref(0, 0), Ppref1(0, 0)),
        (0, 1): theta(Ppref(0, 1), Ppref1(0, 1)),
        (1, 0): theta(Ppref(1, 0), Ppref1(1, 0)),
        (1, 1): theta(Ppref(1, 1), Ppref1(1, 1)),
    }
    return ang
