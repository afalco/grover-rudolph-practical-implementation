# tests/test_angles.py
from __future__ import annotations

import math
import numpy as np

from gr import (
    target_from_prob8,
    angles_3q_asin_child1,   # expected in gr/angles.py (or re-exported in gr/__init__.py)
    STATES_3Q,
)

PROB8 = [1, 2, 3, 4, 4, 3, 2, 1]


def test_target_normalization():
    target = target_from_prob8(PROB8)
    s = sum(target.values())
    assert abs(s - 1.0) < 1e-12
    assert set(target.keys()) == set(STATES_3Q)
    assert all(v >= 0.0 for v in target.values())


def _sum_indices(p: list[float], prefix_bits: list[int]) -> float:
    # MSB->LSB: index bits are (b0,b1,b2) = (q0,q1,q2)
    total = 0.0
    for i in range(8):
        b0 = (i >> 2) & 1
        b1 = (i >> 1) & 1
        b2 = (i >> 0) & 1
        bits = [b0, b1, b2]
        if bits[: len(prefix_bits)] == prefix_bits:
            total += p[i]
    return total


def test_angles_range_and_consistency_child1():
    """
    For the 'child=1' convention:
      theta(prefix) = 2*asin(sqrt(P(prefix+[1])/P(prefix)))
    hence theta in [0, pi] and it should reproduce the conditional probability by:
      P(1 | prefix) = sin^2(theta/2)
    """
    target = target_from_prob8(PROB8)
    p = [target[f"{i:03b}"] for i in range(8)]
    ang = angles_3q_asin_child1(PROB8)

    # range checks
    for level in [0, 1, 2]:
        for pref, th in ang[level].items():
            assert 0.0 <= th <= math.pi + 1e-12

    # check each node's conditional
    # Level 0: prefix []
    P = _sum_indices(p, [])
    P1 = _sum_indices(p, [1])
    th = ang[0][()]
    cond = 0.0 if P <= 0 else P1 / P
    rec = math.sin(0.5 * th) ** 2
    assert abs(rec - cond) < 1e-12

    # Level 1: prefixes [0], [1]
    for b0 in [0, 1]:
        P = _sum_indices(p, [b0])
        P1 = _sum_indices(p, [b0, 1])
        th = ang[1][(b0,)]
        cond = 0.0 if P <= 0 else P1 / P
        rec = math.sin(0.5 * th) ** 2
        assert abs(rec - cond) < 1e-12

    # Level 2: prefixes [b0,b1]
    for b0 in [0, 1]:
        for b1 in [0, 1]:
            P = _sum_indices(p, [b0, b1])
            P1 = _sum_indices(p, [b0, b1, 1])
            th = ang[2][(b0, b1)]
            cond = 0.0 if P <= 0 else P1 / P
            rec = math.sin(0.5 * th) ** 2
            assert abs(rec - cond) < 1e-12
