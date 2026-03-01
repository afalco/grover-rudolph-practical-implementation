# tests/test_metrics.py
from __future__ import annotations

import numpy as np

from gr import (
    tv_l2_fidelity,   # expected in gr/metrics.py (or re-exported)
    STATES_3Q,
)


def _dist(vals):
    return {STATES_3Q[i]: float(vals[i]) for i in range(8)}


def test_metrics_identity_case():
    p = _dist([0.1, 0.2, 0.0, 0.1, 0.0, 0.3, 0.2, 0.1])
    tv, l2, fid = tv_l2_fidelity(p, p)
    assert abs(tv) < 1e-15
    assert abs(l2) < 1e-15
    assert abs(fid - 1.0) < 1e-15


def test_metrics_basic_properties():
    p = _dist([1, 0, 0, 0, 0, 0, 0, 0])
    q = _dist([0, 1, 0, 0, 0, 0, 0, 0])
    tv, l2, fid = tv_l2_fidelity(p, q)

    # TV distance between two disjoint delta distributions is 1
    assert abs(tv - 1.0) < 1e-15

    # L2 distance is sqrt(2)
    assert abs(l2 - np.sqrt(2)) < 1e-15

    # Classical fidelity is 0 when supports are disjoint
    assert abs(fid - 0.0) < 1e-15


def test_metrics_symmetry():
    p = _dist([0.05, 0.1, 0.15, 0.2, 0.2, 0.15, 0.1, 0.05])
    q = _dist([0.125] * 8)
    tv1, l21, f1 = tv_l2_fidelity(p, q)
    tv2, l22, f2 = tv_l2_fidelity(q, p)
    assert abs(tv1 - tv2) < 1e-15
    assert abs(l21 - l22) < 1e-15
    assert abs(f1 - f2) < 1e-15
