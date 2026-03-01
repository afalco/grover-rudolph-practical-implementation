# tests/test_ucry.py
from __future__ import annotations

import numpy as np

from gr import (
    ucry_coeffs_from_thetas,  # expected in gr/ucry.py (or re-exported)
    ladder_sign_matrix,       # expected in gr/ucry.py (or re-exported)
)

# A deterministic example (not symmetric)
TH00, TH01, TH10, TH11 = 0.3, 1.1, -0.7, 0.9


def test_ladder_sign_matrix_invertible():
    for ladder in ["A", "B"]:
        S = ladder_sign_matrix(ladder)
        assert S.shape == (4, 4)
        det = float(np.linalg.det(S))
        assert abs(det) > 1e-9


def test_ucry_coefficients_reconstruct_thetas():
    """
    The UCRy ladder is linear: theta_vec = S * a
    where a = [a0,a1,a2,a3]^T.
    This test verifies that the coefficients returned by
    ucry_coeffs_from_thetas solve that linear system for each ladder.
    """
    theta_vec = np.array([TH00, TH01, TH10, TH11], dtype=float)

    for ladder in ["A", "B"]:
        S = ladder_sign_matrix(ladder)
        a0, a1, a2, a3 = ucry_coeffs_from_thetas(TH00, TH01, TH10, TH11, ladder=ladder)
        a = np.array([a0, a1, a2, a3], dtype=float)

        recon = S @ a
        assert np.allclose(recon, theta_vec, atol=1e-12, rtol=0.0)
