# gr/__init__.py
from __future__ import annotations

from .constants import STATES_3Q, NQ_3Q
from .angles import normalize_prob8, target_from_prob8, angles_3q_asin_child1
from .circuit import build_gr_circuit_3q, add_identity_safe_tail
from .metrics import (
    normalize_dict,
    tv_l2_fidelity,
    per_qubit_marginals,
    print_distribution,
    print_marginals,
    print_compare,
)
from .backends import (
    run_sim_probs,
    run_nmr_probs_robust,
    run_nmr_repeated_avg,
    load_nmr_credentials_from_env,
)
from .readout import (
    calibrate_readout_matrix_8x8,
    mitigate_readout,
)

from .ucry import (
    ladder_sign_matrix,
    ucry_coeffs_from_thetas,
    apply_ucry_2ctrl,
)

__all__ = [
    "STATES_3Q",
    "NQ_3Q",
    "normalize_prob8",
    "target_from_prob8",
    "angles_3q_asin_child1",
    "build_gr_circuit_3q",
    "add_identity_safe_tail",
    "normalize_dict",
    "tv_l2_fidelity",
    "per_qubit_marginals",
    "print_distribution",
    "print_marginals",
    "print_compare",
    "run_sim_probs",
    "run_nmr_probs_robust",
    "run_nmr_repeated_avg",
    "load_nmr_credentials_from_env",
    "calibrate_readout_matrix_8x8",
    "mitigate_readout",
    "ladder_sign_matrix",
    "ucry_coeffs_from_thetas",
    "apply_ucry_2ctrl",
]
