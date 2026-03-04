from __future__ import annotations

from typing import Dict, List, Tuple
import numpy as np

from .constants import STATES_3Q


def normalize_dict(d: Dict[str, float]) -> Dict[str, float]:
    s = float(sum(d.values()))
    if s <= 0.0:
        return {k: 0.0 for k in STATES_3Q}
    return {k: float(d.get(k, 0.0)) / s for k in STATES_3Q}


def dict_to_vec(d: Dict[str, float]) -> np.ndarray:
    d = normalize_dict(d)
    return np.array([float(d.get(k, 0.0)) for k in STATES_3Q], dtype=float)


def tv_l2_fidelity(a: Dict[str, float], b: Dict[str, float]) -> Tuple[float, float, float]:
    """
    Total variation distance: 0.5 * ||p - q||_1
    L2 distance: ||p - q||_2
    Classical fidelity: (sum_i sqrt(p_i q_i))^2
    """
    p = dict_to_vec(a)
    q = dict_to_vec(b)
    tv = 0.5 * float(np.sum(np.abs(p - q)))
    l2 = float(np.sqrt(np.sum((p - q) ** 2)))
    fid = float((np.sum(np.sqrt(np.maximum(p, 0.0) * np.maximum(q, 0.0)))) ** 2)
    return tv, l2, fid


def per_qubit_marginals(probs: Dict[str, float]) -> List[Tuple[float, float]]:
    probs = normalize_dict(probs)
    out: List[Tuple[float, float]] = []
    for qi in range(3):  # MSB->LSB
        p0 = sum(v for bs, v in probs.items() if bs[qi] == "0")
        p1 = sum(v for bs, v in probs.items() if bs[qi] == "1")
        out.append((float(p0), float(p1)))
    return out


def print_distribution(title: str, probs: Dict[str, float], sig: int = 12) -> None:
    probs = normalize_dict(probs)
    fmt = f"{{:.{sig}g}}"
    print(f"\n{title}")
    print("State order:", STATES_3Q)
    print("Probs:", [fmt.format(probs[s]) for s in STATES_3Q])
    dom = max(probs.items(), key=lambda kv: kv[1])
    print(f"Dominant: {dom[0]} (p={dom[1]:.6f})")


def print_marginals(title: str, probs: Dict[str, float], sig: int = 12) -> None:
    fmt = f"{{:.{sig}g}}"
    print(f"\n{title} — per-qubit marginals (MSB->LSB):")
    for i, (p0, p1) in enumerate(per_qubit_marginals(probs)):
        print(f"  q{i}: P(0)={fmt.format(p0)}  P(1)={fmt.format(p1)}")


def print_compare(a: Dict[str, float], b: Dict[str, float], name_a: str = "A", name_b: str = "B") -> None:
    tv, l2, fid = tv_l2_fidelity(a, b)
    print(f"\n=== Comparison ({name_a} vs {name_b}) ===")
    print("TV:", tv)
    print("L2:", l2)
    print("Fidelity:", fid)
