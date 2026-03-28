import os
import math
import time
import random
from typing import Dict, List, Optional, Tuple

import numpy as np
from spinqit import Circuit, X, Ry, CX, get_compiler, NMRConfig
from spinqit import get_basic_simulator, BasicSimulatorConfig

# ============================================================
# Minimal smoke-test for Triangulum using the same conventions
# as definitive_gr.py, but restricted to one FULL run.
# ============================================================

# --- Test vector (right-skew example; replace if needed) ---
PROB8 = [1, 3, 5, 7, 9, 11, 13, 15]

# --- Execution knobs ---
RUN_SIM = True
RUN_NMR = True
SHOTS_SIM = 50000
SHOTS_NMR = 2048
COOLDOWN_S = 2.0

# --- Same hardware-style knobs as the experimental code ---
CLIP_THETA_CMD: Optional[float] = None
EPS_TAIL = 1e-3
UCRY_LADDER = "B"       # "A" or "B"
STAGE = "full"          # "l0", "l01", "full"
TASK_NAME = "GR_test_full_B_n1"

# --- Optional bit-order canonicalization for comparison ---
SIMULATOR_ORDER = "MSB"  # simulator typically returns 000..111 in MSB->LSB convention
HARDWARE_ORDER = "MSB"   # Triangulum export may effectively be LSB->MSB

# --- Robust retries for NMR ---
NMR_MAX_TRIES = 4
NMR_BASE_SLEEP = 2.0
NMR_JITTER = 0.35

STATES = [format(i, "03b") for i in range(8)]


def normalize_prob8(prob8: List[float]) -> List[float]:
    q = [max(0.0, float(x)) for x in prob8]
    s = float(sum(q))
    if s <= 0:
        raise ValueError("prob8 must have positive sum")
    return [x / s for x in q]


def target_from_prob8(prob8: List[float]) -> Dict[str, float]:
    p = normalize_prob8(prob8)
    return {format(i, "03b"): p[i] for i in range(8)}


def normalize_dict(d: Dict[str, float]) -> Dict[str, float]:
    s = float(sum(d.values()))
    if s <= 0:
        return {k: 0.0 for k in STATES}
    return {k: float(d.get(k, 0.0)) / s for k in STATES}


def reverse_bitstring_keys(d: Dict[str, float]) -> Dict[str, float]:
    return {k[::-1]: float(v) for k, v in d.items()}


def canonicalize_probs(d: Dict[str, float], source_order: str) -> Dict[str, float]:
    source_order = source_order.upper().strip()
    if source_order not in ("MSB", "LSB"):
        raise ValueError("source_order must be 'MSB' or 'LSB'")
    d = {str(k): float(v) for k, v in d.items() if str(k) in STATES}
    if source_order == "LSB":
        d = reverse_bitstring_keys(d)
    return normalize_dict({s: d.get(s, 0.0) for s in STATES})


def print_distribution(title: str, probs: Dict[str, float], sig: int = 10) -> None:
    probs = normalize_dict(probs)
    fmt = f"{{:.{sig}g}}"
    print(f"\n{title}")
    print("State order:", STATES)
    print("Probs:", [fmt.format(probs[s]) for s in STATES])


def compare(a: Dict[str, float], b: Dict[str, float], name_a: str = "A", name_b: str = "B") -> Tuple[float, float, float]:
    a = normalize_dict(a)
    b = normalize_dict(b)
    p = np.array([a[s] for s in STATES], dtype=float)
    q = np.array([b[s] for s in STATES], dtype=float)
    tv = 0.5 * float(np.sum(np.abs(p - q)))
    l2 = float(np.sqrt(np.sum((p - q) ** 2)))
    fidelity = float((np.sum(np.sqrt(np.maximum(p, 0) * np.maximum(q, 0)))) ** 2)
    print(f"\n=== Comparison ({name_a} vs {name_b}) ===")
    print("TV:", tv)
    print("L2:", l2)
    print("Fidelity:", fidelity)
    return tv, l2, fidelity


# ---------------- Angle synthesis ----------------

def angles_3q_asin_child1(prob8: List[float]):
    p = normalize_prob8(prob8)
    arr = np.array(p, dtype=float)
    a = arr.reshape(2, 2, 2)

    def theta_from_pair(x0: float, x1: float) -> float:
        s = x0 + x1
        if s <= 0:
            return 0.0
        r = max(0.0, min(1.0, x1 / s))
        return 2.0 * math.asin(math.sqrt(r))

    ang0 = {(): theta_from_pair(float(np.sum(a[0, :, :])), float(np.sum(a[1, :, :])))}
    ang1 = {
        (0,): theta_from_pair(float(np.sum(a[0, 0, :])), float(np.sum(a[0, 1, :]))),
        (1,): theta_from_pair(float(np.sum(a[1, 0, :])), float(np.sum(a[1, 1, :]))),
    }
    ang2 = {
        (0, 0): theta_from_pair(float(a[0, 0, 0]), float(a[0, 0, 1])),
        (0, 1): theta_from_pair(float(a[0, 1, 0]), float(a[0, 1, 1])),
        (1, 0): theta_from_pair(float(a[1, 0, 0]), float(a[1, 0, 1])),
        (1, 1): theta_from_pair(float(a[1, 1, 0]), float(a[1, 1, 1])),
    }
    return {0: ang0, 1: ang1, 2: ang2}


def scale_clip(theta: float, gain: float, clip_cmd: Optional[float]) -> float:
    x = gain * theta
    if clip_cmd is not None and abs(x) > clip_cmd:
        x = math.copysign(clip_cmd, x)
    return x


def add_identity_safe_tail(circ: Circuit, q):
    circ << (Ry, q[0], EPS_TAIL)
    circ << (Ry, q[0], -EPS_TAIL)


def apply_patterned_cry_on_q0(circ: Circuit, q, target: int, theta_cmd: float, q0_value: int, clip_cmd: Optional[float]):
    t = theta_cmd
    if clip_cmd is not None and abs(t) > clip_cmd:
        t = math.copysign(clip_cmd, t)
    if q0_value == 0:
        circ << (X, q[0])
    circ << (CX, (q[0], q[target]))
    circ << (Ry, q[target], t / 2.0)
    circ << (CX, (q[0], q[target]))
    circ << (Ry, q[target], -t / 2.0)
    if q0_value == 0:
        circ << (X, q[0])


def ladder_sign_matrix(ladder: str) -> np.ndarray:
    ladder = ladder.upper()
    if ladder not in ("A", "B"):
        raise ValueError("ladder must be 'A' or 'B'")
    cx_controls = [1, 0, 1, 0] if ladder == "A" else [0, 1, 0, 1]
    rows = []
    for c0, c1 in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        bits = [c0, c1]
        parity = 0
        signs = [+1]
        for k in range(3):
            ctrl = cx_controls[k]
            if bits[ctrl] == 1:
                parity ^= 1
            signs.append(+1 if parity == 0 else -1)
        rows.append(signs)
    return np.array(rows, dtype=float)


def ucry_coeffs_from_thetas(theta00: float, theta01: float, theta10: float, theta11: float, ladder: str):
    S = ladder_sign_matrix(ladder)
    theta_vec = np.array([theta00, theta01, theta10, theta11], dtype=float)
    a = np.linalg.solve(S, theta_vec)
    return float(a[0]), float(a[1]), float(a[2]), float(a[3])


def apply_ucry_2ctrl(circ: Circuit, c0, c1, t, theta00_cmd: float, theta01_cmd: float, theta10_cmd: float, theta11_cmd: float, ladder: str, clip_cmd: Optional[float]):
    a0, a1, a2, a3 = ucry_coeffs_from_thetas(theta00_cmd, theta01_cmd, theta10_cmd, theta11_cmd, ladder=ladder)

    def clip(x: float) -> float:
        if clip_cmd is None:
            return x
        if abs(x) > clip_cmd:
            return math.copysign(clip_cmd, x)
        return x

    ladder = ladder.upper()
    circ << (Ry, t, clip(a0))
    if ladder == "A":
        circ << (CX, (c1, t))
        circ << (Ry, t, clip(a1))
        circ << (CX, (c0, t))
        circ << (Ry, t, clip(a2))
        circ << (CX, (c1, t))
        circ << (Ry, t, clip(a3))
        circ << (CX, (c0, t))
    else:
        circ << (CX, (c0, t))
        circ << (Ry, t, clip(a1))
        circ << (CX, (c1, t))
        circ << (Ry, t, clip(a2))
        circ << (CX, (c0, t))
        circ << (Ry, t, clip(a3))
        circ << (CX, (c1, t))


def build_gr_circuit_3q(prob8: List[float], depth: str, ladder: str, ensure_nmr_attrs: bool) -> Circuit:
    depth = depth.lower()
    if depth not in ("l0", "l01", "full"):
        raise ValueError("depth must be one of: 'l0', 'l01', 'full'")

    ang = angles_3q_asin_child1(prob8)
    circ = Circuit()
    q = circ.allocateQubits(3)

    circ << (Ry, q[0], scale_clip(ang[0][()], 1.0, CLIP_THETA_CMD))
    if depth == "l0":
        if ensure_nmr_attrs:
            add_identity_safe_tail(circ, q)
        return circ

    apply_patterned_cry_on_q0(circ, q, target=1, theta_cmd=scale_clip(ang[1][(0,)], 1.0, CLIP_THETA_CMD), q0_value=0, clip_cmd=CLIP_THETA_CMD)
    apply_patterned_cry_on_q0(circ, q, target=1, theta_cmd=scale_clip(ang[1][(1,)], 1.0, CLIP_THETA_CMD), q0_value=1, clip_cmd=CLIP_THETA_CMD)

    if depth == "l01":
        if ensure_nmr_attrs:
            add_identity_safe_tail(circ, q)
        return circ

    apply_ucry_2ctrl(
        circ, q[0], q[1], q[2],
        scale_clip(ang[2][(0, 0)], 1.0, CLIP_THETA_CMD),
        scale_clip(ang[2][(0, 1)], 1.0, CLIP_THETA_CMD),
        scale_clip(ang[2][(1, 0)], 1.0, CLIP_THETA_CMD),
        scale_clip(ang[2][(1, 1)], 1.0, CLIP_THETA_CMD),
        ladder=ladder,
        clip_cmd=CLIP_THETA_CMD,
    )

    if ensure_nmr_attrs:
        add_identity_safe_tail(circ, q)
    return circ


# ---------------- Execution ----------------

def run_sim_probs(circ: Circuit, shots: int) -> Dict[str, float]:
    comp = get_compiler("native")
    exe = comp.compile(circ, 0)
    sim = get_basic_simulator()
    cfg = BasicSimulatorConfig()
    cfg.configure_shots(shots)
    res = sim.execute(exe, cfg)
    out = getattr(res, "probabilities", None) or getattr(res, "counts", None)
    if out is None or len(out) == 0:
        raise RuntimeError("Simulator returned empty probabilities/counts.")
    probs = {k: float(out.get(k, 0.0)) for k in STATES}
    return canonicalize_probs(probs, SIMULATOR_ORDER)


def _get_nmr_backend():
    from spinqit import get_nmr
    return get_nmr()


def _env(name: str) -> str:
    v = os.environ.get(name, "").strip()
    if not v:
        raise RuntimeError(f"Missing environment variable: {name}")
    return v


def _make_nmr_config(shots: int, task_name: str) -> NMRConfig:
    cfg = NMRConfig()
    cfg.configure_shots(shots)
    cfg.configure_ip(_env("SPINQ_IP"))
    cfg.configure_port(int(_env("SPINQ_PORT")))
    cfg.configure_account(_env("SPINQ_USER"), _env("SPINQ_PASS"))
    cfg.configure_task(task_name, task_name)
    return cfg


def run_nmr_probs_robust(circ: Circuit, task_name: str, shots: int) -> Dict[str, float]:
    comp = get_compiler("native")
    exe = comp.compile(circ, 0)
    last_err = None
    for attempt in range(1, NMR_MAX_TRIES + 1):
        try:
            eng = _get_nmr_backend()
            cfg = _make_nmr_config(shots, task_name)
            res = eng.execute(exe, cfg)
            out = getattr(res, "probabilities", None) or getattr(res, "counts", None)
            if out is None or len(out) == 0:
                raise RuntimeError("NMR returned empty probabilities/counts.")
            probs = {k: float(out.get(k, 0.0)) for k in STATES}
            probs = canonicalize_probs(probs, HARDWARE_ORDER)
            if COOLDOWN_S > 0:
                time.sleep(COOLDOWN_S)
            return probs
        except Exception as e:
            last_err = e
            sleep_s = NMR_BASE_SLEEP * (2 ** (attempt - 1))
            sleep_s *= (1.0 + random.uniform(-NMR_JITTER, NMR_JITTER))
            sleep_s = max(0.5, sleep_s)
            print(f"[NMR] attempt {attempt}/{NMR_MAX_TRIES} failed for '{task_name}': {e}")
            if attempt < NMR_MAX_TRIES:
                print(f"[NMR] sleeping {sleep_s:.2f}s then retrying...")
                time.sleep(sleep_s)
    raise RuntimeError(f"NMR job '{task_name}' failed after {NMR_MAX_TRIES} attempts. Last error: {last_err}")


if __name__ == "__main__":
    target = target_from_prob8(PROB8)
    print_distribution("Target", target)

    circ_sim = build_gr_circuit_3q(PROB8, depth=STAGE, ladder=UCRY_LADDER, ensure_nmr_attrs=False)
    circ_nmr = build_gr_circuit_3q(PROB8, depth=STAGE, ladder=UCRY_LADDER, ensure_nmr_attrs=True)

    if RUN_SIM:
        sim_probs = run_sim_probs(circ_sim, SHOTS_SIM)
        print_distribution("SIM", sim_probs)
        compare(target, sim_probs, "Target", "SIM")
    else:
        sim_probs = None

    if RUN_NMR:
        nmr_probs = run_nmr_probs_robust(circ_nmr, task_name=TASK_NAME, shots=SHOTS_NMR)
        print_distribution("NMR raw", nmr_probs)
        compare(target, nmr_probs, "Target", "NMR")
        if sim_probs is not None:
            compare(sim_probs, nmr_probs, "SIM", "NMR")
