# definitive_gr.py (ANONYMIZED + OPERATIVE for SpinQ Triangulum NMR)
# - No real credentials are committed.
# - Uses env vars for NMR connection:
#     SPINQ_IP, SPINQ_PORT, SPINQ_USER, SPINQ_PASS
# - Runs staged GR (L0/L01/FULL) on SIM and (optionally) on NMR with robust retries.
# - Uses only {Ry, X, CX} (ancilla-free), and a ladder-consistent UCRy solver (A/B).
#
# Quick start (Windows PowerShell):
#   setx SPINQ_IP "10.30.227.5"
#   setx SPINQ_PORT "55444"
#   setx SPINQ_USER "user1"
#   setx SPINQ_PASS "********"
#   python definitive_gr.py
#
# Quick start (bash):
#   export SPINQ_IP="10.30.227.5"
#   export SPINQ_PORT="55444"
#   export SPINQ_USER="user1"
#   export SPINQ_PASS="********"
#   python definitive_gr.py

import os
import math
import time
import random
from typing import Dict, List, Tuple, Optional

import numpy as np

from spinqit import Circuit, X, Ry, CX, get_compiler, NMRConfig
from spinqit import get_basic_simulator, BasicSimulatorConfig

# ===================== USER CONFIG =====================

# Target weights (8 entries for 3 qubits): order is |q0 q1 q2> = 000..111 (MSB->LSB)
PROB8 = [1, 2, 3, 4, 4, 3, 2, 1]

# Execution knobs
RUN_SIM = True
# RUN_NMR = False  -> simulation-only mode (no hardware connection, no SPINQ_* needed)
# RUN_NMR = True   -> hardware mode (runs on SpinQ Triangulum NMR; requires SPINQ_IP/PORT/USER/PASS)
RUN_NMR = False       

SHOTS_SIM = 200_000
SHOTS_NMR = 2048
REPEATS_NMR = 5
COOLDOWN_S = 2.0

# NMR connection (ANONYMIZED placeholders — env vars are required in practice)
# NOTE: The code uses environment variables exclusively, not these placeholders.
IP = "192.0.2.1"     # RFC 5737 documentation-only IP
PORT = 12345         # documentation-only port
USER = "userX"
PASS = "********"

# Robust retries for NMR calls
NMR_MAX_TRIES = 6
NMR_BASE_SLEEP = 2.0
NMR_JITTER = 0.35

# Optional: minimal Ry gain calibration (single point). Needs NMR.
DO_RY_GAIN_CAL = False
SHOTS_CAL = 2048

# Optional: hardware command clipping (start with None; try e.g. 0.95*pi)
CLIP_THETA_CMD: Optional[float] = None  # e.g. math.pi * 0.95

# Optional: readout mitigation (8×8). Needs NMR; takes extra time.
DO_READOUT_MITIGATION = False
SHOTS_RO = 4096
RIDGE = 1e-3

# Identity-safe tail (helps avoid backend IR attribute edge cases)
EPS_TAIL = 1e-3

# Choose level-2 CX ladder ordering: "A" or "B"
UCRY_LADDER = "B"

# ===================== INTERNALS =====================

STATES = [format(i, "03b") for i in range(8)]  # 000..111 (MSB->LSB)


# --------------------- Utilities: formatting & metrics ---------------------

def dict_to_vec(d: Dict[str, float]) -> np.ndarray:
    return np.array([float(d.get(s, 0.0)) for s in STATES], dtype=float)

def vec_to_dict(v: np.ndarray) -> Dict[str, float]:
    return {STATES[i]: float(v[i]) for i in range(8)}

def normalize_dict(d: Dict[str, float]) -> Dict[str, float]:
    s = float(sum(d.values()))
    if s <= 0:
        return {k: 0.0 for k in STATES}
    return {k: float(d.get(k, 0.0)) / s for k in STATES}

def print_distribution(title: str, probs: Dict[str, float], sig: int = 12) -> None:
    probs = normalize_dict(probs)
    fmt = f"{{:.{sig}g}}"
    print(f"\n{title}")
    print("State order:", STATES)
    print("Probs:", [fmt.format(probs[s]) for s in STATES])
    dom = max(probs.items(), key=lambda kv: kv[1])
    print(f"Dominant: {dom[0]} (p={dom[1]:.6f})")

def per_qubit_marginals(probs: Dict[str, float]) -> List[Tuple[float, float]]:
    probs = normalize_dict(probs)
    out = []
    for qi in range(3):
        p0 = sum(v for b, v in probs.items() if b[qi] == "0")
        p1 = sum(v for b, v in probs.items() if b[qi] == "1")
        out.append((p0, p1))
    return out

def print_marginals(title: str, probs: Dict[str, float], sig: int = 12) -> None:
    fmt = f"{{:.{sig}g}}"
    print(f"\n{title} — per-qubit marginals (MSB->LSB):")
    for i, (p0, p1) in enumerate(per_qubit_marginals(probs)):
        print(f"  q{i}: P(0)={fmt.format(p0)}  P(1)={fmt.format(p1)}")

def compare(a: Dict[str, float], b: Dict[str, float], name_a="A", name_b="B") -> Tuple[float, float, float]:
    a = normalize_dict(a)
    b = normalize_dict(b)
    p = dict_to_vec(a)
    q = dict_to_vec(b)
    tv = 0.5 * float(np.sum(np.abs(p - q)))
    l2 = float(np.sqrt(np.sum((p - q) ** 2)))
    fidelity = float((np.sum(np.sqrt(np.maximum(p, 0) * np.maximum(q, 0)))) ** 2)
    print(f"\n=== Comparison ({name_a} vs {name_b}) ===")
    print("TV:", tv)
    print("L2:", l2)
    print("Fidelity:", fidelity)
    return tv, l2, fidelity


# --------------------- Target ---------------------

def normalize_prob8(prob8: List[float]) -> List[float]:
    q = [max(0.0, float(x)) for x in prob8]
    s = float(sum(q))
    if s <= 0:
        raise ValueError("prob8 must have positive sum")
    return [x / s for x in q]

def target_from_prob8(prob8: List[float]) -> Dict[str, float]:
    p = normalize_prob8(prob8)
    return {format(i, "03b"): p[i] for i in range(8)}


# --------------------- Identity-safe tail ---------------------

def add_identity_safe_tail(circ: Circuit, q):
    circ << (Ry, q[0], EPS_TAIL)
    circ << (Ry, q[0], -EPS_TAIL)


# --------------------- Simulator execution ---------------------

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
    return normalize_dict(probs)


# --------------------- NMR execution (env vars, robust) ---------------------

def _get_nmr_backend():
    from spinqit import get_nmr
    return get_nmr()

def _nmr_env(var: str) -> Optional[str]:
    v = os.environ.get(var, None)
    if v is None:
        return None
    v = str(v).strip()
    return v if v else None

def _load_nmr_credentials() -> Tuple[str, int, str, str]:
    ip = _nmr_env("SPINQ_IP")
    port_s = _nmr_env("SPINQ_PORT")
    user = _nmr_env("SPINQ_USER")
    pw = _nmr_env("SPINQ_PASS")

    missing = [k for k, v in [("SPINQ_IP", ip), ("SPINQ_PORT", port_s), ("SPINQ_USER", user), ("SPINQ_PASS", pw)] if not v]
    if missing:
        raise RuntimeError(
            "Missing NMR configuration in environment variables: "
            + ", ".join(missing)
            + "\nSet them before running with RUN_NMR=True."
        )
    try:
        port = int(port_s)  # type: ignore[arg-type]
    except Exception as e:
        raise RuntimeError(f"SPINQ_PORT must be an integer, got '{port_s}'.") from e
    return ip, port, user, pw  # type: ignore[return-value]

def _make_nmr_config(shots: int, name: str) -> NMRConfig:
    ip, port, user, pw = _load_nmr_credentials()
    cfg = NMRConfig()
    cfg.configure_shots(shots)
    cfg.configure_ip(ip)
    cfg.configure_port(port)
    cfg.configure_account(user, pw)
    cfg.configure_task(name, name)
    return cfg

def run_nmr_probs_robust(circ: Circuit, name: str, shots: int) -> Dict[str, float]:
    comp = get_compiler("native")
    exe = comp.compile(circ, 0)

    last_err: Optional[Exception] = None
    for attempt in range(1, NMR_MAX_TRIES + 1):
        try:
            eng = _get_nmr_backend()
            cfg = _make_nmr_config(shots, name)
            res = eng.execute(exe, cfg)

            out = getattr(res, "probabilities", None) or getattr(res, "counts", None)
            if out is None or len(out) == 0:
                raise RuntimeError("NMR returned empty probabilities/counts.")

            probs = {k: float(out.get(k, 0.0)) for k in STATES}
            probs = normalize_dict(probs)

            if COOLDOWN_S > 0:
                time.sleep(COOLDOWN_S)
            return probs

        except Exception as e:
            last_err = e
            sleep_s = NMR_BASE_SLEEP * (2 ** (attempt - 1))
            sleep_s *= (1.0 + random.uniform(-NMR_JITTER, NMR_JITTER))
            sleep_s = max(0.5, sleep_s)
            print(f"[NMR] attempt {attempt}/{NMR_MAX_TRIES} failed for '{name}': {e}")
            if attempt < NMR_MAX_TRIES:
                print(f"[NMR] sleeping {sleep_s:.2f}s then retrying...")
                time.sleep(sleep_s)

    raise RuntimeError(f"NMR job '{name}' failed after {NMR_MAX_TRIES} attempts. Last error: {last_err}")

def run_nmr_repeated_avg(circ: Circuit, base_name: str, shots: int, repeats: int) -> Dict[str, float]:
    outs = []
    for r in range(repeats):
        outs.append(run_nmr_probs_robust(circ, name=f"{base_name}_r{r+1}", shots=shots))
    avg = {s: 0.0 for s in STATES}
    for d in outs:
        for s in STATES:
            avg[s] += d[s]
    for s in STATES:
        avg[s] /= max(1, len(outs))
    return normalize_dict(avg)


# --------------------- Grover–Rudolph angles (asin / child=1 convention) ---------------------

def sum_indices_prob(p: List[float], prefix_bits: List[int]) -> float:
    total = 0.0
    for i in range(8):
        b0 = (i >> 2) & 1
        b1 = (i >> 1) & 1
        b2 = (i >> 0) & 1
        bits = [b0, b1, b2]  # MSB->LSB = q0 q1 q2
        if bits[:len(prefix_bits)] == prefix_bits:
            total += p[i]
    return total

def angles_3q_asin_child1(prob8: List[float]) -> Dict[int, Dict[Tuple[int, ...], float]]:
    p = normalize_prob8(prob8)

    def theta(P_parent: float, P_child1: float) -> float:
        if P_parent <= 0:
            return 0.0
        x = min(1.0, max(0.0, P_child1 / P_parent))
        return 2.0 * math.asin(math.sqrt(x))

    ang: Dict[int, Dict[Tuple[int, ...], float]] = {}

    P = sum_indices_prob(p, [])
    P1 = sum_indices_prob(p, [1])
    ang[0] = {(): theta(P, P1)}

    P0 = sum_indices_prob(p, [0]); P01 = sum_indices_prob(p, [0, 1])
    P1p = sum_indices_prob(p, [1]); P11 = sum_indices_prob(p, [1, 1])
    ang[1] = {(0,): theta(P0, P01), (1,): theta(P1p, P11)}

    def Ppref(a, b): return sum_indices_prob(p, [a, b])
    def Ppref1(a, b): return sum_indices_prob(p, [a, b, 1])
    ang[2] = {
        (0, 0): theta(Ppref(0, 0), Ppref1(0, 0)),
        (0, 1): theta(Ppref(0, 1), Ppref1(0, 1)),
        (1, 0): theta(Ppref(1, 0), Ppref1(1, 0)),
        (1, 1): theta(Ppref(1, 1), Ppref1(1, 1)),
    }
    return ang


# --------------------- Gate building blocks (Ry + CX only) ---------------------

def scale_clip(theta: float, k: float, clip_cmd: Optional[float]) -> float:
    theta_cmd = theta / (k if k > 1e-6 else 1.0)
    if clip_cmd is not None and abs(theta_cmd) > clip_cmd:
        theta_cmd = math.copysign(clip_cmd, theta_cmd)
    return theta_cmd

def apply_cry_1ctrl_decomp(circ: Circuit, c, t, theta_cmd: float, clip_cmd: Optional[float]) -> None:
    half = 0.5 * theta_cmd
    if clip_cmd is not None and abs(half) > clip_cmd:
        half = math.copysign(clip_cmd, half)
    circ << (Ry, t, half)
    circ << (CX, (c, t))
    circ << (Ry, t, -half)
    circ << (CX, (c, t))

def apply_patterned_cry_on_q0(circ: Circuit, q, target: int, theta_cmd: float, q0_value: int, clip_cmd: Optional[float]) -> None:
    if q0_value == 0:
        circ << (X, q[0])
    apply_cry_1ctrl_decomp(circ, c=q[0], t=q[target], theta_cmd=theta_cmd, clip_cmd=clip_cmd)
    if q0_value == 0:
        circ << (X, q[0])


# --------------------- Level-2 UCRy coefficients via sign-system ---------------------

def ladder_sign_matrix(ladder: str) -> np.ndarray:
    ladder = ladder.upper()
    if ladder not in ("A", "B"):
        raise ValueError("ladder must be 'A' or 'B'")
    cx_controls = [1, 0, 1, 0] if ladder == "A" else [0, 1, 0, 1]
    rows = []
    for c0, c1 in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        bits = [c0, c1]
        parity = 0
        signs = []
        signs.append(+1 if parity == 0 else -1)
        for k in range(3):
            ctrl = cx_controls[k]
            if bits[ctrl] == 1:
                parity ^= 1
            signs.append(+1 if parity == 0 else -1)
        rows.append(signs)
    return np.array(rows, dtype=float)

def ucry_coeffs_from_thetas(theta00: float, theta01: float, theta10: float, theta11: float, ladder: str) -> Tuple[float, float, float, float]:
    S = ladder_sign_matrix(ladder)
    theta_vec = np.array([theta00, theta01, theta10, theta11], dtype=float)
    a = np.linalg.solve(S, theta_vec)
    return float(a[0]), float(a[1]), float(a[2]), float(a[3])

def apply_ucry_2ctrl(circ: Circuit, c0, c1, t,
                     theta00_cmd: float, theta01_cmd: float, theta10_cmd: float, theta11_cmd: float,
                     ladder: str, clip_cmd: Optional[float]) -> None:
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


# --------------------- Build Grover–Rudolph circuit ---------------------

def build_gr_circuit_3q(prob8: List[float],
                        k_gain: List[float],
                        clip_cmd: Optional[float],
                        depth: str,
                        ladder: str,
                        ensure_nmr_attrs: bool) -> Circuit:
    depth = depth.lower()
    if depth not in ("l0", "l01", "full"):
        raise ValueError("depth must be one of: 'L0','L01','full'")

    ang = angles_3q_asin_child1(prob8)
    circ = Circuit()
    q = circ.allocateQubits(3)

    theta0_cmd = scale_clip(ang[0][()], k_gain[0], clip_cmd)
    circ << (Ry, q[0], theta0_cmd)
    if depth == "l0":
        if ensure_nmr_attrs:
            add_identity_safe_tail(circ, q)
        return circ

    th0 = scale_clip(ang[1][(0,)], k_gain[1], clip_cmd)
    th1 = scale_clip(ang[1][(1,)], k_gain[1], clip_cmd)
    apply_patterned_cry_on_q0(circ, q, target=1, theta_cmd=th0, q0_value=0, clip_cmd=clip_cmd)
    apply_patterned_cry_on_q0(circ, q, target=1, theta_cmd=th1, q0_value=1, clip_cmd=clip_cmd)

    if depth == "l01":
        if ensure_nmr_attrs:
            add_identity_safe_tail(circ, q)
        return circ

    th00 = scale_clip(ang[2][(0, 0)], k_gain[2], clip_cmd)
    th01 = scale_clip(ang[2][(0, 1)], k_gain[2], clip_cmd)
    th10 = scale_clip(ang[2][(1, 0)], k_gain[2], clip_cmd)
    th11 = scale_clip(ang[2][(1, 1)], k_gain[2], clip_cmd)
    apply_ucry_2ctrl(circ, q[0], q[1], q[2], th00, th01, th10, th11, ladder=ladder, clip_cmd=clip_cmd)

    if ensure_nmr_attrs:
        add_identity_safe_tail(circ, q)
    return circ


# --------------------- Optional: gain calibration ---------------------

def calibrate_ry_gain_single_point(qi: int, shots: int) -> float:
    c = Circuit()
    q = c.allocateQubits(3)
    c << (Ry, q[qi], math.pi / 2)
    add_identity_safe_tail(c, q)
    probs = run_nmr_probs_robust(c, name=f"CAL_RY_q{qi}", shots=shots)
    p1 = sum(v for b, v in probs.items() if b[qi] == "1")
    p1 = max(0.0, min(1.0, p1))
    k = (4.0 / math.pi) * math.asin(math.sqrt(p1))
    return max(1e-6, float(k))


# --------------------- Optional: readout mitigation ---------------------

def circ_prep_basis(bits: str) -> Circuit:
    c = Circuit()
    q = c.allocateQubits(3)
    if bits[0] == "1": c << (X, q[0])
    if bits[1] == "1": c << (X, q[1])
    if bits[2] == "1": c << (X, q[2])
    add_identity_safe_tail(c, q)
    return c

def calibrate_readout_matrix_8x8(shots: int) -> np.ndarray:
    M = np.zeros((8, 8), dtype=float)
    for j, prep in enumerate(STATES):
        probs = run_nmr_probs_robust(circ_prep_basis(prep), name=f"RO_prep_{prep}", shots=shots)
        for i, meas in enumerate(STATES):
            M[i, j] = probs.get(meas, 0.0)
    return M

def mitigate_readout(probs_meas: Dict[str, float], M: np.ndarray, ridge: float) -> Dict[str, float]:
    y = dict_to_vec(normalize_dict(probs_meas))
    A = M + ridge * np.eye(8)
    p = np.linalg.solve(A, y)
    p = np.maximum(p, 0.0)
    s = float(p.sum())
    if s > 0:
        p /= s
    return vec_to_dict(p)


# ===================== MAIN =====================

def run_stage(depth: str, k_gain: List[float], ladder: str, target_ref: Dict[str, float], Mfull: Optional[np.ndarray]) -> None:
    circ_sim = build_gr_circuit_3q(PROB8, k_gain=[1.0, 1.0, 1.0], clip_cmd=None,
                                   depth=depth, ladder=ladder, ensure_nmr_attrs=False)
    circ_nmr = build_gr_circuit_3q(PROB8, k_gain=k_gain, clip_cmd=CLIP_THETA_CMD,
                                   depth=depth, ladder=ladder, ensure_nmr_attrs=True)

    print("\n" + "=" * 90)
    print(f"STAGE {depth.upper()}  (ladder {ladder.upper()})")
    print("=" * 90)

    sim_probs = None
    if RUN_SIM:
        sim_probs = run_sim_probs(circ_sim, shots=SHOTS_SIM)
        print_distribution(f"SIM ideal — {depth.upper()}", sim_probs)
        print_marginals(f"SIM — {depth.upper()}", sim_probs)

    if RUN_NMR:
        nmr_avg = run_nmr_repeated_avg(circ_nmr, base_name=f"GR_{depth.upper()}_{ladder}", shots=SHOTS_NMR, repeats=REPEATS_NMR)
        print_distribution(f"NMR raw avg — {depth.upper()}", nmr_avg)
        print_marginals(f"NMR raw — {depth.upper()}", nmr_avg)

        if sim_probs is not None:
            compare(sim_probs, nmr_avg, name_a=f"SIM {depth.upper()}", name_b=f"NMR {depth.upper()} raw avg")

        if depth.lower() == "full":
            compare(target_ref, nmr_avg, name_a="Target", name_b="NMR FULL raw avg")

        if Mfull is not None:
            mitig = mitigate_readout(nmr_avg, Mfull, ridge=RIDGE)
            print_distribution(f"NMR mitigated — {depth.upper()} (ridge={RIDGE})", mitig)
            print_marginals(f"NMR mitigated — {depth.upper()}", mitig)
            if sim_probs is not None:
                compare(sim_probs, mitig, name_a=f"SIM {depth.upper()}", name_b=f"NMR {depth.upper()} mitigated")
            if depth.lower() == "full":
                compare(target_ref, mitig, name_a="Target", name_b="NMR FULL mitigated")


def main():
    target = target_from_prob8(PROB8)
    print_distribution("Target distribution (normalized)", target)
    print_marginals("Target", target)

    k_gain = [1.0, 1.0, 1.0]
    if RUN_NMR and DO_RY_GAIN_CAL:
        print("\n=== Calibrating Ry gains (single-point) ===")
        k_gain = []
        for qi in range(3):
            k = calibrate_ry_gain_single_point(qi, shots=SHOTS_CAL)
            k_gain.append(k)
            print(f"q{qi}: k = {k:.6g}")
        print("Calibrated gains:", k_gain)

    Mfull = None
    if RUN_NMR and DO_READOUT_MITIGATION:
        print("\n=== Calibrating 8x8 readout matrix (Mfull) ===")
        Mfull = calibrate_readout_matrix_8x8(shots=SHOTS_RO)
        print("Done. Condition number:", float(np.linalg.cond(Mfull)))

    for depth in ["L0", "L01", "full"]:
        run_stage(depth, k_gain=k_gain, ladder=UCRY_LADDER, target_ref=target, Mfull=Mfull)

    if RUN_SIM:
        circ_check = build_gr_circuit_3q(PROB8, k_gain=[1.0, 1.0, 1.0], clip_cmd=None,
                                         depth="full", ladder=UCRY_LADDER, ensure_nmr_attrs=False)
        sim_full = run_sim_probs(circ_check, shots=SHOTS_SIM)
        print("\n" + "=" * 90)
        print("SANITY CHECK (SIM FULL vs Target)")
        print("=" * 90)
        compare(target, sim_full, name_a="Target", name_b="SIM FULL")


if __name__ == "__main__":
    main()
