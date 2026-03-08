# gr/backends.py
from __future__ import annotations

import os
import random
import time
from typing import Dict, Optional, Tuple, List

from spinqit import NMRConfig, get_compiler
from spinqit import get_basic_simulator, BasicSimulatorConfig

from .constants import STATES_3Q


# ======================= Normalization =======================

def _normalize_out(out: Dict[str, float]) -> Dict[str, float]:
    s = float(sum(out.values()))
    if s <= 0.0:
        return {k: 0.0 for k in STATES_3Q}
    return {k: float(out.get(k, 0.0)) / s for k in STATES_3Q}


# ======================= Env helpers =======================

def _env(var: str) -> Optional[str]:
    v = os.environ.get(var, None)
    if v is None:
        return None
    v = str(v).strip()
    return v if v else None


# ======================= Bitorder remapping (NMR only) =======================
#
# Goal: map backend-reported bitstrings to the repo canonical order:
#   canonical: MSB->LSB = |q0 q1 q2> with keys in STATES_3Q = ["000","001",...,"111"].
#
# Controlled by env var:
#   SPINQ_BITORDER
#
# Supported values:
#   - "MSB->LSB", "msb", "msb_to_lsb", "012"  -> identity
#   - "LSB->MSB", "lsb", "lsb_to_msb"         -> reverse (bs[::-1])
#   - permutation like "201" meaning:
#         new[0] = old[2], new[1] = old[0], new[2] = old[1]
#

def _parse_bitorder_spec(spec_raw: Optional[str]) -> Optional[str]:
    if spec_raw is None:
        return None
    s = spec_raw.strip()
    if not s:
        return None

    low = s.lower()
    if low in ("msb", "msb->lsb", "msb_to_lsb", "msb2lsb"):
        return "012"  # identity
    if low in ("lsb", "lsb->msb", "lsb_to_msb", "lsb2msb"):
        return "LSB->MSB"
    if s == "MSB->LSB":
        return "012"
    if s == "LSB->MSB":
        return "LSB->MSB"

    # "012" etc: explicit permutation
    if len(s) == 3 and all(c in "012" for c in s):
        # validate it's a permutation
        if sorted(s) != ["0", "1", "2"]:
            raise ValueError(
                f"Invalid SPINQ_BITORDER='{spec_raw}': must be a permutation of 0,1,2 (e.g. '012','201')."
            )
        return s

    raise ValueError(
        f"Invalid SPINQ_BITORDER='{spec_raw}'. "
        "Use 'MSB->LSB'/'012', 'LSB->MSB', or a permutation like '201'."
    )


def _remap_bitstring(bs: str, spec: str) -> str:
    """
    Remap a 3-bit string according to spec.

    spec:
      - "012" -> identity
      - "LSB->MSB" -> reverse
      - "201" -> new[i] = old[int(spec[i])]
    """
    if len(bs) != 3:
        return bs

    if spec == "012":
        return bs
    if spec == "LSB->MSB":
        return bs[::-1]

    # permutation case
    idx = [int(c) for c in spec]  # e.g. "201" -> [2,0,1]
    return "".join(bs[i] for i in idx)


def _apply_bitorder_map(raw_out: Dict[str, float]) -> Dict[str, float]:
    """
    Map backend bitstrings to repo canonical order (MSB->LSB = |q0 q1 q2>).
    If SPINQ_BITORDER is not set, return raw_out unchanged.
    """
    spec_raw = _env("SPINQ_BITORDER")
    spec = _parse_bitorder_spec(spec_raw)
    if spec is None:
        return raw_out

    mapped: Dict[str, float] = {k: 0.0 for k in STATES_3Q}

    for bs, v in raw_out.items():
        if not isinstance(bs, str):
            continue
        if len(bs) != 3 or any(c not in "01" for c in bs):
            # ignore unexpected keys
            continue
        new_bs = _remap_bitstring(bs, spec)
        mapped[new_bs] = mapped.get(new_bs, 0.0) + float(v)

    return mapped


# ======================= SIM =======================

def run_sim_probs(circ, shots: int) -> Dict[str, float]:
    """Run on SpinQit basic simulator and return probabilities dict over STATES_3Q."""
    comp = get_compiler("native")
    exe = comp.compile(circ, 0)
    sim = get_basic_simulator()
    cfg = BasicSimulatorConfig()
    cfg.configure_shots(shots)
    res = sim.execute(exe, cfg)
    out = getattr(res, "probabilities", None) or getattr(res, "counts", None)
    if out is None or len(out) == 0:
        raise RuntimeError("Simulator returned empty probabilities/counts.")
    probs = {k: float(out.get(k, 0.0)) for k in STATES_3Q}
    return _normalize_out(probs)


# ======================= NMR =======================

def _get_nmr_backend():
    from spinqit import get_nmr
    return get_nmr()


def load_nmr_credentials_from_env() -> Tuple[str, int, str, str]:
    """
    Read NMR connection from env vars:
      SPINQ_IP, SPINQ_PORT, SPINQ_USER, SPINQ_PASS
    """
    ip = _env("SPINQ_IP")
    port_s = _env("SPINQ_PORT")
    user = _env("SPINQ_USER")
    pw = _env("SPINQ_PASS")

    missing = [k for k, v in [
        ("SPINQ_IP", ip),
        ("SPINQ_PORT", port_s),
        ("SPINQ_USER", user),
        ("SPINQ_PASS", pw),
    ] if not v]
    if missing:
        raise RuntimeError("Missing NMR env vars: " + ", ".join(missing))

    try:
        port = int(port_s)  # type: ignore[arg-type]
    except Exception as e:
        raise RuntimeError(f"SPINQ_PORT must be int, got '{port_s}'.") from e

    return ip, port, user, pw  # type: ignore[return-value]


def make_nmr_config(shots: int, name: str) -> NMRConfig:
    ip, port, user, pw = load_nmr_credentials_from_env()
    cfg = NMRConfig()
    cfg.configure_shots(shots)
    cfg.configure_ip(ip)
    cfg.configure_port(port)
    cfg.configure_account(user, pw)
    cfg.configure_task(name, name)
    return cfg


def run_nmr_probs_robust(
    circ,
    name: str,
    shots: int,
    *,
    max_tries: int = 6,
    base_sleep: float = 2.0,
    jitter: float = 0.35,
    cooldown_s: float = 2.0,
) -> Dict[str, float]:
    """
    Execute a circuit on NMR with retries + exponential backoff + jitter.
    Returns normalized probabilities dict over STATES_3Q (after optional bitorder remap).
    """
    comp = get_compiler("native")
    exe = comp.compile(circ, 0)

    last_err: Optional[Exception] = None
    for attempt in range(1, max_tries + 1):
        try:
            eng = _get_nmr_backend()
            cfg = make_nmr_config(shots, name)
            res = eng.execute(exe, cfg)

            out = getattr(res, "probabilities", None) or getattr(res, "counts", None)
            if out is None or len(out) == 0:
                raise RuntimeError("NMR returned empty probabilities/counts.")

            # Keep backend keys first; then remap (if enabled) into canonical STATES_3Q
            raw_backend = {str(k): float(out.get(k, 0.0)) for k in out.keys()}
            raw_mapped = _apply_bitorder_map(raw_backend)

            # If no bitorder mapping was applied, raw_mapped might not be restricted to STATES_3Q.
            # Normalize on canonical support.
            probs = {k: float(raw_mapped.get(k, 0.0)) for k in STATES_3Q}
            probs = _normalize_out(probs)

            if cooldown_s > 0:
                time.sleep(cooldown_s)
            return probs

        except Exception as e:
            last_err = e
            sleep_s = base_sleep * (2 ** (attempt - 1))
            sleep_s *= (1.0 + random.uniform(-jitter, jitter))
            sleep_s = max(0.5, sleep_s)
            print(f"[NMR] attempt {attempt}/{max_tries} failed for '{name}': {e}")
            if attempt < max_tries:
                print(f"[NMR] sleeping {sleep_s:.2f}s then retrying...")
                time.sleep(sleep_s)

    raise RuntimeError(f"NMR job '{name}' failed after {max_tries} attempts. Last error: {last_err}")


def run_nmr_repeated_avg(
    circ,
    base_name: str,
    shots: int,
    repeats: int,
    *,
    max_tries: int = 6,
    base_sleep: float = 2.0,
    jitter: float = 0.35,
    cooldown_s: float = 2.0,
) -> Dict[str, float]:
    """Run the same NMR circuit multiple times and return the averaged distribution."""
    outs: List[Dict[str, float]] = []
    for r in range(repeats):
        outs.append(
            run_nmr_probs_robust(
                circ,
                name=f"{base_name}_r{r+1}",
                shots=shots,
                max_tries=max_tries,
                base_sleep=base_sleep,
                jitter=jitter,
                cooldown_s=cooldown_s,
            )
        )

    avg = {s: 0.0 for s in STATES_3Q}
    for d in outs:
        for s in STATES_3Q:
            avg[s] += d[s]
    for s in STATES_3Q:
        avg[s] /= max(1, len(outs))

    return _normalize_out(avg)
