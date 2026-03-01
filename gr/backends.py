# gr/backends.py
from __future__ import annotations

import os
import random
import time
from typing import Dict, Optional, Tuple

from spinqit import NMRConfig, get_compiler
from spinqit import get_basic_simulator, BasicSimulatorConfig

from .constants import STATES_3Q


def _normalize_out(out: Dict[str, float]) -> Dict[str, float]:
    s = float(sum(out.values()))
    if s <= 0.0:
        return {k: 0.0 for k in STATES_3Q}
    return {k: float(out.get(k, 0.0)) / s for k in STATES_3Q}


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


def _get_nmr_backend():
    from spinqit import get_nmr
    return get_nmr()


def _env(var: str) -> Optional[str]:
    v = os.environ.get(var, None)
    if v is None:
        return None
    v = str(v).strip()
    return v if v else None


def load_nmr_credentials_from_env() -> Tuple[str, int, str, str]:
    """
    Read NMR connection from env vars:
      SPINQ_IP, SPINQ_PORT, SPINQ_USER, SPINQ_PASS
    """
    ip = _env("SPINQ_IP")
    port_s = _env("SPINQ_PORT")
    user = _env("SPINQ_USER")
    pw = _env("SPINQ_PASS")
    missing = [k for k, v in [("SPINQ_IP", ip), ("SPINQ_PORT", port_s), ("SPINQ_USER", user), ("SPINQ_PASS", pw)] if not v]
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
    Returns normalized probabilities dict over STATES_3Q.
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
            probs = {k: float(out.get(k, 0.0)) for k in STATES_3Q}
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
    """
    Run the same NMR circuit multiple times and return the averaged distribution.
    """
    outs = []
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
