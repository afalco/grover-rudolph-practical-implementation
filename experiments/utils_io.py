# experiments/utils_io.py
from __future__ import annotations

import json
import os
import time
from pathlib import Path
from typing import Dict, Optional

import numpy as np

from gr import STATES_3Q


def ensure_dir(path: str) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def timestamp() -> str:
    return time.strftime("%Y%m%d-%H%M%S")


def _ensure_parent(file_path: str) -> None:
    Path(file_path).expanduser().resolve().parent.mkdir(parents=True, exist_ok=True)


def normalize_probs(probs: Dict[str, float]) -> Dict[str, float]:
    """Return a normalized probs dict over STATES_3Q (missing keys -> 0)."""
    out = {s: float(probs.get(s, 0.0)) for s in STATES_3Q}
    ssum = float(sum(out.values()))
    if ssum > 0:
        out = {k: v / ssum for k, v in out.items()}
    return out


def save_probs_json(path: str, probs: Dict[str, float], *, normalize: bool = True) -> None:
    _ensure_parent(path)
    out = normalize_probs(probs) if normalize else {s: float(probs.get(s, 0.0)) for s in STATES_3Q}
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, sort_keys=True)


def load_probs_json(path: str, *, normalize: bool = False) -> Dict[str, float]:
    with open(path, "r", encoding="utf-8") as f:
        d = json.load(f)
    out = {s: float(d.get(s, 0.0)) for s in STATES_3Q}
    return normalize_probs(out) if normalize else out


def save_matrix(path_npy: str, M: np.ndarray, path_csv: Optional[str] = None) -> None:
    _ensure_parent(path_npy)
    np.save(path_npy, M)
    if path_csv is not None:
        _ensure_parent(path_csv)
        np.savetxt(path_csv, M, delimiter=",")


def load_matrix(path: str) -> np.ndarray:
    """
    Load a matrix from .npy or .csv.
    """
    p = Path(path)
    if p.suffix.lower() == ".npy":
        return np.load(str(p))
    if p.suffix.lower() == ".csv":
        return np.loadtxt(str(p), delimiter=",")
    raise ValueError(f"Unsupported matrix file type: {p.suffix} (expected .npy or .csv)")
