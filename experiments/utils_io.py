---

## `experiments/utils_io.py`

```python
# experiments/utils_io.py
from __future__ import annotations

import json
import os
import time
from typing import Dict, Optional

import numpy as np

from gr import STATES_3Q


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def timestamp() -> str:
    return time.strftime("%Y%m%d-%H%M%S")


def save_probs_json(path: str, probs: Dict[str, float]) -> None:
    out = {s: float(probs.get(s, 0.0)) for s in STATES_3Q}
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, sort_keys=True)


def load_probs_json(path: str) -> Dict[str, float]:
    with open(path, "r", encoding="utf-8") as f:
        d = json.load(f)
    return {s: float(d.get(s, 0.0)) for s in STATES_3Q}


def save_matrix(path_npy: str, M: np.ndarray, path_csv: Optional[str] = None) -> None:
    np.save(path_npy, M)
    if path_csv is not None:
        np.savetxt(path_csv, M, delimiter=",")


def load_matrix(path_npy: str) -> np.ndarray:
    return np.load(path_npy)
