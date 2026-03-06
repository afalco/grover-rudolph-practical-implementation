# experiments/08_full_gr_log_to_csv.py
#
# Runs Grover–Rudolph (3q) in stages (L0/L01/FULL), optionally on SIM and/or NMR,
# optionally with 8x8 readout mitigation, and logs EVERYTHING to CSV:
#   - per-run distributions (raw + mitigated)
#   - averaged distributions (raw + mitigated)
#   - TV/L2/Fidelity vs SIM ideal and vs Target
#   - per-qubit marginals
#   - metadata: timestamps, ladder, shots, repeats, ridge, etc.
#
# Output:
#   artifacts/gr_log_<timestamp>.csv
#   artifacts/gr_log_<timestamp>.jsonl   (one JSON per record, full distributions)
#
# Hardware mode requires env vars: SPINQ_IP, SPINQ_PORT, SPINQ_USER, SPINQ_PASS

from __future__ import annotations

import csv
import json
import os
import time
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Tuple

import sys
from pathlib import Path
# Make "gr/" importable even if this script is launched from inside "experiments/".
try:
    import gr  # noqa: F401
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from gr import (
    STATES_3Q,
    target_from_prob8,
    build_gr_circuit_3q,
    run_sim_probs,
    run_nmr_probs_robust,
    run_nmr_repeated_avg,
    tv_l2_fidelity,
    per_qubit_marginals,
    calibrate_readout_matrix_8x8,
    mitigate_readout,
)

# ----------------------- CONFIG -----------------------

PROB8 = [1, 2, 3, 4, 4, 3, 2, 1]

RUN_SIM = True
RUN_NMR = False

# Stages to run
STAGES = ["L0", "L01", "full"]

# UCRy ladder selection: "A" or "B"
UCRY_LADDER = "B"

# Shots / repeats
SHOTS_SIM = 200000
SHOTS_NMR = 2048
REPEATS_NMR = 5

# Optional: clipping (usually None first)
CLIP_THETA_CMD: Optional[float] = None  # e.g. 0.95 * np.pi

# Optional readout mitigation
DO_READOUT_MITIGATION = False
SHOTS_RO = 4096
RIDGE = 1e-3

# NMR robust parameters (overrides defaults in gr/backends if desired)
NMR_MAX_TRIES = 6
NMR_BASE_SLEEP = 2.0
NMR_JITTER = 0.35
COOLDOWN_S = 2.0

# Output folder
OUTDIR = "artifacts"

# ----------------------- HELPERS -----------------------


def ts() -> str:
    return time.strftime("%Y%m%d-%H%M%S")


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def normalize_dict_local(d: Dict[str, float]) -> Dict[str, float]:
    s = float(sum(d.values()))
    if s <= 0.0:
        return {k: 0.0 for k in STATES_3Q}
    return {k: float(d.get(k, 0.0)) / s for k in STATES_3Q}


def dict_to_list(d: Dict[str, float]) -> List[float]:
    d = normalize_dict_local(d)
    return [float(d.get(k, 0.0)) for k in STATES_3Q]


def marginals_to_flat(m: List[Tuple[float, float]]) -> Dict[str, float]:
    # m[i] = (p0,p1) for qi
    out = {}
    for i, (p0, p1) in enumerate(m):
        out[f"q{i}_p0"] = float(p0)
        out[f"q{i}_p1"] = float(p1)
    return out


def safe_cond_number(M: np.ndarray) -> float:
    try:
        return float(np.linalg.cond(M))
    except Exception:
        return float("nan")


@dataclass
class LogRecord:
    timestamp: str
    backend: str  # "SIM" or "NMR"
    stage: str    # "L0"|"L01"|"FULL"
    ladder: str   # "A"|"B"

    shots: int
    repeats: int
    clip_theta_cmd: str

    readout_mitigation: bool
    ridge: float
    mfull_cond: float

    # Metrics vs target and vs sim-ideal(stage)
    tv_vs_target: float
    l2_vs_target: float
    fid_vs_target: float

    tv_vs_sim: float
    l2_vs_sim: float
    fid_vs_sim: float

    # Distributions (serialized separately in JSONL); in CSV we store as compact list strings
    probs_raw_list: str
    probs_mitig_list: str

    # Per-qubit marginals for raw/mitig
    q0_p0_raw: float
    q0_p1_raw: float
    q1_p0_raw: float
    q1_p1_raw: float
    q2_p0_raw: float
    q2_p1_raw: float

    q0_p0_mitig: float
    q0_p1_mitig: float
    q1_p0_mitig: float
    q1_p1_mitig: float
    q2_p0_mitig: float
    q2_p1_mitig: float


def make_record(
    *,
    backend: str,
    stage: str,
    ladder: str,
    shots: int,
    repeats: int,
    clip_theta_cmd: Optional[float],
    target: Dict[str, float],
    sim_ref: Dict[str, float],
    raw: Dict[str, float],
    mitig: Optional[Dict[str, float]],
    mfull_cond: float,
    readout_mitigation: bool,
    ridge: float,
) -> LogRecord:
    raw_n = normalize_dict_local(raw)
    mitig_n = normalize_dict_local(mitig) if mitig is not None else {k: 0.0 for k in STATES_3Q}

    # metrics
    tv_t, l2_t, fid_t = tv_l2_fidelity(target, mitig_n if mitig is not None else raw_n)
    tv_s, l2_s, fid_s = tv_l2_fidelity(sim_ref, mitig_n if mitig is not None else raw_n)

    # marginals
    m_raw = per_qubit_marginals(raw_n)
    m_mit = per_qubit_marginals(mitig_n)

    flat_raw = marginals_to_flat(m_raw)
    flat_mit = marginals_to_flat(m_mit)

    stage_u = stage.upper()
    if stage_u == "FULL":
        stage_u = "FULL"
    # In our scripts stage values are "L0","L01","full"; standardize:
    stage_u = stage_u if stage_u in ("L0", "L01", "FULL") else stage_u

    return LogRecord(
        timestamp=ts(),
        backend=backend,
        stage=stage_u,
        ladder=ladder.upper(),
        shots=int(shots),
        repeats=int(repeats),
        clip_theta_cmd=str(clip_theta_cmd) if clip_theta_cmd is not None else "None",
        readout_mitigation=bool(readout_mitigation),
        ridge=float(ridge),
        mfull_cond=float(mfull_cond),

        tv_vs_target=float(tv_t),
        l2_vs_target=float(l2_t),
        fid_vs_target=float(fid_t),

        tv_vs_sim=float(tv_s),
        l2_vs_sim=float(l2_s),
        fid_vs_sim=float(fid_s),

        probs_raw_list=str(dict_to_list(raw_n)),
        probs_mitig_list=str(dict_to_list(mitig_n)),

        q0_p0_raw=float(flat_raw["q0_p0"]),
        q0_p1_raw=float(flat_raw["q0_p1"]),
        q1_p0_raw=float(flat_raw["q1_p0"]),
        q1_p1_raw=float(flat_raw["q1_p1"]),
        q2_p0_raw=float(flat_raw["q2_p0"]),
        q2_p1_raw=float(flat_raw["q2_p1"]),

        q0_p0_mitig=float(flat_mit["q0_p0"]),
        q0_p1_mitig=float(flat_mit["q0_p1"]),
        q1_p0_mitig=float(flat_mit["q1_p0"]),
        q1_p1_mitig=float(flat_mit["q1_p1"]),
        q2_p0_mitig=float(flat_mit["q2_p0"]),
        q2_p1_mitig=float(flat_mit["q2_p1"]),
    )


def append_jsonl(path: str, obj: dict) -> None:
    with open(path, "a", encoding="utf-8") as f:
        f.write(json.dumps(obj) + "\n")


# ----------------------- MAIN -----------------------


def main() -> None:
    ensure_dir(OUTDIR)
    out_csv = os.path.join(OUTDIR, f"gr_log_{ts()}.csv")
    out_jsonl = os.path.join(OUTDIR, f"gr_log_{ts()}.jsonl")

    target = target_from_prob8(PROB8)

    # Optional readout matrix
    Mfull = None
    mcond = float("nan")
    if RUN_NMR and DO_READOUT_MITIGATION:
        print("Calibrating 8x8 readout matrix (Mfull)...")
        Mfull = calibrate_readout_matrix_8x8(shots=SHOTS_RO, base_name="RO_FULL")
        mcond = safe_cond_number(Mfull)
        print("Done. Condition number:", mcond)

    # Open CSV and write header
    with open(out_csv, "w", newline="", encoding="utf-8") as fcsv:
        writer = csv.DictWriter(fcsv, fieldnames=list(asdict(LogRecord(
            timestamp="",
            backend="",
            stage="",
            ladder="",
            shots=0,
            repeats=0,
            clip_theta_cmd="",
            readout_mitigation=False,
            ridge=0.0,
            mfull_cond=0.0,
            tv_vs_target=0.0,
            l2_vs_target=0.0,
            fid_vs_target=0.0,
            tv_vs_sim=0.0,
            l2_vs_sim=0.0,
            fid_vs_sim=0.0,
            probs_raw_list="",
            probs_mitig_list="",
            q0_p0_raw=0.0,
            q0_p1_raw=0.0,
            q1_p0_raw=0.0,
            q1_p1_raw=0.0,
            q2_p0_raw=0.0,
            q2_p1_raw=0.0,
            q0_p0_mitig=0.0,
            q0_p1_mitig=0.0,
            q1_p0_mitig=0.0,
            q1_p1_mitig=0.0,
            q2_p0_mitig=0.0,
            q2_p1_mitig=0.0,
        )).keys()))
        writer.writeheader()

        # Stage loop
        for stage in STAGES:
            stage_norm = stage if stage.lower() != "full" else "full"

            # SIM reference for this stage (always ideal circuit)
            sim_ref = None
            if RUN_SIM:
                circ_sim = build_gr_circuit_3q(
                    PROB8,
                    k_gain=[1.0, 1.0, 1.0],
                    clip_cmd=None,
                    depth=stage_norm,
                    ladder=UCRY_LADDER,
                    ensure_nmr_attrs=False,
                )
                sim_ref = run_sim_probs(circ_sim, shots=SHOTS_SIM)

                rec_sim = make_record(
                    backend="SIM",
                    stage=stage,
                    ladder=UCRY_LADDER,
                    shots=SHOTS_SIM,
                    repeats=1,
                    clip_theta_cmd=None,
                    target=target,
                    sim_ref=sim_ref,
                    raw=sim_ref,
                    mitig=None,
                    mfull_cond=mcond,
                    readout_mitigation=False,
                    ridge=RIDGE,
                )
                writer.writerow(asdict(rec_sim))
                append_jsonl(out_jsonl, {
                    "record": asdict(rec_sim),
                    "probs_raw": normalize_dict_local(sim_ref),
                    "probs_mitig": None,
                })

            # NMR runs for this stage
            if RUN_NMR:
                if sim_ref is None:
                    # still need a sim reference for metrics; compute quickly
                    circ_sim = build_gr_circuit_3q(
                        PROB8, [1.0, 1.0, 1.0], None, depth=stage_norm,
                        ladder=UCRY_LADDER, ensure_nmr_attrs=False,
                    )
                    sim_ref = run_sim_probs(circ_sim, shots=SHOTS_SIM)

                circ_nmr = build_gr_circuit_3q(
                    PROB8,
                    k_gain=[1.0, 1.0, 1.0],
                    clip_cmd=CLIP_THETA_CMD,
                    depth=stage_norm,
                    ladder=UCRY_LADDER,
                    ensure_nmr_attrs=True,
                )

                # Per-repeat raw results
                raw_runs: List[Dict[str, float]] = []
                for r in range(REPEATS_NMR):
                    raw_r = run_nmr_probs_robust(
                        circ_nmr,
                        name=f"GR_{stage.upper()}_{UCRY_LADDER}_r{r+1}",
                        shots=SHOTS_NMR,
                        max_tries=NMR_MAX_TRIES,
                        base_sleep=NMR_BASE_SLEEP,
                        jitter=NMR_JITTER,
                        cooldown_s=COOLDOWN_S,
                    )
                    raw_runs.append(raw_r)

                    mitig_r = mitigate_readout(raw_r, Mfull, ridge=RIDGE) if (Mfull is not None) else None

                    rec_r = make_record(
                        backend="NMR",
                        stage=stage,
                        ladder=UCRY_LADDER,
                        shots=SHOTS_NMR,
                        repeats=1,
                        clip_theta_cmd=CLIP_THETA_CMD,
                        target=target,
                        sim_ref=sim_ref,
                        raw=raw_r,
                        mitig=mitig_r if DO_READOUT_MITIGATION else None,
                        mfull_cond=mcond,
                        readout_mitigation=bool(DO_READOUT_MITIGATION),
                        ridge=RIDGE,
                    )
                    writer.writerow(asdict(rec_r))
                    append_jsonl(out_jsonl, {
                        "record": asdict(rec_r),
                        "probs_raw": normalize_dict_local(raw_r),
                        "probs_mitig": normalize_dict_local(mitig_r) if mitig_r is not None else None,
                    })

                # Averaged raw + averaged mitigated
                avg_raw = {s: 0.0 for s in STATES_3Q}
                for d in raw_runs:
                    for s in STATES_3Q:
                        avg_raw[s] += float(d.get(s, 0.0))
                for s in STATES_3Q:
                    avg_raw[s] /= max(1, len(raw_runs))
                avg_raw = normalize_dict_local(avg_raw)

                avg_mitig = mitigate_readout(avg_raw, Mfull, ridge=RIDGE) if (Mfull is not None) else None

                rec_avg = make_record(
                    backend="NMR",
                    stage=stage,
                    ladder=UCRY_LADDER,
                    shots=SHOTS_NMR,
                    repeats=REPEATS_NMR,
                    clip_theta_cmd=CLIP_THETA_CMD,
                    target=target,
                    sim_ref=sim_ref,
                    raw=avg_raw,
                    mitig=avg_mitig if DO_READOUT_MITIGATION else None,
                    mfull_cond=mcond,
                    readout_mitigation=bool(DO_READOUT_MITIGATION),
                    ridge=RIDGE,
                )
                writer.writerow(asdict(rec_avg))
                append_jsonl(out_jsonl, {
                    "record": asdict(rec_avg),
                    "probs_raw": normalize_dict_local(avg_raw),
                    "probs_mitig": normalize_dict_local(avg_mitig) if avg_mitig is not None else None,
                    "note": "AVERAGE_OVER_REPEATS",
                })

    print("Saved CSV:", out_csv)
    print("Saved JSONL:", out_jsonl)
    print("Tip: JSONL contains full dicts for plotting; CSV is for tables.")


if __name__ == "__main__":
    main()
