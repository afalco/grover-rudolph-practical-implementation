# experiments/

This folder contains **reproducible experiment scripts** built on top of the `gr/` package.  
Each script is intentionally small and single-purpose so you can isolate effects and keep a clean
“paper → code → hardware” trace.

---

## Quick start

From the repo root:

```bash
python experiments/00_smoke_sim_only.py
```

If you get import issues, make sure you run from the repo root and that `gr/__init__.py` exists.
You can also run modules with:

```bash
python -m experiments.00_smoke_sim_only
```

---

## Hardware setup (SpinQ Triangulum NMR)

Hardware runs require environment variables:

- `SPINQ_IP`
- `SPINQ_PORT`
- `SPINQ_USER`
- `SPINQ_PASS`

Linux/macOS:

```bash
export SPINQ_IP="192.0.2.1"
export SPINQ_PORT="12345"
export SPINQ_USER="userX"
export SPINQ_PASS="********"
```

Windows PowerShell:

```powershell
setx SPINQ_IP "192.0.2.1"
setx SPINQ_PORT "12345"
setx SPINQ_USER "userX"
setx SPINQ_PASS "********"
```

Restart the terminal after `setx`.

### RUN_NMR toggle

Most scripts contain a flag like:

- `RUN_NMR = False` → **simulation-only** (no hardware calls, no env vars required)
- `RUN_NMR = True`  → **hardware mode** (NMR backend used, env vars required)

---

## Output conventions

- Bitstrings are always in **MSB→LSB** order: `|q0 q1 q2⟩`.
- Printed state order is:
  `['000','001','010','011','100','101','110','111']`.

---

## Scripts

### 00_smoke_sim_only.py
Sanity check: **Target vs SIM FULL** should yield TV ≈ 0 and fidelity ≈ 1.  
Use as a pre-commit/CI guard.

### 01_bit_order_check.py
Applies `X` to `q0`, `q1`, `q2` and checks the dominant bitstring in SIM and optionally NMR.  
Confirms the convention `|q0 q1 q2⟩` = MSB→LSB.

### 02_trivial_suite.py
Minimal diagnostic circuits:
- identity-safe circuit
- `Ry(pi/2)` on each qubit
- Bell state on `(q0,q1)`
Prints per-qubit marginals (useful for leakage/crosstalk checks).

### 03_ladder_full_A_vs_B.py
Runs **FULL Grover–Rudolph** with UCRy ladder **A** and **B**.  
Compares distributions (TV/L2/fidelity) in SIM and optionally NMR.

### 04_stage_ablation_L0_L01_FULL.py
Runs `L0`, `L01`, `FULL` and reports metrics vs SIM ideal.  
Best first diagnostic to see at which depth performance degrades.

### 05_readout_calibration_8x8.py
Calibrates the full **8×8 readout confusion matrix** `Mfull`.  
Saves:
- `artifacts/Mfull_<timestamp>.npy`
- `artifacts/Mfull_<timestamp>.csv`

### 06_readout_mitigation_apply.py
Loads a saved `Mfull` and applies mitigation:
- raw → mitigated (ridge inversion)
- optional mixing (raw/mitigated)
Prints comparisons vs Target.

### 07_cx_only_stress_test.py
Stress test: circuits that are (ideally) identity but contain many `CX` operations.  
Helps detect sensitivity to entangling depth even when logical operations cancel.

### 08_full_gr_log_to_csv.py
Full logging pipeline:
- staged runs (`L0/L01/FULL`)
- per-run + averaged results
- optional readout mitigation
Outputs:
- `artifacts/gr_log_<timestamp>.csv`
- `artifacts/gr_log_<timestamp>.jsonl` (full distributions per record)

---

## Recommended run order

1. `00_smoke_sim_only.py` (verify math + circuit in SIM)
2. `01_bit_order_check.py` (confirm qubit/bitstring mapping)
3. `04_stage_ablation_L0_L01_FULL.py` (staged hardware diagnostics)
4. `03_ladder_full_A_vs_B.py` (choose ladder A or B)
5. `05_readout_calibration_8x8.py` + `06_readout_mitigation_apply.py` (optional)
6. `08_full_gr_log_to_csv.py` (final logged runs)

---

## Notes on stability

- If you see intermittent NMR connection failures, increase cooldown in `gr/backends.py`
  or reduce `REPEATS_NMR` and run multiple sessions.
- If `FULL` degrades relative to `L01`, try:
  - switching ladder (`03_ladder_full_A_vs_B.py`)
  - enabling clipping (`CLIP_THETA_CMD`, e.g. `0.95*pi`)
  - calibrating and applying readout mitigation (`05/06`)

