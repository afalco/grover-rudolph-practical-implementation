# experiments/

This folder contains reproducible, single-purpose experiment scripts built on top of the `gr/` Python package.
Each script is intentionally small so you can isolate effects and maintain a clean “paper → code → hardware” trace.

## Quick start

Run from the repo root.

Recommended (module mode, most robust):

```bash
python -m experiments.00_smoke_sim_only
```

Direct execution also works if you are in the repo root:

```bash
python experiments/00_smoke_sim_only.py
```

## Import troubleshooting

If you see import errors such as `ModuleNotFoundError: gr`, verify:

1. You are in the repo root, the folder containing `gr/` and `experiments/`.
2. `gr/__init__.py` exists.
3. Prefer module execution: `python -m experiments.<script_name_without_py>`.

## Hardware setup (SpinQ Triangulum NMR)

Hardware runs typically require these environment variables:

- `SPINQ_IP`
- `SPINQ_PORT`
- `SPINQ_ACCOUNT`
- `SPINQ_PASSWORD`

For backward compatibility, some older code/comments may still mention `SPINQ_USER` and `SPINQ_PASS`. In this repository, the preferred names are:

- `SPINQ_ACCOUNT`
- `SPINQ_PASSWORD`

### Linux / macOS

```bash
export SPINQ_IP="192.0.2.1"
export SPINQ_PORT="12345"
export SPINQ_ACCOUNT="userX"
export SPINQ_PASSWORD="********"
```

### Windows PowerShell

```powershell
$env:SPINQ_IP = "192.0.2.1"
$env:SPINQ_PORT = "12345"
$env:SPINQ_ACCOUNT = "userX"
$env:SPINQ_PASSWORD = "********"
```

## Bit-order convention

For this Grover–Rudolph workflow, the experimentally validated comparison setting is:

```python
SIMULATOR_ORDER = "MSB"
HARDWARE_ORDER = "MSB"
```

So the effective comparison convention is **MSB for both simulator and hardware**.

That same canonical state order is used throughout the experiment scripts:

- bitstrings are interpreted as `|q0 q1 q2⟩`,
- printed state order is always `['000','001','010','011','100','101','110','111']`,
- target, simulator, hardware, and readout-mitigation artifacts should be interpreted in that same canonical order.

## Output conventions

- Generated outputs go to `artifacts/` (CSV/JSONL/PNG/PDF depending on script).
- Readout-calibration matrices and mitigated distributions are interpreted in the same canonical state order as the rest of the repository.
- If you compare NMR against SIM or target, do so in the canonical MSB state order used above.

## Scripts

### `00_smoke_sim_only.py`

Sanity check: Target vs SIM FULL should yield TV ≈ 0 and Fidelity ≈ 1.

Use as a pre-commit / CI guard.

### `01_bit_order_check.py`

Applies `X` to `q0`, `q1`, `q2` and checks the dominant bitstring in SIM and optionally NMR.

Confirms the convention `|q0 q1 q2⟩` in canonical MSB ordering.

### `02_trivial_suite.py`

Minimal diagnostic circuits:

- identity-safe circuit
- `Ry(π/2)` on each qubit
- Bell state on `(q0,q1)`

Prints per-qubit marginals, useful for readout bias and crosstalk checks.

### `03_ladder_full_A_vs_B.py`

Runs FULL Grover–Rudolph with UCRy ladder A and B.

Compares distributions, TV, L2, and fidelity in SIM and optionally NMR.

### `04_stage_ablation_L0_L01_FULL.py`

Runs `L0`, `L01`, and `FULL` and reports metrics versus SIM ideal and/or target.

Best first diagnostic to see at which depth performance degrades.

### `05_readout_calibration_8x8.py`

Calibrates the `8×8` readout confusion matrix `Mfull`.

Saves timestamped artifacts such as:

- `artifacts/Mfull_<timestamp>.npy`
- `artifacts/Mfull_<timestamp>.csv`

The matrix is indexed in the same canonical state order used across the repo.

### `06_readout_mitigation_apply.py`

Loads a saved `Mfull` and applies readout mitigation:

- raw → mitigated via ridge-regularized inversion
- optional mixing between raw and mitigated

Prints comparisons versus target and/or SIM ideal.

All raw and mitigated distributions are interpreted in the same canonical state order.

### `07_cx_only_stress_test.py`

Stress test: circuits that are ideally identity but contain many `CX` operations.

Helps detect sensitivity to entangling depth even when logical operations cancel.

### `08_full_gr_log_to_csv.py`

Full logging pipeline:

- staged runs (`L0`, `L01`, `FULL`)
- per-run plus averaged results
- optional readout mitigation
- metrics versus SIM ideal and versus target

Outputs:

- `artifacts/gr_log_<timestamp>.csv`
- `artifacts/gr_log_<timestamp>.jsonl`

## Recommended run order

1. `00_smoke_sim_only.py` — verify math and circuit in SIM
2. `01_bit_order_check.py` — confirm qubit/bitstring mapping
3. `04_stage_ablation_L0_L01_FULL.py` — staged hardware diagnostics
4. `03_ladder_full_A_vs_B.py` — choose ladder A or B
5. `05_readout_calibration_8x8.py` + `06_readout_mitigation_apply.py` — optional
6. `08_full_gr_log_to_csv.py` — final logged runs

## Notes on stability and debugging

If you see intermittent NMR connection failures:

- increase cooldown in `gr/backends.py` or the script’s cooldown parameter,
- reduce `REPEATS_NMR` and run multiple sessions.

If `FULL` degrades substantially relative to `L01`, try:

- switching UCRy ladder (`03_ladder_full_A_vs_B.py`),
- enabling command clipping (`CLIP_THETA_CMD`, e.g. `0.95 * π`),
- recalibrating the readout matrix (`05_readout_calibration_8x8.py`),
- applying mitigation and optionally mixing (`06_readout_mitigation_apply.py`).

## Expected outputs (quick validation checklist)

This section gives rough expected behavior to confirm the software stack is working end-to-end.

Exact values on NMR will vary with calibration, drift, and shot noise, but the qualitative patterns should match.

### Common conventions

- State order is always: `['000','001','010','011','100','101','110','111']`
- Bitstrings are MSB→LSB in the effective comparison convention used in this workflow.
- For TV/L2/Fidelity:
  - TV is in `[0,1]`, smaller is better.
  - L2 is in `[0,√2]`, smaller is better.
  - Classical fidelity is in `[0,1]`, larger is better.

### `00_smoke_sim_only.py`

Goal: validate the Grover–Rudolph math and circuit build in simulation.

Expected (SIM):

- Target vs SIM FULL:
  - TV ≈ 0, typically `< 1e-9`
  - L2 ≈ 0
  - Fidelity ≈ 1

If it fails:

- angle convention mismatch (`asin`/`acos`) or wrong conditional branch,
- wrong bit-ordering when mapping `prob8` to `|q0 q1 q2⟩`,
- bug in UCRy ladder coefficients or ordering.

### `01_bit_order_check.py`

Goal: verify the mapping between qubit indices and output bitstrings.

Expected (SIM):

- `X` on `q0` → dominant `100`
- `X` on `q1` → dominant `010`
- `X` on `q2` → dominant `001`

Expected (NMR):

- same dominant bitstrings as SIM, but with leakage to other states,
- dominant probability often in the ~0.8–0.95 range, device-dependent.

If dominant states differ, you have a real bit-order or qubit-index mismatch in the backend labeling.

### `02_trivial_suite.py`

Goal: basic health check for readout bias and single-qubit control.

Includes:

- identity-safe circuit
- `Ry(π/2)` on each qubit
- Bell on `(q0,q1)`

Expected (SIM):

- Identity → `000: 1.0`
- `Ry(π/2)` on `qi` → exactly 50/50 on that qubit, others fixed at 0
- Bell `(q0,q1)` → `000: 0.5`, `110: 0.5`

Expected (NMR):

- Identity dominated by `000` but not exactly 1.0
- `Ry(π/2)` gives a near 50/50 marginal on that qubit
- Bell shows correlation with peaks at `000` and `110`, but with extra mass elsewhere

### `03_ladder_full_A_vs_B.py`

Goal: compare the two UCRy ladder orderings on the FULL circuit.

Expected (SIM):

- Ladder A and B should both match target if both implementations are correct.

Expected (NMR):

- A and B can differ due to pulse scheduling or CX sensitivity.
- Compare TV/L2/Fidelity versus SIM ideal or target and pick the ladder that improves metrics.

### `04_stage_ablation_L0_L01_FULL.py`

Goal: pinpoint where performance degrades as depth increases.

Stages:

- `L0`: only level-0 rotation on `q0`
- `L01`: adds level-1 controlled rotations on `q1`
- `FULL`: adds level-2 UCRy on `q2`

Expected (NMR):

- often `L0` best, `L01` slightly worse, `FULL` worst, though this is device-dependent.

### `05_readout_calibration_8x8.py`

Goal: estimate the full readout confusion matrix `Mfull`.

Expected:

- diagonal entries typically largest,
- `cond(Mfull)` moderate,
- if `cond(Mfull)` is very large, inversion becomes unstable.

### `06_readout_mitigation_apply.py`

Goal: apply readout mitigation to a raw distribution.

Expected:

- mitigation may reduce obvious readout bias,
- but may also worsen metrics if `Mfull` is stale or drift is significant.

### `07_cx_only_stress_test.py`

Goal: test sensitivity to entangling depth using circuits that ideally cancel to identity.

Expected (SIM):

- output exactly `000: 1.0`

Expected (NMR):

- output remains biased toward `000`, but can drift with more CX layers.

### `08_full_gr_log_to_csv.py`

Goal: full logging run for paper-quality reporting.

Expected (SIM):

- near-perfect match to target for FULL.

Expected (NMR):

- non-trivial deviation from SIM/target,
- mitigation may slightly help or may not, depending on drift and `Mfull`.

Tip: prefer averaging across repeats and reporting mean ± std from the CSV.
