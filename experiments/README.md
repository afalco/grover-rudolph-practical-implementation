# experiments/

Here’s a tightened, more “repo-quality” version. I kept your structure, but:
	•	makes module-style invocation the default (avoids import issues)
	•	clarifies hardware env vars + Windows note
	•	adds artifacts conventions and “what to do when things go wrong”
	•	fixes small wording and the truncated last bullet (“mitigation”)
	•	keeps everything single-source-of-truth and actionable

You can paste this as /experiments/README.md.

⸻

experiments/

This folder contains reproducible, single-purpose experiment scripts built on top of the gr/ Python package.
Each script is intentionally small so you can isolate effects and maintain a clean “paper → code → hardware” trace.

⸻

Quick start

Run from the repo root.

Recommended (module mode, most robust):

```bash
python -m experiments.00_smoke_sim_only
```

Direct execution also works if you are in the repo root:

```bash
python experiments/00_smoke_sim_only.py
```

Import troubleshooting

If you see import errors (e.g., ModuleNotFoundError: gr), verify:
	1.	You are in the repo root (the folder containing gr/ and experiments/).
	2.	gr/__init__.py exists.
	3.	Prefer module execution: python -m experiments.<script_name_without_py>.

⸻

Hardware setup (SpinQ Triangulum NMR)

Hardware runs require these environment variables:
	•	SPINQ_IP
	•	SPINQ_PORT
	•	SPINQ_USER
	•	SPINQ_PASS

Linux/macOS

```bash
export SPINQ_IP="192.0.2.1"
export SPINQ_PORT="12345"
export SPINQ_USER="userX"
export SPINQ_PASS="********"
```

Windows PowerShell

```poweshell
setx SPINQ_IP "192.0.2.1"
setx SPINQ_PORT "12345"
setx SPINQ_USER "userX"
setx SPINQ_PASS "********"
```

Restart the terminal after setx.

RUN_NMR toggle

Most scripts expose a flag like:

RUN_NMR = False  # simulation-only
# RUN_NMR = True # enable hardware runs (requires SPINQ_* env vars)

	•	RUN_NMR = False → simulation-only (no hardware calls, no env vars required)
	•	RUN_NMR = True  → hardware mode (NMR backend used, env vars required)

⸻

Output conventions
	•	Bitstrings are always MSB→LSB: |q0 q1 q2⟩.
	•	Printed state order is always:
['000','001','010','011','100','101','110','111'].
	•	Generated outputs go to artifacts/ (CSV/JSONL/PNG/PDF depending on script).

⸻

Scripts

00_smoke_sim_only.py

Sanity check: Target vs SIM FULL should yield TV ≈ 0 and Fidelity ≈ 1.
Use as a pre-commit / CI guard.

01_bit_order_check.py

Applies X to q0, q1, q2 and checks the dominant bitstring in SIM and optionally NMR.
Confirms the convention |q0 q1 q2⟩ = MSB→LSB.

02_trivial_suite.py

Minimal diagnostic circuits:
	•	identity-safe circuit
	•	Ry(π/2) on each qubit
	•	Bell state on (q0,q1)

Prints per-qubit marginals (useful for readout bias / crosstalk checks).

03_ladder_full_A_vs_B.py

Runs FULL Grover–Rudolph with UCRy ladder A and B.
Compares distributions (TV/L2/Fidelity) in SIM and optionally NMR.

04_stage_ablation_L0_L01_FULL.py

Runs L0, L01, FULL and reports metrics vs SIM ideal (and/or target).
Best first diagnostic to see at which depth performance degrades.

05_readout_calibration_8x8.py

Calibrates the 8×8 readout confusion matrix Mfull.
Saves (timestamped):
	•	artifacts/Mfull_<timestamp>.npy
	•	artifacts/Mfull_<timestamp>.csv

06_readout_mitigation_apply.py

Loads a saved Mfull and applies readout mitigation:
	•	raw → mitigated (ridge-regularized inversion)
	•	optional mixing (raw/mitigated)

Prints comparisons vs Target (and/or SIM ideal).

07_cx_only_stress_test.py

Stress test: circuits that are (ideally) identity but contain many CX operations.
Helps detect sensitivity to entangling depth even when logical operations cancel.

08_full_gr_log_to_csv.py

Full logging pipeline:
	•	staged runs (L0/L01/FULL)
	•	per-run + averaged results
	•	optional readout mitigation
	•	metrics vs SIM ideal and vs Target

Outputs:
	•	artifacts/gr_log_<timestamp>.csv (tables / summaries)
	•	artifacts/gr_log_<timestamp>.jsonl (full records, distributions, metadata)

⸻

Recommended run order
	1.	00_smoke_sim_only.py (verify math + circuit in SIM)
	2.	01_bit_order_check.py (confirm qubit/bitstring mapping)
	3.	04_stage_ablation_L0_L01_FULL.py (staged hardware diagnostics)
	4.	03_ladder_full_A_vs_B.py (choose ladder A or B)
	5.	05_readout_calibration_8x8.py + 06_readout_mitigation_apply.py (optional)
	6.	08_full_gr_log_to_csv.py (final logged runs)

⸻

Notes on stability and debugging
	•	If you see intermittent NMR connection failures:
	•	increase cooldown in gr/backends.py (or the script’s cooldown parameter)
	•	reduce REPEATS_NMR and run multiple sessions
	•	If FULL degrades substantially relative to L01, try:
	•	switching UCRy ladder (03_ladder_full_A_vs_B.py)
	•	enabling command clipping (CLIP_THETA_CMD, e.g. 0.95 * π)
	•	recalibrating the readout matrix (05_readout_calibration_8x8.py)
	•	applying mitigation and optionally mixing (06_readout_mitigation_apply.py)

