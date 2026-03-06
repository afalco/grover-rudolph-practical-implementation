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

Below is a drop-in “Expected outputs” section you can append to /experiments/README.md. I’ve kept it practical: what files appear, what numbers should look like, and what deviations mean.

⸻

Expected outputs (quick validation checklist)

This section gives rough expected behavior to confirm the software stack is working end-to-end.
Exact values on NMR will vary with calibration, drift, and shot noise, but the qualitative patterns should match.

Common conventions
	•	State order is always:
['000','001','010','011','100','101','110','111']
	•	Bitstrings are MSB→LSB: |q0 q1 q2⟩.
	•	For TV/L2/Fidelity:
	•	TV is in [0,1], smaller is better.
	•	L2 is in [0,\sqrt{2}], smaller is better.
	•	Classical fidelity is in [0,1], larger is better.

⸻

00_smoke_sim_only.py

Goal: validate the Grover–Rudolph math + circuit build in simulation.

Expected (SIM):
	•	Target vs SIM FULL:
	•	TV ≈ 0 (typically < 1e-9)
	•	L2 ≈ 0
	•	Fidelity ≈ 1

If it fails (SIM not matching target):
	•	angle convention mismatch (asin/acos) or wrong conditional branch
	•	wrong bit-ordering when mapping prob8 to |q0 q1 q2⟩
	•	bug in UCRy ladder coefficients / ordering

⸻

01_bit_order_check.py

Goal: verify the mapping between qubit indices and output bitstrings.

Expected (SIM):
	•	X on q0 → dominant 100
	•	X on q1 → dominant 010
	•	X on q2 → dominant 001

Expected (NMR):
	•	same dominant bitstrings as SIM, but with leakage to other states
	•	dominant probability often in the ~0.8–0.95 range (device-dependent)

If dominant states differ:
	•	you have a true bit-order / qubit index mismatch in the backend labeling

⸻

02_trivial_suite.py

Goal: basic health check for readout bias and single-qubit control.

Includes:
	•	identity-safe circuit
	•	Ry(π/2) on each qubit
	•	Bell on (q0,q1)

Expected (SIM):
	•	Identity → 000: 1.0
	•	Ry(π/2) on qi → exactly 50/50 on that qubit, others fixed at 0
	•	Bell (q0,q1) → 000: 0.5, 110: 0.5

Expected (NMR):
	•	Identity dominated by 000 but not exactly 1.0 (readout noise)
	•	Ry(π/2) gives a near 50/50 marginal on that qubit (often 0.45–0.55)
	•	Bell shows correlation (peaks at 000 and 110), but with extra mass elsewhere

If identity shows large population away from 000:
	•	readout bias / drift; mitigate with 05_readout_calibration_8x8.py

⸻

03_ladder_full_A_vs_B.py

Goal: compare the two UCRy ladder orderings on the FULL circuit.

Expected (SIM):
	•	Ladder A and B should both match target if both implementations are correct.
	•	If SIM differs between A and B: one ladder implementation is wrong.

Expected (NMR):
	•	A and B can differ due to pulse scheduling / CX sensitivity.
	•	Compare TV/L2/Fidelity vs SIM ideal (or vs target); pick the ladder that improves metrics.

Output:
	•	prints side-by-side distributions and metrics
	•	optionally writes artifacts (if enabled in script)

⸻

04_stage_ablation_L0_L01_FULL.py

Goal: pinpoint where performance degrades as depth increases.

Stages:
	•	L0: only level-0 rotation on q0
	•	L01: adds level-1 controlled rotations on q1
	•	FULL: adds level-2 UCRy on q2

Expected (SIM):
	•	each stage should match its ideal stage target (by construction)

Expected (NMR):
	•	typical trend: L0 best, L01 slightly worse, FULL worst
(but not always; depends on CX quality and ladder choice)

If FULL collapses to almost-uniform:
	•	strong indicator of CX/entangler sensitivity or a broken 2-control implementation
	•	try ladder switch, clipping, and run 07_cx_only_stress_test.py

⸻

05_readout_calibration_8x8.py

Goal: estimate the full readout confusion matrix M_{\mathrm{full}}.

Outputs (timestamped):
	•	artifacts/Mfull_<timestamp>.csv
	•	artifacts/Mfull_<timestamp>.npy

Expected:
	•	diagonal entries typically largest (not necessarily close to 1.0)
	•	cond(Mfull) should be moderate (often ~1–10; device-dependent)
	•	too large (e.g. > 50–100) → inversion becomes unstable

If cond(Mfull) is huge:
	•	increase shots, repeat calibration, or use stronger ridge in mitigation

⸻

06_readout_mitigation_apply.py

Goal: apply readout mitigation to a raw distribution.

Expected:
	•	mitigated distribution may reduce obvious readout bias (e.g., identity closer to 000)
	•	but mitigation can sometimes worsen metrics if M is stale / drifting

Rule of thumb:
	•	if mitigation improves identity and trivial circuits but worsens deep circuits,
the limitation is likely gate noise, not readout.

⸻

07_cx_only_stress_test.py

Goal: test sensitivity to entangling depth using circuits that ideally cancel to identity.

Expected (SIM):
	•	output is exactly 000: 1.0

Expected (NMR):
	•	output stays biased toward 000, but can drift with more CX layers
	•	strong drift away from 000 indicates CX pulse imperfections accumulating

If CX-only “identity” becomes highly non-uniform:
	•	CX calibration is the bottleneck; expect FULL GR to degrade too

⸻

08_full_gr_log_to_csv.py

Goal: full logging run for paper-quality reporting.

Outputs (timestamped):
	•	artifacts/gr_log_<timestamp>.csv
	•	artifacts/gr_log_<timestamp>.jsonl

Expected (SIM):
	•	near-perfect match to target for FULL
	•	per-stage SIM behaves deterministically given the build

Expected (NMR):
	•	non-trivial deviation from SIM/target
	•	TV often ~0.1–0.3 depending on depth and calibration
	•	mitigation may slightly help or may not, depending on drift and M_{\mathrm{full}}

Tip: prefer averaging across repeats and reporting mean ± std from the CSV.

⸻

If you paste this into /experiments/README.md, I can also tailor the numeric ranges (TV/L2/Fid per stage) to your actual historical runs (e.g., the ones you posted: L0/L01/FULL averages and their typical TVs).
