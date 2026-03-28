# Grover–Rudolph Practical Implementation

This repository provides a **paper-to-hardware** implementation of the **Grover–Rudolph state preparation algorithm** for **3 qubits**, including:

- **Ideal simulator verification** (exact distribution matching)
- **Execution on a real quantum device**: *SpinQ Triangulum (NMR)*
- **Ancilla-free transpilation** into the native gate set `{Ry(·), X, CNOT}`
- **Staged diagnostics** (`L0` / `L01` / `FULL`) to localize hardware error accumulation
- **Metrics**: Total Variation (TV), L2 distance, and classical fidelity
- Optional: **8×8 readout calibration and mitigation**

## Reference paper (arXiv)

This codebase is intended as a practical companion to:

**A Rigorous and Self–Contained Proof of the Grover–Rudolph State Preparation Algorithm**  
arXiv:2601.17930

## What this repo demonstrates

### 1) Correctness in the ideal model

For a target distribution `p(x)` over 3-bit strings `x ∈ {0,1}^3`, Grover–Rudolph prepares the quantum state

\[
|\psi\rangle = \sum_{x\in\{0,1\}^3} \sqrt{p(x)}\,|x\rangle.
\]

The script `src/definitive_gr.py` includes a sanity check:

- **SIM FULL vs Target**: TV ≈ 0 and Fidelity ≈ 1

This is the computational counterpart of the paper’s formal correctness guarantee.

### 2) A hardware-friendly compilation (no nested multi-controlled parametric gates)

Real backends may be sensitive to nested or parametric multi-controlled gates.

This repository implements the key circuit-theoretic idea:

- each Grover–Rudolph “level” is a **uniformly controlled rotation** (UCRy),
- the circuit is compiled **ancilla-free** using only **Ry + CNOT + X**.

This is essential for stable execution on the SpinQ Triangulum NMR device.

### 3) Real-device behavior and diagnostics

Hardware runs are reported with:

- staged execution (`L0`, `L01`, `FULL`)
- repeated runs and averaging
- quantitative metrics and per-qubit marginals
- optional readout mitigation via an 8×8 confusion matrix

## Paper → Code mapping

The following map is intended to make the repository easy to audit as a faithful implementation of arXiv:2601.17930.

### A) Dyadic probability tree and conditional masses → angle routines

**Paper concept:** probability mass aggregation over dyadic prefixes, producing conditional probabilities along a binary tree.

**Code:** `src/definitive_gr.py` and modularized versions in `gr/angles.py`

- `sum_indices_prob(...)`: computes prefix masses `P(prefix)`
- `angles_3q_asin_child1(...)`: builds the GR angle dictionary for 3 qubits

**Key implementation choice:** the unambiguous **child = 1** convention

\[
\theta(\text{prefix}) = 2\arcsin\sqrt{\frac{P(\text{prefix}+1)}{P(\text{prefix})}},
\]

so that applying `Ry(θ)` on `|0⟩` yields `P(1) = sin²(θ/2)`.

This prevents branch and bit swaps that often arise from mixing `asin` and `acos` conventions.

### B) Stage-by-stage GR circuit construction → `build_gr_circuit_3q(...)`

**Paper concept:** inductive circuit construction by levels; each level rotates the next qubit conditioned on the prefix.

**Code:** `src/definitive_gr.py`

- `build_gr_circuit_3q(..., depth="L0"|"L01"|"full", ...)`

Stages:

- `L0`: only the top rotation on `q0`
- `L01`: adds level-1 controlled rotations on `q1`
- `FULL`: adds level-2 uniformly controlled rotation on `q2`

### C) Ancilla-free transpilation to `{Ry, X, CNOT}`

**Paper concept:** uniformly controlled `Ry` blocks can be decomposed using Gray-code or ladder constructions and a linear transform of angles, requiring only single-qubit rotations and CNOTs.

**Code:** `src/definitive_gr.py`

- level 1 (single control):
  - `apply_cry_1ctrl_decomp(...)`: CRy decomposition using `Ry + CNOT`
- level 2 (two controls):
  - `apply_ucry_2ctrl(...)`: 2-control UCRy using only `Ry + CNOT`
  - `ladder_sign_matrix(...)` and `ucry_coeffs_from_thetas(...)`: build and solve the ladder sign system `S a = θ`

**Why this matters:** it avoids nested `ControlledGate(ControlledGate(Ry))` patterns that can compile poorly or fail on constrained backends.

### D) Empirical hardware evaluation → NMR backend + staged tests + metrics

**Paper motivation:** a mathematically correct algorithm must still be validated under hardware imperfections.

**Code:** `src/definitive_gr.py`

- NMR execution:
  - `run_nmr_probs_robust(...)`: retries + exponential backoff + cooldown
  - `run_nmr_repeated_avg(...)`: averaging across repeated runs
- Analysis:
  - `compare(...)`: TV / L2 / Fidelity
  - `per_qubit_marginals(...)`: marginal probabilities per qubit
- Optional:
  - `calibrate_readout_matrix_8x8(...)` and `mitigate_readout(...)`

## Repository layout

```text
Grover-Rudolph-Practical-Implementation/
├── README.md
├── LICENSE
├── CITATION.cff
├── environment.yml
├── requirements.txt
├── .gitignore
├── gr/                         # core library (importable package)
│   ├── __init__.py
│   ├── circuit.py              # build_gr_circuit_3q + staged variants
│   ├── angles.py               # Grover–Rudolph angle tree (asin convention, etc.)
│   ├── ucry.py                 # UCRy ladder (A/B) coefficients + decompositions
│   ├── backends.py             # simulator + NMR runners (robust retries)
│   ├── metrics.py              # TV/L2/Fidelity + marginals
│   ├── readout.py              # 8×8 Mfull calibration + mitigation
│   ├── constants.py            # STATES_3Q and shared constants
│   └── utils.py                # small helpers
├── src/                        # CLI / entry points
│   └── definitive_gr.py        # main runner script (SIM + optional NMR + metrics)
├── experiments/                # diagnostic and reproducibility scripts
│   ├── README.md
│   ├── 00_smoke_sim_only.py
│   ├── 01_bit_order_check.py
│   ├── 02_trivial_suite.py
│   ├── 03_ladder_full_A_vs_B.py
│   ├── 04_stage_ablation_L0_L01_FULL.py
│   ├── 05_readout_calibration_8x8.py
│   ├── 06_readout_mitigation_apply.py
│   ├── 07_cx_only_stress_test.py
│   └── 08_full_gr_log_to_csv.py
├── docs/
│   ├── HARDWARE_TRIANGULUM.md
│   ├── UCRY_LADDER.md
│   ├── TROUBLESHOOTING.md
│   └── REPRODUCIBILITY.md
├── artifacts/                  # generated outputs (ignored by git)
└── tests/
    ├── test_angles.py
    ├── test_ucry.py
    └── test_metrics.py
```

## Create the environment

```bash
conda env create -f environment.yml
conda activate spinq-gr
```

A pip fallback can also be used if appropriate for your setup.

## Bit-order convention

The canonical state ordering used throughout this repository is:

- `|q0 q1 q2⟩ = 000, 001, ..., 111` (**MSB → LSB**)

This is the reference convention used for:

- target probability vectors,
- simulator comparisons,
- reported metrics,
- CSV and JSON artifacts.

### Simulator vs hardware

In practice, the repository assumes:

- the simulator is already aligned with the canonical ordering;
- the SpinQ Triangulum NMR backend exports bitstrings effectively in reversed order.

Therefore, when using Triangulum hardware outputs, the measured bitstrings must be remapped back to the canonical repo order before comparison.

### Backend remapping

The modular backend layer supports this through the environment variable:

```bash
export SPINQ_BITORDER=LSB->MSB
```

This means that a hardware-reported bitstring is reversed before being interpreted in the canonical repository convention.

If needed, the identity convention can also be used:

```bash
export SPINQ_BITORDER=MSB->LSB
```

### Recommended Triangulum setting

For Triangulum runs, the recommended setting is:

```bash
export SPINQ_BITORDER=LSB->MSB
```

This ensures that hardware outputs are mapped consistently to the canonical ordering used everywhere else in the repository.



## Triangulum connection parameters

When running hardware experiments on the SpinQ Triangulum backend, the connection
parameters are typically provided through environment variables so that scripts do
not need to hard-code credentials.

The most common parameters are:

- `SPINQ_IP`
- `SPINQ_PORT`
- `SPINQ_ACCOUNT`
- `SPINQ_PASSWORD`

Depending on the script or backend wrapper, these may be read directly from the
environment or passed as command-line arguments after being defined.

### Bash / zsh

In Bash, zsh, or a similar Unix shell, define them as:

```bash
export SPINQ_IP=<TRIANGULUM_IP>
export SPINQ_PORT=55444
export SPINQ_ACCOUNT=<ACCOUNT>
export SPINQ_PASSWORD=<PASSWORD>
export SPINQ_BITORDER=LSB->MSB
```

Example:

```bash
export SPINQ_IP=192.168.1.25
export SPINQ_PORT=55444
export SPINQ_ACCOUNT=my_user
export SPINQ_PASSWORD='my_secret_password'
export SPINQ_BITORDER=LSB->MSB
```

You can then run the hardware workflow in the same shell session.

### PowerShell

In PowerShell, use:

```powershell
$env:SPINQ_IP = "<TRIANGULUM_IP>"
$env:SPINQ_PORT = "55444"
$env:SPINQ_ACCOUNT = "<ACCOUNT>"
$env:SPINQ_PASSWORD = "<PASSWORD>"
$env:SPINQ_BITORDER = "LSB->MSB"
```

Example:

```powershell
$env:SPINQ_IP = "192.168.1.25"
$env:SPINQ_PORT = "55444"
$env:SPINQ_ACCOUNT = "my_user"
$env:SPINQ_PASSWORD = "my_secret_password"
$env:SPINQ_BITORDER = "LSB->MSB"
```

### Temporary vs persistent variables

The commands above define the variables only for the current shell session.

If you want them to persist:

- in Bash/zsh, add the `export ...` lines to `~/.bashrc`, `~/.zshrc`, or the
  shell startup file you use;
- in PowerShell, add the `$env:...` assignments to your PowerShell profile if
  appropriate for your workflow.

### Security note

Avoid committing credentials to the repository or hard-coding them inside Python
files. Environment variables are preferable for local execution, and a separate
private credentials file can also be used if your local workflow requires it.

## Verification utility: bit-order calibration

To verify the backend bit-order convention experimentally, use:

```bash
python calibrate_bit_order.py --backend sim --shots 1024 --outdir artifacts
```

or on Triangulum:

```bash
python calibrate_bit_order.py \
  --backend triangulum \
  --ip <TRIANGULUM_IP> \
  --port 55444 \
  --account <ACCOUNT> \
  --password <PASSWORD> \
  --shots 1024 \
  --outdir artifacts
```

This utility runs simple calibration circuits with known `X` flips and reports the effective measurement bit-order of the backend.

## Notes

The smoke-test script `gr_triangulum_smoke_test.py` already reflects the experimentally verified convention:

- simulator output treated as canonical,
- hardware output treated as reversed and canonized before comparison.

The recommended next step is to keep `src/definitive_gr.py` and any NMR-facing scripts aligned with the same backend remapping policy through `gr/backends.py`.
