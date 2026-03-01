# Grover–Rudolph Practical Implementation

This repository provides a **paper-to-hardware** implementation of the **Grover–Rudolph state preparation algorithm** for **3 qubits**, including:

- **Ideal simulator verification** (exact distribution matching)
- **Execution on a real quantum device**: *SpinQ Triangulum (NMR)*
- **Ancilla-free transpilation** into the native gate set `{Ry(·), X, CNOT}`
- **Staged diagnostics** (L0 / L01 / FULL) to localize hardware error accumulation
- **Metrics**: Total Variation (TV), L2 distance, and classical fidelity
- Optional: **8×8 readout calibration and mitigation**

## Reference paper (arXiv)

This codebase is intended as a practical companion to:

**A Rigorous and Self–Contained Proof of the Grover–Rudolph State Preparation Algorithm**  
arXiv:2601.17930  
https://arxiv.org/abs/2601.17930

## What this repo demonstrates

### 1) Correctness in the ideal model
For a target distribution `p(x)` over 3-bit strings `x ∈ {0,1}^3`, Grover–Rudolph prepares the quantum state

\[
|\psi\rangle=\sum_{x\in\{0,1\}^3}\sqrt{p(x)}\,|x\rangle.
\]

The script `src/definitive_gr.py` includes a **sanity check**:

- **SIM FULL vs Target**: TV ≈ 0 and Fidelity ≈ 1

This is the computational counterpart of the paper’s formal correctness guarantee.

### 2) A hardware-friendly compilation (no nested multi-controlled parametric gates)
Real backends may be sensitive to nested/parametric multi-controlled gates. This repository implements the key circuit-theoretic idea:

- Each Grover–Rudolph “level” is a **uniformly controlled rotation** (UCRy)
- The circuit is compiled **ancilla-free** using only **Ry + CNOT + X**

This is essential for stable execution on the SpinQ Triangulum NMR device.

### 3) Real-device behavior and diagnostics
Hardware runs are reported with:

- staged execution (`L0`, `L01`, `FULL`)
- repeated runs + averaging
- quantitative metrics + per-qubit marginals
- optional readout mitigation via an 8×8 confusion matrix

## Paper → Code mapping (precise)

The following map is intended to make the repository easy to audit as a faithful implementation of arXiv:2601.17930.

### A) Dyadic probability tree and conditional masses → `angles.py` / angle routines
**Paper concept:** probability mass aggregation over dyadic prefixes, producing conditional probabilities along a binary tree.  
**Code:** `src/definitive_gr.py` (functionality is factored in `src/gr/angles.py` if you split modules)
- `sum_indices_prob(...)`: computes prefix masses \(P(\text{prefix})\)
- `angles_3q_asin_child1(...)`: builds the GR angle dictionary for 3 qubits

**Key implementation choice:** the **unambiguous “child=1” convention**
\[
\theta(\text{prefix}) = 2\arcsin\sqrt{\frac{P(\text{prefix}+1)}{P(\text{prefix})}},
\]
so that applying \(R_y(\theta)\) on \(|0\rangle\) yields \(P(1)=\sin^2(\theta/2)\).  
This prevents branch/bit swaps that often arise from mixing `asin` vs `acos` conventions.

---

### B) Stage-by-stage GR circuit construction → `build_gr_circuit_3q(...)`
**Paper concept:** inductive circuit construction by levels; each level rotates the next qubit conditioned on the prefix.  
**Code:** `src/definitive_gr.py`
- `build_gr_circuit_3q(..., depth="L0"|"L01"|"full", ...)`

Stages:
- `L0`: only the top rotation on `q0`
- `L01`: adds level-1 controlled rotations on `q1`
- `FULL`: adds level-2 uniformly controlled rotation on `q2`

---

### C) Ancilla-free transpilation to `{Ry, X, CNOT}` → `apply_cry_1ctrl_decomp` and UCRy
**Paper concept:** uniformly controlled \(R_y\) blocks can be decomposed using Gray-code / ladder constructions and a linear transform of angles, requiring only single-qubit rotations and CNOTs.  
**Code:** `src/definitive_gr.py`
- Level 1 (single control):
  - `apply_cry_1ctrl_decomp(...)`: CRy decomposition using **Ry + CNOT**
- Level 2 (two controls):
  - `apply_ucry_2ctrl(...)`: 2-control UCRy using only **Ry + CNOT**
  - `ladder_sign_matrix(...)` and `ucry_coeffs_from_thetas(...)`: build and solve the ladder’s sign system \(S a=\theta\)

**Why this matters:** it avoids nested `ControlledGate(ControlledGate(Ry))` patterns that can compile poorly or fail on constrained backends.

---

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
  README.md
  LICENSE
  CITATION.cff

  paper/                 # link-only reference (e.g., arXiv link + citation)
    README.md

  environment.yml
  requirements.txt       # optional (pip fallback)
  .gitignore

  gr/                    # core library (importable package)
    __init__.py
    circuit.py           # build_gr_circuit_3q + staged variants
    angles.py            # Grover–Rudolph angle tree (asin convention, etc.)
    ucry.py              # UCRy ladder (A/B) coefficients + decompositions
    backends.py          # simulator + NMR runners (robust retries)
    metrics.py           # TV/L2/Fidelity + marginals
    readout.py           # 8×8 Mfull calibration + mitigation
    constants.py         # STATES_3Q and shared constants
    utils.py             # small helpers (normalization, formatting)

  src/                   # CLI/entry points (thin wrappers)
    definitive_gr.py     # main runner script (SIM + optional NMR + metrics)

  experiments/           # diagnostic and reproducibility scripts
    README.md
    00_smoke_sim_only.py
    01_bit_order_check.py
    02_trivial_suite.py
    03_ladder_full_A_vs_B.py
    04_stage_ablation_L0_L01_FULL.py
    05_readout_calibration_8x8.py
    06_readout_mitigation_apply.py
    07_cx_only_stress_test.py
    08_full_gr_log_to_csv.py
    utils_io.py

  docs/                  # documentation (human-readable)
    HARDWARE_TRIANGULUM.md
    UCRY_LADDER.md
    TROUBLESHOOTING.md    # optional
    REPRODUCIBILITY.md    # optional

  artifacts/             # generated outputs (ignored by git)
    .gitkeep

  tests/                 # optional unit tests
    test_angles.py
    test_ucry.py
    test_metrics.py

  .github/               # optional automation
    workflows/
      ci.yml
```

### Create the environment

```bash
conda env create -f environment.yml
conda activate spinq-gr
```


