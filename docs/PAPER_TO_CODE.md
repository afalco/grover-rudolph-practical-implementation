# docs/PAPER_TO_CODE.md

## Scope

This document provides a **precise mapping** between the mathematical objects and constructions in:

**A Rigorous and Self–Contained Proof of the Grover–Rudolph State Preparation Algorithm**  
arXiv:2601.17930 — https://arxiv.org/abs/2601.17930

and the implementation in this repository (notably `src/definitive_gr.py`).

The goal is to make the repository auditable as a faithful “paper-to-hardware” realization.

---

## 1) Dyadic probability tree and prefix masses

### Paper object
The proof is organized around a **binary (dyadic) tree** over prefixes. Each node corresponds to a prefix $b_0\cdots b_{k-1}$ and carries a **mass**

$$
P(b_0\cdots b_{k-1}) = \sum_{x \in \{0,1\}^{3-k}} p(b_0\cdots b_{k-1}x).
$$

### Code mapping
- `sum_indices_prob(p, prefix_bits)` computes \(P(\text{prefix})\) by summing over all basis indices whose MSB prefix matches `prefix_bits`.
- Basis convention is **MSB→LSB**:
  $$|q_0 q_1 q_2\rangle \leftrightarrow \text{bitstring} \,b_0 b_1 b_2.$$

---

## 2) Angle map (conditional masses → rotation angles)

### Paper object
At each node (prefix), Grover–Rudolph defines a rotation on the next qubit so that the probability of branching to the next bit equals the conditional probability induced by the masses.

### Code convention (asin / child=1)
This repository uses the convention:

$$
\theta(\text{prefix}) = 2\,\arcsin\sqrt{\frac{P(\text{prefix}+1)}{P(\text{prefix})}},
\quad\Rightarrow\quad
\Pr(1 \mid \text{prefix}) = \sin^2(\theta/2).
$$

This avoids ambiguity between `asin` vs `acos` conventions and prevents common “branch swap” bugs.

### Code mapping
- `angles_3q_asin_child1(prob8)` returns the angles dictionary:
  - level 0: angle for $q_0$,
  - level 1: angles for $q_1$ conditioned on $q_0$,
  - level 2: angles for $q_2$ conditioned on $(q_0,q_1)$.

---

## 3) Circuit construction by levels (L0 / L01 / FULL)

### Paper object
The circuit is constructed inductively: at each level $k$, the algorithm applies a rotation on qubit $q_k$ conditioned on the prefix $q_0\cdots q_{k-1}$. The proof shows that after level $k$, the correct marginal distribution over prefixes of length \(k\) is obtained.

### Code mapping
- `build_gr_circuit_3q(..., depth="L0"|"L01"|"full", ...)`

Stages:
- **L0**: apply the level-0 rotation on $q_0$.
- **L01**: L0 plus level-1 controlled rotations on \(q_1\) conditioned on $q_0$.
- **FULL**: L01 plus the level-2 uniformly controlled rotation on $q_2$ conditioned on $(q_0,q_1)$.

These staged runs are also used as **hardware diagnostics**.

---

## 4) Ancilla-free transpilation into {Ry, X, CNOT}

### Paper object
A central practical message is that the Grover–Rudolph circuit can be implemented without ancillas using uniformly controlled rotations decomposed into single-qubit rotations and CNOT ladders (Gray-code based).

### Code mapping

#### 4.1 Single-control Ry (level 1)
- `apply_cry_1ctrl_decomp(...)` implements CRy using only Ry + CNOT:

$$\mathrm{CR}_y(\theta)=
R_y(\theta/2)\;\mathrm{CNOT}\;R_y(-\theta/2)\;\mathrm{CNOT}.$$

0-controls are handled with an X-sandwich on the control qubit.

#### 4.2 Two-control UCRy (level 2)
Level 2 requires four conditional angles:

$$
\theta_{00}, \theta_{01}, \theta_{10}, \theta_{11}.
$$

The code builds a Ry–CNOT ladder with **four Ry angles** $a_0,a_1,a_2,a_3$, whose effect depends on the parity toggled by the CNOT pattern. The code constructs the associated sign system:

$$
S a = \theta,
$$

and solves for \(a\) numerically.

- `ladder_sign_matrix(ladder)` builds $S$ for ladder type **A** or **B**.
- `ucry_coeffs_from_thetas(...)` solves for $(a_0,a_1,a_2,a_3)$.
- `apply_ucry_2ctrl(...)` applies the ladder using only Ry + CNOT.

**Ladder A/B**: The paper’s decomposition can be realized with different Gray-like orders. On real hardware, the choice can matter due to pulse scheduling/crosstalk. This repo exposes `UCRY_LADDER = "A"|"B"` to allow empirical comparison.

---

## 5) Ideal verification vs hardware evaluation

### Paper object
The proof guarantees the ideal correctness of the circuit. Hardware experiments provide an empirical measure of how close the prepared distribution is to the target.

### Code mapping
- Ideal simulator:
  - `run_sim_probs(...)` produces the ideal output distribution.
- Hardware (NMR):
  - `run_nmr_probs_robust(...)`: robust execution with retries/backoff.
  - `run_nmr_repeated_avg(...)`: averaging across repeated runs.

### Metrics
- `compare(a,b)` prints:
  - Total variation distance (TV)
  - L2 distance
  - classical fidelity $(\sum \sqrt{p_i q_i})^2$
- `print_marginals(...)` prints per-qubit marginals to detect systematic bias.

---

## 6) Optional modules: gain calibration and readout mitigation

These are not part of the theoretical core proof but support real-device analysis.

### 6.1 Minimal Ry gain calibration (single-point)
- `calibrate_ry_gain_single_point(qi)` estimates an effective amplitude scaling $k_i$ for each target qubit from an `Ry(pi/2)` experiment.

### 6.2 Readout mitigation (8×8)
- `calibrate_readout_matrix_8x8(...)` estimates a full readout confusion matrix.
- `mitigate_readout(...)` applies ridge-regularized inversion:
  $$
  p_{\text{true}} \approx (M + \lambda I)^{-1} p_{\text{meas}}.
  $$

---

## Reproducibility checklist

To reproduce a paper-to-hardware figure/table from this repo:

1. Fix `PROB8` and `UCRY_LADDER`.
2. Run staged execution: `L0`, `L01`, `FULL`.
3. Average over `REPEATS_NMR` runs with cooldown.
4. Report TV/L2/fidelity and per-qubit marginals.
5. (Optional) Calibrate/readout-mitigate and report raw vs mitigated.
