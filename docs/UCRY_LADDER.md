# docs/UCRY_LADDER.md

## Scope

This note documents the **2-control uniformly controlled $R_y$** layer used in the FULL 3-qubit Grover–Rudolph circuit, and explains:

- why we avoid nested parametric controlled gates,
- how the $R_y$–CNOT ladder is constructed,
- what “ladder A” vs “ladder B” means,
- how to validate ladder correctness in simulation,
- why the ladder choice can matter on real NMR hardware.

---

## 1) What is a 2-control UCRy?

At level 2 of Grover–Rudolph (for 3 qubits), we need a rotation on the target $q_2$ conditioned on the two control bits $(q_0,q_1)$.

Concretely, we require four conditional angles:

$$
\theta_{00},\ \theta_{01},\ \theta_{10},\ \theta_{11},
$$

such that the applied operation on $q_2$ equals:

- $R_y(\theta_{00})$ if $(q_0,q_1)=(0,0)$,
- $R_y(\theta_{01})$ if $(q_0,q_1)=(0,1)$,
- $R_y(\theta_{10})$ if $(q_0,q_1)=(1,0)$,
- $R_y(\theta_{11})$ if $(q_0,q_1)=(1,1)$.

This is a standard uniformly controlled rotation with 2 controls (often denoted UCRy).

---

## 2) Why not use a nested controlled gate?

Some constrained backends (or certain compiler paths) can struggle with a doubly-controlled parametric gate such as a “CCRy” expressed via nested controlled gates. Typical failure modes include:

- incorrect parameter binding,
- gate decomposition bugs,
- unexpected flattening of the output distribution.

To mitigate this, we implement the level-2 operation using **only** the native gate set:

- single-qubit $R_y$,
- CNOT,
- (and potentially $X$ elsewhere for 0-controls, though level 2 itself uses only CNOT+Ry).

---

## 3) Ladder form: four $R_y$ rotations separated by CNOTs

We implement UCRy with a fixed template:

$$
R_y(a_0)\ \cdot\ \text{CNOT}\ \cdot\ R_y(a_1)\ \cdot\ \text{CNOT}\ \cdot\ R_y(a_2)\ \cdot\ \text{CNOT}\ \cdot\ R_y(a_3)\ \cdot\ \text{CNOT},
$$

where the CNOTs are between one of the control qubits and the target $q_2$.

The key point: depending on the control values, the CNOTs toggle the computational basis on $q_2$, which effectively changes the sign of subsequent $R_y$ contributions on $q_2$. Therefore, for each control pattern $(q_0,q_1)$, the net effective rotation angle is a signed sum of $(a_0,a_1,a_2,a_3)$.

This yields a linear system:

$$
S a = \theta,
$$

where:

- $a = (a_0,a_1,a_2,a_3)^T$,
- $\theta = (\theta_{00},\theta_{01},\theta_{10},\theta_{11})^T$,
- $S \in \{\pm 1\}^{4\times 4}$ is a sign matrix determined solely by the CNOT ordering.

In code:

- `ladder_sign_matrix(ladder)` builds $S$,
- `ucry_coeffs_from_thetas(...)` solves $a = S^{-1}\theta$ numerically.

---

## 4) Ladder A vs Ladder B

The ladder is defined by the sequence of which control qubit drives each CNOT.

### Ladder A
CNOT control sequence is:

- $q_1, q_0, q_1, q_0$

i.e. `CX(q1,q2)`, `CX(q0,q2)`, `CX(q1,q2)`, `CX(q0,q2)`.

### Ladder B
CNOT control sequence is:

- $q_0, q_1, q_0, q_1$

i.e. `CX(q0,q2)`, `CX(q1,q2)`, `CX(q0,q2)`, `CX(q1,q2)`.

Both ladders are valid decompositions, but they induce **different** sign matrices $S$, and therefore different $(a_0,a_1,a_2,a_3)$ for the same target conditional angles $\theta$.

A correct implementation must compute the coefficients using the sign matrix consistent with the chosen ladder.

---

## 5) How to validate ladder correctness (simulation)

A robust validation is:

1. Fix a target distribution `PROB8`.
2. Compute the ideal Grover–Rudolph FULL output distribution with your chosen ladder.
3. Verify that the ideal simulator output matches the target (up to sampling noise).

In addition, you can run a “FULL A vs FULL B” test in simulation:

- if the solver and ladder definition are consistent, both ladder A and ladder B should produce the **same** ideal distribution (they implement the same conditional angles, just decomposed differently).

If simulation outputs differ between ladders, then:
- either the ladder CNOT ordering in code does not match the sign matrix used to compute $(a_0,\dots,a_3)$,
- or the conditional angles passed to UCRy are incorrect.

---

## 6) Why ladder choice can matter on NMR hardware

Even if ladder A and ladder B are mathematically equivalent, on real NMR hardware they can differ due to:

- pulse scheduling,
- crosstalk and accumulated phase errors,
- different ordering of entangling pulses,
- device drift during longer sequences.

Therefore we expose `UCRY_LADDER = "A"|"B"` to let experiments select the empirically more stable option.

A recommended workflow is:

1. Run staged tests L0 and L01 to ensure the front part is stable.
2. Run FULL with ladder A and ladder B (same shots and repeats).
3. Compare:
   - `SIM FULL` vs `NMR FULL` for each ladder,
   - `NMR FULL A` vs `NMR FULL B` directly using TV/L2/fidelity.

---

## 7) Implementation note: solving the sign system

We compute $a$ by solving $S a = \theta$ using standard linear algebra (exact for the $\pm 1$ matrix):

- this avoids hard-coding coefficient formulas that can silently become wrong when the ladder is changed,
- and makes the ladder definition explicit and auditable.

In `definitive_gr.py`:

- `ladder_sign_matrix(ladder)` constructs $S$ by tracking parity flips induced by each CNOT for each control pattern.
- `ucry_coeffs_from_thetas(...)` solves for $a$.

This ensures the decomposition remains correct as long as:

1. the CNOT sequence in `apply_ucry_2ctrl(...)` matches the ladder definition,
2. the sign matrix uses the same ladder sequence.

---

## 8) Summary checklist

To ensure a correct UCRy layer:

- confirm conditional angles $\theta_{00},\theta_{01},\theta_{10},\theta_{11}$ correspond to the intended Grover–Rudolph tree convention,
- confirm the ladder CNOT sequence matches the sign-matrix builder,
- confirm `SIM FULL` matches the target distribution,
- only then compare ladder A vs B on hardware (performance/stability, not correctness).
