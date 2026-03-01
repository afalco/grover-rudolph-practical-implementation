# docs/HARDWARE_TRIANGULUM.md

## Scope

This note documents practical considerations when running `src/definitive_gr.py` on the **SpinQ Triangulum (NMR)** backend, focusing on:

- bitstring conventions and bit-order validation,
- common backend/network failure modes and mitigations,
- staged execution (L0 / L01 / FULL) as a diagnostic tool,
- readout bias and optional readout mitigation.

The objective is to make runs reproducible and to help interpret discrepancies between ideal simulation and hardware output.

---

## 1) Bit order and state labeling

### 1.1 Convention used in this repository

All probability vectors and outputs are interpreted in **MSB→LSB** order:

- bitstring $b_0 b_1 b_2$ corresponds to qubits $(q_0,q_1,q_2)$,
- basis ordering is $\{|000\rangle,|001\rangle,\dots,|111\rangle\}$.

This matches the reporting format used by SpinQit in our tests.

### 1.2 How to validate the bit order on hardware

Run a “single X” test for each qubit and check the dominant state:

- apply $X$ to $q_0$ $\Rightarrow$ dominant state should be `100`,
- apply $X$ to $q_1$ $\Rightarrow$ dominant state should be `010`,
- apply $X$ to $q_2$ $\Rightarrow$ dominant state should be `001`.

If your output matches these three checks, then the mapping `q0 q1 q2 ↔ MSB LSB` is confirmed.

---

## 2) Why we include an “identity-safe tail”

Some backends can throw internal graph/IR attribute errors when a circuit is “too trivial” (or when certain IR fields are omitted). To avoid this, we optionally append a tiny canceling sequence:

- $R_y(\varepsilon)$ followed by $R_y(-\varepsilon)$ on $q_0$.

This preserves the logical circuit but helps ensure the backend generates a well-formed internal representation.

In code this is controlled by `ensure_nmr_attrs=True`, and implemented via:

- `add_identity_safe_tail(circ, q)`.

---

## 3) Connection / backend instability and robust execution

### 3.1 Typical symptoms

You may see log lines such as:

- `Error sending message: invalid state`
- handshake timeouts / reconnect loops
- sporadic job failures that disappear on retry

These are typically not algorithmic errors, but connectivity/session instability.

### 3.2 Mitigations used in this repo

The hardware runner implements:

- retries with exponential backoff,
- randomized jitter on retry delays,
- cooldown between successful jobs.

These controls reduce the probability of overwhelming the backend and improve repeatability.

Key parameters:

- `NMR_MAX_TRIES`
- `NMR_BASE_SLEEP`
- `NMR_JITTER`
- `COOLDOWN_S`

---

## 4) Staged execution as a diagnostic tool

The Grover–Rudolph circuit is naturally staged:

- **L0**: only the rotation on $q_0$,
- **L01**: adds controlled rotations on $q_1$ conditioned on $q_0$,
- **FULL**: adds a 2-control uniformly controlled rotation on $q_2$ conditioned on $(q_0,q_1)$.

### 4.1 What to expect ideally

If the algorithm is correct, the simulator should match the target distribution up to sampling noise. In particular:

- after L0, only states `000` and `100` are populated,
- after L01, only the four states with $q_2=0$ are populated,
- after FULL, all 8 states can be populated according to the target.

### 4.2 How to interpret hardware deviations

A useful heuristic is: “where does the deviation start?”

- if L0 is already far off, suspect state initialization/readout bias or drift,
- if L0 is acceptable but L01 degrades, suspect single-control operations / calibration,
- if L0 and L01 look reasonable but FULL degrades, suspect depth/crosstalk in the final layer (UCRy).

---

## 5) Readout bias and optional readout mitigation

### 5.1 Observed readout bias

In NMR outputs, even “identity” or “basis-prep” circuits can yield non-zero probability mass on unintended bitstrings. This is consistent with:

- state preparation imperfections,
- measurement/readout bias,
- device drift.

A common check is to run a basis state such as $|000\rangle$ and inspect the measured distribution.

### 5.2 Full 8×8 readout calibration

We can calibrate an 8×8 confusion matrix $M$ by preparing each computational basis state and measuring the observed distribution:

$$
M_{i,j} \approx \Pr(\text{meas}=i \mid \text{prep}=j),
$$

where $i,j \in \{000,\dots,111\}$.

We then mitigate a measured distribution $p_{\text{meas}}$ via ridge-regularized inversion:

$$
p_{\text{true}} \approx (M + \lambda I)^{-1} p_{\text{meas}}.
$$

### 5.3 Important caveat

Mitigation can sometimes **worsen** TV/L2 metrics if:

- $M$ is noisy,
- the inversion amplifies noise,
- the run-to-run variance dominates.

A robust workflow is:

1. average raw distributions over several runs,
2. mitigate the averaged distribution,
3. compare raw vs mitigated vs a convex mix.

---

## 6) Practical recommendations for experiments

1. **Confirm bit order** once per environment (SIM + NMR).
2. Use staged runs **L0 / L01 / FULL** routinely.
3. Increase `REPEATS_NMR` and keep `COOLDOWN_S` non-zero.
4. If FULL collapses or becomes unstable, try:
   - switching UCRy ladder A/B,
   - enabling mild clipping `CLIP_THETA_CMD`,
   - lowering shots per job but increasing repeats (stability vs variance tradeoff).
5. If using readout mitigation, always report:
   - condition number $\kappa(M)$,
   - ridge $\lambda$,
   - raw vs mitigated metrics.

---

## 7) Reproducibility metadata

For each reported run, record:

- date/time,
- SpinQit version and Python version,
- `SHOTS_NMR`, `REPEATS_NMR`, `COOLDOWN_S`,
- ladder choice A/B,
- clipping value if enabled,
- whether readout mitigation was used, with $\kappa(M)$ and $\lambda$.

This makes hardware results interpretable and comparable over time.
