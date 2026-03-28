# docs/HARDWARE_TRIANGULUM.md

## Scope

This note documents practical considerations when running `src/definitive_gr.py` on the **SpinQ Triangulum (NMR)** backend, focusing on:

- bitstring conventions and bit-order validation,
- connection parameters and shell setup,
- common backend/network failure modes and mitigations,
- staged execution (`L0` / `L01` / `FULL`) as a diagnostic tool,
- readout bias and optional readout mitigation.

The objective is to make runs reproducible and to help interpret discrepancies between ideal simulation and hardware output.

---

## 1) Bit order and state labeling

### 1.1 Canonical convention used in this repository

All probability vectors and reported outputs are interpreted in the **canonical repo order**

- bitstring `b0 b1 b2` corresponds to qubits `(q0, q1, q2)`,
- basis ordering is `|000⟩, |001⟩, ..., |111⟩`,
- this is the **MSB → LSB** convention.

This is the reference convention used for:

- target distributions,
- simulator comparisons,
- TV / L2 / fidelity metrics,
- JSON / CSV artifacts.

### 1.2 Simulator vs Triangulum hardware

For the present Grover–Rudolph workflow, the experimentally validated comparison
setting is:

```python
SIMULATOR_ORDER = "MSB"
HARDWARE_ORDER = "MSB"
```


Operationally, the repository assumes:

- the simulator is already aligned with the canonical ordering;
- the SpinQ Triangulum backend is currently best interpreted directly in the canonical ordering for this workflow.

So, for Triangulum hardware in this workflow, measured bitstrings should be compared directly in the same canonical repo order used for simulation.

The recommended backend setting is therefore:

```bash
export SPINQ_BITORDER=MSB->LSB
```

In PowerShell:

```powershell
$env:SPINQ_BITORDER = "MSB->LSB"
```

### 1.3 How to validate the bit order on hardware

Use the calibration utility:

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

The script applies known `X` flips and infers whether the backend reports effectively as:

- `q0 q1 q2` (canonical), or
- `q2 q1 q0` (reversed).

A manual spot check is also possible:

- apply `X` to `q0`,
- apply `X` to `q1`,
- apply `X` to `q2`,

and inspect the dominant reported bitstring after canonization.

---

## 2) Triangulum connection parameters

When running hardware experiments, it is best to define the connection parameters through environment variables rather than hard-coding them in Python files.

Typical variables are:

- `SPINQ_IP`
- `SPINQ_PORT`
- `SPINQ_ACCOUNT`
- `SPINQ_PASSWORD`
- `SPINQ_BITORDER`

### 2.1 Bash / zsh

```bash
export SPINQ_IP=<TRIANGULUM_IP>
export SPINQ_PORT=55444
export SPINQ_ACCOUNT=<ACCOUNT>
export SPINQ_PASSWORD=<PASSWORD>
export SPINQ_BITORDER=MSB->LSB
```

Example:

```bash
export SPINQ_IP=192.168.1.25
export SPINQ_PORT=55444
export SPINQ_ACCOUNT=my_user
export SPINQ_PASSWORD='my_secret_password'
export SPINQ_BITORDER=MSB->LSB
```

### 2.2 PowerShell

```powershell
$env:SPINQ_IP = "<TRIANGULUM_IP>"
$env:SPINQ_PORT = "55444"
$env:SPINQ_ACCOUNT = "<ACCOUNT>"
$env:SPINQ_PASSWORD = "<PASSWORD>"
$env:SPINQ_BITORDER = "MSB->LSB"
```

Example:

```powershell
$env:SPINQ_IP = "192.168.1.25"
$env:SPINQ_PORT = "55444"
$env:SPINQ_ACCOUNT = "my_user"
$env:SPINQ_PASSWORD = "my_secret_password"
$env:SPINQ_BITORDER = "MSB->LSB"
```

### 2.3 Security note

Do not commit credentials to the repository, and avoid hard-coding them in scripts. Environment variables are the preferred local workflow.

---

## 3) Why we include an “identity-safe tail”

Some backends can throw internal graph / IR attribute errors when a circuit is too trivial, or when certain IR fields are omitted.

To reduce this risk, we may append a tiny canceling sequence:

- `Ry(ε)` followed by `Ry(-ε)` on a qubit.

This preserves the logical circuit while helping the backend generate a well-formed internal representation.

In code this is controlled by `ensure_nmr_attrs=True`, and implemented via:

- `add_identity_safe_tail(circ, q)`.

---

## 4) Connection / backend instability and robust execution

### 4.1 Typical symptoms

You may see log lines such as:

- `Error sending message: invalid state`
- handshake timeouts / reconnect loops
- sporadic job failures that disappear on retry

These are typically not algorithmic errors, but connectivity or session instability.

### 4.2 Mitigations used in this repository

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

## 5) Staged execution as a diagnostic tool

The Grover–Rudolph circuit is naturally staged:

- **L0**: only the rotation on `q0`,
- **L01**: adds controlled rotations on `q1` conditioned on `q0`,
- **FULL**: adds a 2-control uniformly controlled rotation on `q2` conditioned on `(q0, q1)`.

### 5.1 What to expect ideally

If the algorithm is correct, the simulator should match the target distribution up to sampling noise. In particular:

- after `L0`, only states `000` and `100` are populated,
- after `L01`, only the four states with `q2 = 0` are populated,
- after `FULL`, all 8 states can be populated according to the target.

### 5.2 How to interpret hardware deviations

A useful heuristic is: **where does the deviation start?**

- if `L0` is already far off, suspect state initialization, readout bias, or drift;
- if `L0` is acceptable but `L01` degrades, suspect single-control operations or calibration;
- if `L0` and `L01` look reasonable but `FULL` degrades, suspect depth / crosstalk in the final UCRy layer.

---

## 6) Readout bias and optional mitigation

### 6.1 Observed readout bias

In NMR outputs, even identity or basis-preparation circuits can yield non-zero probability mass on unintended bitstrings. This is consistent with:

- state preparation imperfections,
- measurement / readout bias,
- device drift.

A standard check is to run a basis state such as `|000⟩` and inspect the measured distribution.

### 6.2 Full 8×8 readout calibration

We can calibrate an `8×8` confusion matrix `M` by preparing each computational basis state and measuring the observed distribution:

```math
M_{i,j} \approx \Pr(\text{meas}=i \mid \text{prep}=j),
```

where `i,j ∈ {000, ..., 111}`.

We then mitigate a measured distribution `p_meas` via ridge-regularized inversion:

```math
p_{\text{true}} \approx (M + \lambda I)^{-1} p_{\text{meas}}.
```

### 6.3 Important caveat

Mitigation can sometimes worsen TV / L2 metrics if:

- `M` is noisy,
- the inversion amplifies noise,
- run-to-run variance dominates.

A robust workflow is:

1. average raw distributions over several runs,
2. mitigate the averaged distribution,
3. compare raw vs mitigated vs a convex mix.

---

## 7) Practical recommendations for experiments

1. Confirm bit order once per environment.
2. Set `SPINQ_BITORDER=MSB->LSB` for Triangulum in this workflow unless you have strong evidence to override it.
3. Use staged runs `L0 / L01 / FULL` routinely.
4. Increase `REPEATS_NMR` and keep `COOLDOWN_S` non-zero.
5. If `FULL` collapses or becomes unstable, try:
   - switching UCRy ladder A/B,
   - enabling mild clipping `CLIP_THETA_CMD`,
   - lowering shots per job but increasing repeats.
6. If using readout mitigation, always report:
   - condition number `κ(M)`,
   - ridge `λ`,
   - raw vs mitigated metrics.

---

## 8) Reproducibility metadata

For each reported run, record:

- date / time,
- SpinQit version and Python version,
- `SHOTS_NMR`, `REPEATS_NMR`, `COOLDOWN_S`,
- ladder choice A/B,
- clipping value if enabled,
- whether readout mitigation was used, with `κ(M)` and `λ`,
- bit-order convention used,
- whether raw Triangulum outputs were canonized from `LSB->MSB`.

This makes hardware results interpretable and comparable over time.
