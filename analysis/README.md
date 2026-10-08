# `analysis/` — two robustness analyses for the compilation-sensitivity paper

Two questions the manuscript raises but does not answer from the campaign
tables alone. Each is a single self-contained script with no dependency beyond
NumPy and SciPy (the second also imports the exact simulator from the
`gr-triangulum-verification` checkout). Both write their full output to JSON
beside the script.

```
python3 01_serial_dependence.py  --json serial_dependence.json
python3 02_angle_noise_model.py  --json angle_noise_model.json
```

---

## 1. Serial dependence in the paired differences

`01_serial_dependence.py`

**The question.** The October intervention campaign reports five paired
contrasts over 25 repetitions. Pairing buys tolerance to arbitrary dependence
*within* a repetition — which is what the undocumented shot count puts in
doubt — but the paired *t* and the signed-rank test still assume the 25
differences are independent, or at least exchangeable, *across* repetitions.
The acquisition order was counterbalanced rather than randomised, so nothing
in the design guarantees it. A drift slow on the scale of one repetition but
fast on the scale of the session would violate it.

**What it does.** For each contrast: autocorrelation at lags 1–5 with Bartlett
bands, Ljung–Box, Durbin–Watson, a trend in repetition index, and three
dependence-robust alternatives to the iid interval — a Newey–West HAC
standard error, a moving-block bootstrap (block length `n^(1/3)` ≈ 3), and a
block sign-flip permutation test.

**What it found** (20 000 replicates, seed 20261007):

| contrast | mean | acf(1) | Ljung–Box *p* | trend *p* | HAC/iid | verdict |
|---|---|---|---|---|---|---|
| A−B, Control | +0.0395 | +0.175 | 0.561 | 0.844 | 1.16 | unchanged |
| A−B, Reduced | −0.0229 | +0.209 | 0.675 | **0.003** | 1.15 | **changes** |
| Reduction, ladder A | +0.0140 | +0.338 | 0.127 | 0.613 | 1.26 | unchanged |
| Reduction, ladder B | +0.0765 | +0.134 | 0.692 | **0.012** | 1.13 | unchanged |
| Ladder × reduction | +0.0625 | +0.194 | 0.524 | 0.055 | 1.17 | unchanged |

No autocorrelation exceeds the Bartlett band (±0.392), and no Ljung–Box test
rejects — though with 25 points that is weak evidence of absence, and the
Durbin–Watson of 1.03 on the ladder-A reduction is low enough to be worth
noting. Standard errors inflate by 13–26 % under HAC, which is modest.

The conclusions hold:

| contrast | iid 95 % CI | HAC 95 % CI | block bootstrap 95 % CI | block-flip *p* |
|---|---|---|---|---|
| A−B, Control | [+0.0165, +0.0626] | [+0.0141, +0.0650] | [+0.0111, +0.0619] | 0.023 |
| A−B, Reduced | [−0.0431, −0.0027] | [−0.0449, −0.0009] | [−0.0431, **+0.0001**] | 0.050 |
| Reduction, ladder A | [−0.0131, +0.0411] | [−0.0183, +0.0464] | [−0.0132, +0.0489] | 0.465 |
| Reduction, ladder B | [+0.0573, +0.0957] | [+0.0558, +0.0972] | [+0.0538, +0.0963] | 0.009 |
| **Ladder × reduction** | **[+0.0306, +0.0944]** | **[+0.0269, +0.0980]** | **[+0.0254, +0.0928]** | **0.004** |

**Conclusion.** The interaction — the paper's central result — survives every
dependence-robust treatment: the block bootstrap interval still excludes zero
and the block sign-flip *p* is 0.004. The one contrast that does not is the
A/B difference in the Reduced configuration, whose block-bootstrap interval
reaches +0.0001 and whose block-flip *p* is 0.050. That is the same contrast
that fails a Holm correction across the five rows, so the sign-reversal claim
is weakly supported for two independent reasons and the manuscript reports it
as an estimate rather than as a finding.

A secondary observation: the A−B Reduced and ladder-B reduction series carry
significant trends in repetition index (*p* = 0.003 and 0.012), so the
instrument was drifting within the session even though the contrasts are
robust to it.

---

## 2. Does an independent per-gate angle error predict the measured effects?

`02_angle_noise_model.py`

**The question.** The simplest mechanical account of compilation sensitivity
is that each commanded rotation is delivered with an angle error. An earlier
draft of the manuscript argued that no model of this family could produce an
A/B difference *at all*, because the perturbation is exchangeable over gates
and the two ladders command identical angle multisets. **That argument is
wrong**: a perturbation at a given position propagates through the CNOT
pattern that follows it, and the two ladders differ in exactly that pattern.
The question is quantitative, so it is simulated here rather than argued.

**What it does.** Evolves the exact compiled circuits with
`Ry(φ) → Ry(φ + ε)` on every rotation, for the four cells of the October
campaign — ladders A and B × Control (21 gates, tail included) and Reduced
(18 gates) — with the tail placed as the instrument's record has it: adjacent
at the end for ladder A, split by the final CNOT for ladder B. Each
configuration is checked to reproduce the target exactly at ε = 0 before any
noise is added. Two error structures are run: `independent` (ε redrawn for
every gate of every run) and `systematic` (one ε per gate shared across the
four cells within a run, redrawn between runs), the second being the structure
the campaign's own angle measurement actually supports.

**What it found** (20 000 trials per cell per σ; measured r.m.s. angle error
at the level-1 gate is 3.73°):

| quantity | σ=1° | σ=2° | σ=3.73° | σ=6° | σ=10° | measured |
|---|---|---|---|---|---|---|
| F, Control A | 0.9992 | 0.9967 | 0.9884 | 0.9705 | 0.9206 | 0.9669 |
| F, Control B | 0.9992 | 0.9967 | 0.9885 | 0.9704 | 0.9208 | 0.8481 |
| A−B, Control | +0.0000 | −0.0000 | −0.0001 | +0.0001 | −0.0002 | **+0.0395** |
| A−B, Reduced | +0.0000 | +0.0000 | +0.0001 | +0.0002 | +0.0001 | **−0.0229** |
| Reduction, ladder A | +0.0002 | +0.0009 | +0.0032 | +0.0080 | +0.0209 | +0.0140 |
| Reduction, ladder B | +0.0002 | +0.0009 | +0.0030 | +0.0079 | +0.0206 | +0.0765 |
| **Ladder × reduction** | +0.0000 | −0.0000 | −0.0002 | −0.0001 | −0.0004 | **+0.0625** |

Monte Carlo standard error on the interaction is 9×10⁻⁵ at σ = 3.73°. The
`systematic` variant gives the same means with a far smaller standard error,
confirming that the result is a symmetry of the model and not Monte Carlo
noise.

**Conclusion, in three parts.**

1. *The model is not degenerate.* With one ε vector shared by both ladders,
   the per-run A−B difference has a standard deviation of 3.0×10⁻⁴ at
   σ = 3.73° and reaches 2.4×10⁻³ in the tail. The perturbation does
   propagate differently through the two ladders, exactly as the exchangeability
   argument denied.
2. *Its mean A/B contrast is nevertheless zero.* At every σ from 1° to 10°,
   and under both error structures, the predicted A−B contrast and the
   predicted ladder × intervention interaction are zero to Monte Carlo
   precision — |·| ≤ 0.0002 against +0.0395 and +0.0625 measured, a factor of
   200 to 300.
3. *Its predicted spread is too small as well.* The run-to-run standard
   deviation it predicts for these differences at σ = 3.73° is 0.010 against
   0.056–0.077 measured, so the model cannot produce the observed contrasts
   even as fluctuations.

The model does reproduce the *direction* of the reduction effect — fewer gates,
fewer error sources, +0.0032 at σ = 3.73° — but symmetrically in the two
ladders, which is precisely the asymmetry the campaign measures. And no single
σ fits both ladders: ladder A's measured fidelity of 0.967 corresponds to
σ ≈ 6°, while ladder B's 0.848 lies below the model's value at σ = 10°.

An independent per-gate angle error is therefore not excluded by an argument
from symmetry, but it is excluded by simulation: it predicts none of the three
compilation effects the paper reports.

---

### Reproducibility

Both scripts are deterministic given `--seed` (default 20261007). Inputs are
`data/ladder_ab/D1_ladder_ab_runs.csv` in this repository and the `grtri`
package of the `gr-triangulum-verification` checkout; neither needs hardware
or the vendor SDK. The JSON outputs beside the scripts are the ones the
figures above were read from.
