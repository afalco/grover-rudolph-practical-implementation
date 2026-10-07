# D1 pilot session, 8 March 2026 — readout-mitigation comparison

This directory holds the run logs and the readout calibration matrix behind
Table `tab:readout-mitigation` of the manuscript

> **Hardware characterisation of exact Grover–Rudolph state preparation on a
> 3-qubit NMR quantum computer** — Falcó, Falcó-Pomares, Latorre & Lin.

The manuscript (Section 5.3) cites these files by name and states that they are
deposited in this repository. They previously existed only under `artifacts/`,
which is excluded by `.gitignore`, so nothing was actually published. They are
copied here, under a tracked path, to make that statement true. The originals
remain in `artifacts/` as working files.

## What the session was

Two acquisitions on the same device, independent of the 700-run benchmark
campaign analysed in Section 6:

| | `gr_log_20260308-124206` | `gr_log_20260308-131538` |
|---|---|---|
| hardware runs | 15 (5 per stage × L0/L01/FULL) | 15 (5 per stage × L0/L01/FULL) |
| ladder | B only | B only |
| shots | 2048 | 2048 |
| readout mitigation | none | applied |

Thirty hardware runs in total. Ladder B only, so the session supports **no A/B
comparison and no structural claim** — it is used in the paper for the single
result below. Each log additionally carries three aggregate rows (`repeats: 5`)
and three simulator reference rows at 200 000 shots.

## Files

| file | content |
|---|---|
| `gr_log_20260308-124206.{jsonl,csv}` | raw session, run level |
| `gr_log_20260308-131538.{jsonl,csv}` | mitigated session, run level |
| `gr_summary_20260308-*.{csv,tex}` | the per-session comparison tables |
| `Mfull_20260308-112738.{npy,csv}` | the 8×8 readout confusion matrix |
| `gr_results_analysis_20260308.{tex,pdf}` | the session's own analysis write-up |

`artifacts/Mfull_latest.npy` is byte-identical to `Mfull_20260308-112738.npy`
and is not duplicated here.

## The result it supports

From `gr_summary_20260308-131538.csv`, FULL stage, against the target law:

| post-processing | TV | ℓ² | fidelity |
|---|---|---|---|
| raw | 0.245 | 0.192 | 0.927 |
| mitigated (ridge = 1e-3) | 0.251 | 0.203 | 0.922 |
| mixed (λ = 0.3) | 0.247 | 0.195 | 0.926 |

Correcting the measurement layer does not rescue the FULL-stage output — it
degrades it marginally, consistently across the three post-processing variants.
The confusion matrix is well conditioned (condition number 1.369, recorded as
`mfull_cond` in every row of the log), so this is not an ill-posed inversion.
This is the experimental basis for the paper reporting raw rather than mitigated
distributions throughout the 700-run campaign.

## Reproducing the matrix

The calibration itself is produced by `experiments/05_readout_calibration_8x8.py`
and applied by `experiments/06_readout_mitigation_apply.py`.
