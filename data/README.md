# Deposited data

| path | content |
|---|---|
| `ladder_ab/` | the October 2026 re-acquisition: 100 runs, run level |
| `pilot_d1_20260308/` | the March 2026 D1 readout-mitigation session |
| `spinqit_export_ladderab_20261005.csv` | the instrument's own job database, exported after the October session |

## `ladder_ab/`

`D1_ladder_ab_runs.jsonl` (primary) and `D1_ladder_ab_runs.csv` (flat). One row
per hardware run: configuration, ladder, repeat, shot count, gate census,
measured and reference distributions over the eight outcomes, single-qubit
marginals, UTC start and finish. 100 rows, 25 per (configuration, ladder) cell.
Produced by `experiments/10_ladder_ab_campaign.py`.

The CSV reproduces the column order of `runs_flat_v2.csv` in the campaign
repository exactly, then appends the extra fields, so the verification scripts
read it unmodified.

## `spinqit_export_ladderab_20261005.csv`

The SpinQuasar job database of the CEU Triangulum, 1139 experiments spanning
30 January to 6 October 2026. It is a **second and independent record** of the
same executions: written by the instrument, not by our pipeline. For each job it
carries the compiled gate list actually executed, the measured populations and
the measured 8×8 density matrix.

### Cross-check against the run record

| check | result |
|---|---|
| runs in our record found in the export | 100 of 100 |
| measured distributions agree | max deviation 5.6×10⁻¹⁷ over all 100 runs |
| gate census, Control | 21 gates (11 Ry, 2 X, 8 two-qubit) on all 50 runs |
| gate census, Reduced | 18 gates (8 Ry, 2 X, 8 two-qubit) on all 50 runs |
| instrument clock offset | +2.000 to +2.001 h relative to UTC, consistently |
| job status | SUCCEED on every campaign job |

The gate census is the useful part: it confirms on the instrument side that the
Reduced configuration really executed 18 gates and not 21, independently of what
our compiler reported.

### Bit-order calibration

Jobs `BITORDER_X_q0/q1/q2`, executed 18:15–18:19 UTC, fifteen minutes before the
campaign, by `calibrate_bit_order.py`. Dominant outcomes `100`, `010`, `001`
respectively, confirming the canonical MSB→LSB convention
(`SPINQ_BITORDER=MSB->LSB`) on this unit with instrument-side evidence. Dominant
probabilities 0.657, 0.805 and 0.864 — the readout is visibly imperfect on q0,
but the assignment is unambiguous.

### Seven orphaned jobs — not part of the dataset

The export contains 107 jobs carrying campaign run names, against 100 recorded
runs. Seven names appear twice:

`D1_control_A_003`, `D1_control_B_003`, `D1_control_B_004`,
`D1_exp1_A_003`, `D1_exp1_A_004`, `D1_exp1_B_003`, `D1_exp1_B_004`

The duplicates are **not retries**. They were executed at 18:29–18:41 instrument
time, two hours and nineteen minutes before the corresponding runs of the
recorded campaign, by an earlier aborted attempt that was pointed at the
directory already holding the simulator dry run. Because the resume logic of
that version keyed on run name alone, it treated the eight simulator rows of
repeats 1 and 2 as complete, skipped them, and began acquisition at repeat 3.
The attempt lost its connection to the instrument after seven jobs and recorded
nothing.

The deposited campaign was then acquired from scratch into a clean directory, so
the dataset here is unaffected: every one of its 100 runs matches the **later**
of the two device jobs, to floating-point identity. The earlier seven are listed
here because they exist in the instrument's record and a reader cross-checking
by job name would otherwise find an unexplained duplicate. Their measured
distributions differ from the recorded ones by up to 0.18 in the largest
outcome, which is the ordinary run-to-run spread of this unit, not evidence of
anything else.

The acquisition script now refuses to continue when a log holds runs from a
different backend, so this cannot recur.
