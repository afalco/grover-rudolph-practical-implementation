# Deposited data

| path | content |
|---|---|
| `ladder_ab/` | the 5 October 2026 intervention campaign: 100 runs, run level |
| `angle_convention/` | the 7 October 2026 angle-representative campaign: 100 runs, run level |
| `pilot_d1_20260308/` | the March 2026 D1 readout-mitigation session |
| `spinqit_export_ladderab_20261005.csv` | the instrument's job database, exported after the 5 October session |
| `spinqit_export_angle_20261007.csv` | the same database, exported after the 7 October session |

## `ladder_ab/`

`D1_ladder_ab_runs.jsonl` (primary) and `D1_ladder_ab_runs.csv` (flat). One row
per hardware run: configuration, ladder, repeat, shot count, gate census,
measured and reference distributions over the eight outcomes, single-qubit
marginals, UTC start and finish. 100 rows, 25 per (configuration, ladder) cell.
Produced by `experiments/10_ladder_ab_campaign.py`.

The CSV reproduces the column order of `runs_flat_v2.csv` in the campaign
repository exactly, then appends the extra fields, so the verification scripts
read it unmodified.

## `angle_convention/`

`D1_ang_runs.jsonl` and `D1_ang_runs.csv`, same schema plus one column,
`undo_angle_deg`. 100 rows, 25 per (convention, ladder) cell, acquired
7 October 2026, 11:03–14:31 UTC, interleaved run by run and paired by
repetition. Produced by `experiments/11_angle_convention_campaign.py`.

The two arms command the two level-1 inverse rotations as `-56.7891` and
`-33.2109` degrees (`direct`) or as `663.2109` and `686.7891` (`folded`). These
are the same gates: a 720° rotation of a spin-½ is the identity, phase
included. Everything else is common to all four cells.

### What the instrument's record confirms

| check | result |
|---|---|
| runs in our record found in the export | 100 of 100 |
| measured distributions agree | max deviation 6.1×10⁻¹⁶ |
| `folded` vs `direct`, differing gate positions | exactly 2, both the level-1 undo rotations |
| angle the device was actually commanded, `folded` | 663.211° and 686.789° — **not wrapped** |
| `direct` vs the 5 October `control` circuit | identical on all 21 gates, both ladders |
| job status | SUCCEED on every campaign job |

The fourth row is the one that makes the experiment interpretable: the device
was asked for a 663° rotation and performed it as such, so the null reported in
the paper is a statement about the hardware and not about an internal
angle-wrapping convention.

The fifth row is what makes the 5→7 October comparison a circuit-exact repeat,
and therefore a longitudinal measurement with nothing varying but the date.

## `spinqit_export_angle_20261007.csv`

1239 experiments spanning 30 January to 7 October 2026 — the full job database
at that date, so it supersedes `spinqit_export_ladderab_20261005.csv` and
contains it. Both are kept: the October 5 analysis in the paper was run against
the earlier export and should remain reproducible against the file it used.

## `pilot_d1_20260308/`

The earlier D1-only session at 2048 shots with a full 8×8 confusion matrix:
run logs raw, mitigated and partially mitigated, the matrix itself, and the
summary tables. This is the readout-mitigation experiment the paper reports as
its reason for using raw distributions throughout.

## Seven orphaned jobs in the 5 October export — not part of the dataset

The export contains 107 jobs carrying 5 October campaign run names, against 100
recorded runs. Seven names appear twice:

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
different backend, so this cannot recur. The 7 October campaign has no such
duplicates: all 100 of its job names occur exactly once.

## Bit-order calibration

Jobs `BITORDER_X_q0/q1/q2`, executed 18:15–18:19 UTC on 5 October, fifteen
minutes before that campaign, by `calibrate_bit_order.py`. Dominant outcomes
`100`, `010`, `001` respectively, confirming the canonical MSB→LSB convention
(`SPINQ_BITORDER=MSB->LSB`) on this unit with instrument-side evidence.
Dominant probabilities 0.657, 0.805 and 0.864 — the readout is visibly
imperfect on q0, but the assignment is unambiguous.

## Reproducing the paper's numbers

Every numerical table and every number quoted in running text in the manuscript
and its supplement is generated from these files by
`gr-triangulum-verification/scripts/90_make_tables.py`:

    python scripts/90_make_tables.py \
        --outdir <overleaf folder>/tables \
        --angle-export ../grover-rudolph-practical-implementation/data/spinqit_export_angle_20261007.csv

Nothing in either document is transcribed by hand.
