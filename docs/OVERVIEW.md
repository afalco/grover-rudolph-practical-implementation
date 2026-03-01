# docs/OVERVIEW.md

## Purpose of this repository

This repository provides a **practical, hardware-facing implementation** of the Grover–Rudolph state preparation algorithm for **3 qubits**, with:

- an **ideal simulation** reference (sanity check against the target distribution),
- an **ancilla-free transpilation** into the gate set **{Ry, X, CNOT}**,
- **execution on a real NMR quantum device** (SpinQ Triangulum),
- **staged diagnostics** (L0 / L01 / FULL) to locate where hardware deviations appear,
- quantitative **metrics**: Total Variation (TV), L2 distance, and classical fidelity.

The implementation is designed as a “paper-to-hardware” companion to:

**A Rigorous and Self–Contained Proof of the Grover–Rudolph State Preparation Algorithm**  
arXiv:2601.17930 — https://arxiv.org/abs/2601.17930

## What the algorithm does (3-qubit case)

Given a target probability distribution \(p(x)\) over bitstrings \(x \in \{0,1\}^3\), Grover–Rudolph prepares a quantum state

\[
|\psi\rangle = \sum_{x \in \{0,1\}^3} \sqrt{p(x)}\,|x\rangle,
\]

so that a computational basis measurement returns \(x\) with probability \(p(x)\).

In this repository, the target distribution is specified as a length-8 vector `PROB8` in the **MSB→LSB** basis order:

\[
|q_0 q_1 q_2\rangle \equiv |000\rangle, |001\rangle, \ldots, |111\rangle.
\]

## Key design choices

### 1) Angle convention: asin / “child=1”
A common source of errors in Grover–Rudolph implementations is mixing conventions such as `asin` vs `acos` and/or swapping which child of the dyadic tree corresponds to \(|0\rangle\) or \(|1\rangle\).

This repo uses a consistent convention:

\[
\theta(\text{prefix}) = 2\,\arcsin\sqrt{\frac{P(\text{prefix}+1)}{P(\text{prefix})}},
\]

so that applying \(R_y(\theta)\) to \(|0\rangle\) yields

\[
\Pr(1 \mid \text{prefix}) = \sin^2(\theta/2).
\]

This convention matches the staged circuit structure and enables an exact simulator sanity check.

### 2) Hardware-friendly compilation: only {Ry, X, CNOT}
Instead of relying on nested parametric controlled gates (which can compile poorly on constrained backends), all controlled rotations are decomposed into **Ry + CNOT (+ X for 0-controls)**:

- Level 1 (single control): CRy decomposition using two CNOTs.
- Level 2 (two controls): a **2-control uniformly controlled Ry (UCRy)** built as a **Ry–CNOT ladder**.

### 3) Staged diagnostics (L0 / L01 / FULL)
The code supports three depths:

- **L0**: only the top rotation on \(q_0\),
- **L01**: adds level-1 rotations on \(q_1\) conditioned on \(q_0\),
- **FULL**: adds the level-2 UCRy on \(q_2\) conditioned on \(q_0,q_1\).

This staging helps determine whether deviations arise from:
- state preparation / readout bias (already visible at L0),
- single-control operations (degradation at L01),
- depth/crosstalk in the final UCRy layer (extra degradation at FULL).

## How to run

### Simulation-only
Set `RUN_NMR = False` and run:

```bash
python src/definitive_gr.py
