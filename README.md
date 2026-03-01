# Grover–Rudolph Practical Implementation

This repository provides a practical, reproducible implementation of the Grover–Rudolph state-preparation method for 3 qubits, including:

- **Ideal simulation** verification (sanity check vs target distribution)
- **Execution on a real quantum device** (SpinQ Triangulum NMR)
- **Staged diagnostics** (L0 / L01 / FULL)
- **Metrics**: Total Variation (TV), L2 distance, classical fidelity
- Optional: **8×8 readout calibration and mitigation**

## Reference paper

The implementation is designed as a “paper-to-hardware” companion for:

- `paper/Grover_Rudolph_proof.pdf`

## What the code demonstrates

1. **Correctness (ideal)**  
   The circuit prepares the target distribution exactly in an ideal simulator (TV≈0, Fidelity≈1).

2. **Hardware compilation (ancilla-free)**  
   The implementation avoids nested multi-controlled parametric gates by decomposing the level-2 uniformly-controlled rotation using **Ry + CX** only.

3. **Real-device behavior**  
   We report raw and (optionally) readout-mitigated results on SpinQ Triangulum NMR, including per-stage performance.

## Repository layout

- `src/definitive_gr.py`: main end-to-end runner (SIM + NMR + metrics)
- `src/gr/`: reusable library modules (angles, UCRy, metrics, backends, readout)
- `experiments/`: focused scripts for specific validations (bit order, ladder tests, staged runs)
- `results/`: saved logs/CSVs/figures from runs
- `docs/`: technical notes and LaTeX documentation

## Installation

### Option A: Conda
```bash
conda env create -f environment.yml
conda activate spinq
