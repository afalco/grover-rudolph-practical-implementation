# Running the unit tests

This repository includes optional unit tests under `tests/` to validate the **math layer**
(angles, UCRy linear system, metrics) independently of any hardware backend.

These tests are **simulation-only** and do **not** require NMR access or `SPINQ_*` environment
variables.

---

## 1) Install test dependencies

From the repo root, make sure your environment is active:

```bash
conda activate spinq-gr
```

Install `pytest` (and your normal requirements):

```bash
pip install -r requirements.txt
pip install pytest
```

---

## 2) Run all tests

From the repo root:

```bash
pytest -q
```

---

## 3) Run a single test file

Examples:

```bash
pytest -q tests/test_angles.py
pytest -q tests/test_ucry.py
pytest -q tests/test_metrics.py
```

---

## 4) Run a single test function

Example:

```bash
pytest -q tests/test_ucry.py::test_ucry_coefficients_reconstruct_thetas
```

---

## Troubleshooting

### `ModuleNotFoundError: No module named 'gr'`

Run `pytest` from the **repository root**, and confirm `gr/` is a Python package:

- `gr/__init__.py` must exist.

If you still have issues, run:

```bash
python -c "import gr; print('gr imported OK')"
```

### SpinQit not available on CI / local env

The unit tests should only import the pure-python pieces of `gr/`.
If your `gr/__init__.py` imports SpinQit unconditionally, consider changing it so that
SpinQit-dependent modules are imported lazily (only when needed for circuit/backend runs).
