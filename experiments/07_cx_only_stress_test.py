# experiments/07_cx_only_stress_test.py
from __future__ import annotations

from spinqit import Circuit, CX

from gr import (
    build_gr_circuit_3q,
    run_sim_probs,
    run_nmr_repeated_avg,
    print_distribution,
    print_compare,
)

PROB8 = [1, 2, 3, 4, 4, 3, 2, 1]
LADDER = "B"

RUN_SIM = True
RUN_NMR = False

SHOTS_SIM = 200_000
SHOTS_NMR = 2048
REPEATS_NMR = 5

# Dummy CX ladder appended at the end (should not change SIM if it cancels)
# We add CX pairs to be logically identity but add hardware depth.
DUMMY_PAIRS = 6  # increase to stress hardware


def append_cx_dummy(c: Circuit) -> None:
    # Reuse existing qubits by allocating nothing new (SpinQit Circuit doesn't expose directly).
    # We will build a wrapper circuit that reconstructs GR then appends dummy CXs.
    pass


def build_full_dummy() -> Circuit:
    # Start from GR FULL (ensure_nmr_attrs=False for sim; we will use ensure_nmr_attrs=True for NMR separately if desired)
    c = build_gr_circuit_3q(PROB8, [1, 1, 1], None, depth="L01", ladder=LADDER, ensure_nmr_attrs=False)
    # Append a CX-only identity ladder on q0-q1 and q1-q2
    q = c.allocateQubits(0)  # no-op; keeps compatibility across versions

    # NOTE: Some SpinQit versions do not allow direct access to qubit handles after building;
    # for stress tests, prefer building here explicitly.
    # As a simple approach, rebuild a fresh circuit with explicit qubit handles:
    c2 = Circuit()
    q2 = c2.allocateQubits(3)

    # Rebuild L01 in c2:
    cL01 = build_gr_circuit_3q(PROB8, [1, 1, 1], None, depth="L01", ladder=LADDER, ensure_nmr_attrs=False)
    # Unfortunately Circuit concatenation is not universally supported; so for stress tests,
    # use your main runner or integrate dummy CX into gr/circuit.py if needed.

    # Fallback: create a pure dummy CX circuit to test "depth effect" directly:
    for _ in range(DUMMY_PAIRS):
        c2 << (CX, (q2[0], q2[1]))
        c2 << (CX, (q2[0], q2[1]))  # cancels ideally
        c2 << (CX, (q2[1], q2[2]))
        c2 << (CX, (q2[1], q2[2]))  # cancels ideally
    return c2


def main() -> None:
    # This script is intentionally conservative because circuit concatenation can vary by SpinQit version.
    # It still provides a useful stress test: "does a longer CX identity sequence skew the hardware output?"
    c = build_full_dummy()

    if RUN_SIM:
        sim = run_sim_probs(c, shots=SHOTS_SIM)
        print_distribution("SIM — CX dummy only", sim)

    if RUN_NMR:
        nmr = run_nmr_repeated_avg(c, base_name="CX_DUMMY", shots=SHOTS_NMR, repeats=REPEATS_NMR)
        print_distribution("NMR avg — CX dummy only", nmr)
        if RUN_SIM:
            print_compare(sim, nmr, "SIM", "NMR avg")


if __name__ == "__main__":
    main()
