# experiments/

This folder contains **reproducible experiment scripts** built on top of the `gr/` package.  
Each script is intentionally small and single-purpose so you can isolate effects and keep a clean
“paper → code → hardware” trace.

---

## Quick start

From the repo root:

```bash
python experiments/00_smoke_sim_only.py

Hardware setup (SpinQ Triangulum NMR)
```
Hardware runs require environment variables:
	•	SPINQ_IP
	•	SPINQ_PORT
	•	SPINQ_USER
	•	SPINQ_PASS

