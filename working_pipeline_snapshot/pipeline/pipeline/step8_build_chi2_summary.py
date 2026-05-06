"""
step8_build_chi2_summary.py
============================
Builds the chi2_summary.csv used by figI4 bar chart (Paper I Fig. 4d).

Inputs (all in results/chains/):
  delta_chi2_per_galaxy.csv       -- dchi2_physical, dchi2_oracle  (step4)
  type95_per_galaxy.csv           -- Type-9.5 baseline             (step4b)
  regime2_complete.csv            -- Regime 2 c200 total           (step9)

Output (in results/chains/):
  chi2_summary.csv

Columns: rule, minus_dchi2, label, bar_label, n_galaxies

Bar values:
  Baseline (Type-9.5) : |sum dchi2_type95|   from type95_per_galaxy.csv
  V200-only model     : |sum dchi2_oracle|    from delta_chi2_per_galaxy.csv
  Full BPB model      : oracle_V200 + |regime2 total|

Usage (from repo root):
    python pipeline/step8_build_chi2_summary.py
"""

import os
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT  = os.path.dirname(SCRIPT_DIR)
CHAINS     = os.path.join(REPO_ROOT, "results", "chains")

# ── Inputs ────────────────────────────────────────────────────────────────────
dchi  = pd.read_csv(os.path.join(CHAINS, "delta_chi2_per_galaxy.csv"))
r2    = pd.read_csv(os.path.join(CHAINS, "regime2_complete.csv"))

T95_PATH = os.path.join(CHAINS, "type95_per_galaxy.csv")
# type95_per_galaxy.csv généré par step4b — non requis si Option A (2 barres)

# ── Compute totals from data ──────────────────────────────────────────────────
oracle_v200 = abs(int(round(dchi["dchi2_oracle"].sum())))
regime2_val = abs(int(round(r2["delta_chi2"].sum())))
full_bpb    = oracle_v200 + regime2_val

print(f"V200-only model     : {oracle_v200}")
print(f"Regime 2 (c200)     : {regime2_val}")
print(f"Full BPB model      : {full_bpb}")

# ── Build table ───────────────────────────────────────────────────────────────
summary = pd.DataFrame([
    {"rule": "physical_v200",
     "minus_dchi2": oracle_v200,
     "label": "V200-only\nmodel",
     "n_galaxies": 38},
    {"rule": "full_bpb",
     "minus_dchi2": full_bpb,
     "label": "Full BPB\nmodel",
     "n_galaxies": 118},
])

# ── Save ──────────────────────────────────────────────────────────────────────
OUT = os.path.join(CHAINS, "chi2_summary.csv")
summary.to_csv(OUT, index=False)
print(f"Saved: {OUT}")
