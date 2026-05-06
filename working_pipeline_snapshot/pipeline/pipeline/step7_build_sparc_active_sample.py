"""
step7_build_sparc_active_sample.py
====================================
Assembles the final per-galaxy summary table for the 38 active SPARC galaxies.

Inputs (all in results/chains/):
  delta_chi2_per_galaxy.csv       -- Δχ² and sign predictions   (step4)
  loo_predictions.csv             -- LOO Areq predictions        (step2)
  hybrid_loo_per_galaxy.csv       -- continuous Areq_pred_physical
  sparc_all_clean_inner_shape_proxy.csv -- c200/cLCDM_DM14 ratios
  mouter_per_galaxy.csv           -- outer residuals              (step1)

Output (in results/chains/):
  sparc_active_sample.csv

Key columns:
  name, dR1, sinner, mouter, Areq, sign_correct_physical,
  dchi2_physical, dchi2_oracle, Areq_pred_physical, Areq_pred_oracle,
  c200_lcdm_ratio

Paper II headline results reproduced here:
  LOO correct     : 29 / 38
  Δχ²_physical    : -4108
  Δχ²_oracle      : -4614

Usage (from repo root):
    python pipeline/step7_build_sparc_active_sample.py
"""

import os
import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT  = os.path.dirname(SCRIPT_DIR)
CHAINS     = os.path.join(REPO_ROOT, "results", "chains")

# ── Inputs ────────────────────────────────────────────────────────────────────
dchi   = pd.read_csv(os.path.join(CHAINS, "delta_chi2_per_galaxy.csv"))
loo    = pd.read_csv(os.path.join(CHAINS, "loo_predictions.csv"))
hybrid = pd.read_csv(os.path.join(CHAINS, "hybrid_loo_per_galaxy.csv"))
sparc  = pd.read_csv(os.path.join(CHAINS, "sparc_all_clean_inner_shape_proxy.csv"))

# ── c200/cLCDM ratio (Dutton & Maccio 2014 reference) ────────────────────────
h = 0.7035
sparc["log10c_LCDM"] = 0.905 - 0.101 * (sparc["logM200"] + np.log10(h) - 12)
sparc["c200_lcdm_ratio"] = 10**sparc["log10c_ref"] / 10**sparc["log10c_LCDM"]

# ── Assemble ──────────────────────────────────────────────────────────────────
df = dchi[["name", "dR1", "Areq", "sign_pred", "sign_true",
           "correct", "dchi2_rule", "dchi2_oracle"]].rename(columns={
    "dchi2_rule":  "dchi2_physical",
    "sign_pred":   "Areq_pred_sign",
    "correct":     "sign_correct_physical",
})

# sinner, mouter from LOO table
df = df.merge(loo[["name", "sinner", "mouter"]], on="name", how="left")

# Continuous Areq_pred_physical from hybrid LOO
df = df.merge(
    hybrid[["name", "A_pred_hybrid_loo"]].rename(
        columns={"A_pred_hybrid_loo": "Areq_pred_physical"}),
    on="name", how="left")

# Oracle = true Areq (upper bound with perfect sign)
df["Areq_pred_oracle"] = df["Areq"]

# c200/cLCDM
df = df.merge(sparc[["name", "c200_lcdm_ratio"]], on="name", how="left")

# ── Validate ──────────────────────────────────────────────────────────────────
n_correct = int(df["sign_correct_physical"].sum())
dchi2_tot = df["dchi2_physical"].sum()
dchi2_orc = df["dchi2_oracle"].sum()

assert len(df) == 38,          f"Expected 38 rows, got {len(df)}"
assert n_correct == 29,        f"Expected 29/38 correct, got {n_correct}/38"
assert abs(dchi2_tot + 4108) < 2, f"Expected Δχ²≈-4108, got {dchi2_tot:.1f}"
assert df.isnull().sum().sum() == 0, "NaN values detected"

print(f"LOO correct     : {n_correct}/38")
print(f"Δχ²_physical    : {dchi2_tot:.1f}")
print(f"Δχ²_oracle      : {dchi2_orc:.1f}")

# ── Save ──────────────────────────────────────────────────────────────────────
OUT = os.path.join(CHAINS, "sparc_active_sample.csv")
df.to_csv(OUT, index=False)
print(f"Saved: {OUT}")
